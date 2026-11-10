"""Run with: python -m unittest python.preprocessing.test_local_prepare."""
import contextlib
import io
import random
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pgenlib

from .prepare import GeneSelectionRequest, PrepareRequest, prepare_chromosomes
from .site import read_public_var_list, write_public_var_list

ANALYSIS = Path(__file__).resolve().parents[1] / 'analysis'
sys.path.insert(0, str(ANALYSIS / 'python_skat'))
from load_preprocessed import global_maf_keep, read_gene_genotype, read_plain_input
from compute_reference import analyze


class LocalPrepareTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / 'chr21'
        self.keys = [f'21:{i}:A:G' for i in range(1, 6)]
        self.dosages = np.zeros((5, 100), dtype=np.int8)
        self.dosages[0, 0] = 1
        self.dosages[1, :2] = 1
        self.dosages[2, :3] = 1
        self.dosages[3, :] = 2
        self.dosages[3, 0] = 0
        self.dosages[4, 0] = 1
        with pgenlib.PgenWriter(str(self.prefix).encode() + b'.pgen', 100, variant_ct=5) as writer:
            for row in self.dosages:
                writer.append_biallelic(row)
        Path(f'{self.prefix}.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\tFILTER\n' + ''.join(
            f'21\t{i+1}\t{key}\tA\tG\t' + ('FAIL' if i == 4 else 'PASS') + '\n'
            for i, key in enumerate(self.keys)))
        Path(f'{self.prefix}.psam').write_text('#IID\n' + ''.join(f's{i}\n' for i in range(100)))
        (self.root / 'genes.tsv').write_text('gene_id\tgene_symbol\tchromosome\torder_index\n' + ''.join(
            f'g{i}\tG{i}\t21\t{i}\n' for i in range(5)))
        (self.root / 'annotation.tsv').write_text('variant_key\tgene_id\tgene_symbol\n' + ''.join(
            f'{key}\tg{i}\tG{i}\n' for i, key in enumerate(self.keys)))
        (self.root / 'pheno.tsv').write_text('IID\ty\n' + ''.join(
            f's{i}\t{np.sin(i) + i / 100}\n' for i in range(100)))
        (self.root / 'cov.tsv').write_text('IID\tPC1\n' + ''.join(f's{i}\t{i}\n' for i in range(100)))
        (self.root / 'anc.tsv').write_text('IID\tanc\n' + ''.join(
            f's{i}\t' + ('EUR' if i < 60 else 'AFR') + '\n' for i in range(100)))
        self.request = PrepareRequest(
            run_dir=self.root / 'a', chromosomes=(21,), genotype=str(self.prefix),
            gene_panel=str(self.root / 'genes.tsv'), annotation=str(self.root / 'annotation.tsv'),
            phenotype=self.root / 'pheno.tsv', covariate=self.root / 'cov.tsv', ancestry=self.root / 'anc.tsv',
            plink2_bin='unused', phenotype_id_column='IID', covariate_id_column='IID',
            covariate_column='', covariate_columns=('PC1',), is_cov_single_column=False,
            ancestry_id_column='IID', ancestry_column='anc', phenotype_columns=('y',),
            ancestries=('EUR', 'ALL'), local_max_maf=0.005, samples_per_cohort=0,
            shared_rate=0.6, gene_selection=GeneSelectionRequest('all', 0, None), party=1,
            mode='both', global_max_maf=(0.01, 0.005), prepared_cache_dir=self.root / 'cache-a',
        )

    def prepare(self, request):
        seeded = random.Random(42)
        with contextlib.redirect_stdout(io.StringIO()), patch('python.preprocessing.pipeline.random.Random', return_value=seeded):
            prepare_chromosomes(request)
        return request.run_dir / 'prepared' / 'local'

    def test_modes_and_independent_selection(self):
        local = self.prepare(self.request)
        joint = self.request.run_dir / 'prepared' / 'ALL' / 'chr21'
        self.assertEqual((joint / 'block_sizes.txt').read_text().split(), ['1', '0', '0', '0', '0'])
        for anc in self.request.ancestries:
            directory = local / anc / 'chr21'
            self.assertEqual((directory / 'block_sizes.txt').read_text().split(), ['1', '1', '0', '1', '0'])
            self.assertEqual((directory / 'A/cohort_size.txt').read_text().strip(), '100')
            self.assertEqual((directory / 'A/allele_counts.txt').read_text().splitlines(), ['1', '2', '', '198', ''])
            data = read_plain_input(directory, local=True)
            self.assertEqual([x.tolist() for x in global_maf_keep(data, .01).keep_public],
                             [[True], [True], [], [True], []])
            self.assertEqual([x.tolist() for x in global_maf_keep(data, .005).keep_public],
                             [[True], [False], [], [False], []])
        # B's standalone run needs neither A's gene list nor its cut-list, even if a stale path is configured.
        b = replace(self.request, run_dir=self.root / 'b', prepared_cache_dir=self.root / 'cache-b',
                    party=2, mode='local', public_var_list=self.root / 'missing.tsv', local_max_maf=.001)
        b_local = self.prepare(b)
        self.assertFalse((b.run_dir / 'public_variants.tsv').exists())
        self.assertFalse((b.run_dir / 'prepared/ALL').exists())
        for name in ('genes.txt', 'block_sizes.txt', 'A/allele_counts.txt', 'A/geno/block.1.bin'):
            self.assertEqual((local / 'ALL/chr21' / name).read_bytes(), (b_local / 'ALL/chr21' / name).read_bytes())
        # In both mode, B's joint input obeys A's cut-list; its local input does not.
        public_dir = self.root / 'public'
        public_dir.mkdir()
        write_public_var_list(public_dir / 'public_variants.tsv',
                              read_public_var_list(self.request.run_dir / 'public_variants.tsv')[:1])
        (public_dir / 'local_cut_variants.tsv').write_bytes((self.request.run_dir / 'local_cut_variants.tsv').read_bytes())
        b_both = replace(b, run_dir=self.root / 'b-both', prepared_cache_dir=self.root / 'cache-b-both',
                         mode='both', public_var_list=public_dir / 'public_variants.tsv')
        both_local = self.prepare(b_both)
        self.assertEqual((both_local / 'ALL/chr21/block_sizes.txt').read_bytes(),
                         (local / 'ALL/chr21/block_sizes.txt').read_bytes())
        joint_only = replace(self.request, run_dir=self.root / 'joint', prepared_cache_dir=self.root / 'cache-joint', mode='joint')
        self.prepare(joint_only)
        self.assertFalse((joint_only.run_dir / 'prepared/local').exists())
        self.assertEqual((joint_only.run_dir / 'prepared/ALL/chr21/A/geno/block.0.bin').read_bytes(),
                         (joint / 'A/geno/block.0.bin').read_bytes())

    def test_same_data_local_and_joint_plain(self):
        local = self.prepare(self.request) / 'ALL/chr21'
        data = read_plain_input(local, local=True)
        joint = self.root / 'plain-joint'
        (joint / 'A/geno').mkdir(parents=True)
        (joint / 'B/geno').mkdir(parents=True)
        (joint / 'B/private').mkdir()
        for name in ('genes.txt', 'block_sizes.txt', 'public_variants.tsv'):
            (joint / name).write_bytes((local / name).read_bytes())
        for side, rows in (('A', slice(0, 50)), ('B', slice(50, 100))):
            np.savetxt(joint / side / 'cov.txt', data.covariates[rows])
            np.savetxt(joint / side / 'pheno.txt', data.phenotypes[rows])
            counts = []
            for i in range(len(data.genes)):
                geno = read_gene_genotype(data, i)[rows].astype(np.int8)
                geno.tofile(joint / side / 'geno' / f'block.{i}.bin')
                counts.append(' '.join(map(str, geno.sum(axis=0))))
                if side == 'B':
                    (joint / 'B/private' / f'block.{i}.bin').touch()
            (joint / side / 'allele_counts.txt').write_text('\n'.join(counts) + '\n')
            (joint / side / 'cohort_size.txt').write_text('50\n')
        (joint / 'B/private_allele_counts.txt').write_text('\n' * len(data.genes))
        for threshold in (.01, .005):
            standalone = global_maf_keep(data, threshold)
            combined = global_maf_keep(read_plain_input(joint), threshold)
            for i in range(len(data.genes)):
                np.testing.assert_array_equal(read_gene_genotype(standalone, i), read_gene_genotype(combined, i))
            local_rows, joint_rows = analyze(standalone), analyze(combined)
            self.assertEqual(local_rows, joint_rows)
            self.assertEqual(local_rows[2]['burden_p'], 1)


if __name__ == '__main__':
    unittest.main()
