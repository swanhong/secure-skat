import csv
import os
import subprocess
import unittest
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event
from unittest.mock import Mock, patch

import numpy as np
import pgenlib

from .plink import pgen_extract
from .prepare import (
    GeneSelectionRequest, PrepareRequest, prepare_cached_blocks,
    prepare_chromosomes, prepared_cache_config,
)


def plink_reference(prefix, sample_ids, variant_keys):
    """Independent PLINK CLI oracle; never used by production preparation."""
    if not variant_keys:
        return np.zeros((len(sample_ids), 0), dtype=np.int8), ()
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        (root / 'keep').write_text('#IID\n' + '\n'.join(sample_ids))
        (root / 'extract').write_text('\n'.join(variant_keys))
        (root / 'alleles').write_text(''.join(f'{key}\t{key.rsplit(":", 1)[-1]}\n' for key in variant_keys))
        subprocess.run([
            os.environ.get('PLINK2', 'plink2'), '--pfile', str(prefix),
            '--keep', str(root / 'keep'), '--extract', str(root / 'extract'),
            '--export', 'A', 'include-alt', '--export-allele', str(root / 'alleles'),
            '--max-alleles', '2', '--threads', '1', '--memory', '640', '--out', str(root / 'geno'),
        ], check=True, capture_output=True)
        with (root / 'geno.raw').open() as file:
            reader = csv.reader(file, delimiter='\t')
            header = next(reader)
            start = header.index('PHENOTYPE') + 1
            keys = tuple(column.split('(', 1)[0].rpartition('_')[0] for column in header[start:])
            rows = {fields[header.index('IID')]: fields[start:] for fields in reader}
        values = [[{'0': 0, '1': 1, '2': 2, 'NA': 0}[value] for value in rows[iid]] for iid in sample_ids]
        return np.array(values, dtype=np.int8), keys


def block_files(directory):
    return {str(path.relative_to(directory)): path.read_bytes()
            for path in directory.rglob('*') if path.is_file() and path.name != 'config.json'}


class PreparedTest(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = self.root = Path(self.temporary.name)
        self.prefix = root / 'chr1'
        self.ids = tuple(f's{index}' for index in range(6))
        self.keys = tuple(f'1:{index + 1}:A:G' for index in range(258)) + ('1:259:A:<INS:ME:ALU>',) + ('1:260:A:C,G',)
        self.genotypes = np.array([[(i + j) % 3 for j in range(6)] for i in range(259)], dtype=np.int8)
        self.genotypes[0, 0] = -9
        with pgenlib.PgenWriter(bytes(self.prefix) + b'.pgen', 6, variant_ct=260, allele_ct_limit=3) as writer:
            writer.append_biallelic_batch(self.genotypes)
            writer.append_alleles(np.array([0, 2] * 6, dtype=np.int32), allele_ct=3)
        Path(f'{self.prefix}.psam').write_text('#IID\n' + '\n'.join(self.ids) + '\n')
        Path(f'{self.prefix}.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
            f'1\t{i + 1}\t{key}\tA\t{key.split(":", 3)[-1]}\n' for i, key in enumerate(self.keys)))
        (root / 'panel1.tsv').write_text('gene_id\tgene_symbol\tchromosome\torder_index\nG1\tG1\t1\t0\nG2\tG2\t1\t1\n')
        (root / 'annotation1.tsv').write_text('variant_key\tgene_id\tgene_symbol\tLoF\tMAF\n' + ''.join(
            f'{key}\tG1\tG1\tHC\t0.001\n' for key in self.keys[:130]) + ''.join(
            f'{key}\tG2\tG2\tHC\t0.001\n' for key in self.keys[120:]))
        (root / 'phenotype.csv').write_text('IID,y,z\n' + ''.join(f'{iid},{i},{i + 2}\n' for i, iid in enumerate(self.ids)))
        (root / 'ancestry.tsv').write_text('IID\tpcs\tancestry\n' + ''.join(f'{iid}\t[{i},1]\tEUR\n' for i, iid in enumerate(self.ids)))
        self.request = PrepareRequest(
            run_dir=root / 'run', chromosomes=(1,), genotype=str(root / 'chr{chromosome}'),
            gene_panel=str(root / 'panel{chromosome}.tsv'), annotation=str(root / 'annotation{chromosome}.tsv'),
            phenotype=root / 'phenotype.csv', covariate=root / 'ancestry.tsv', ancestry=root / 'ancestry.tsv',
            plink2_bin='plink2', phenotype_id_column='IID', covariate_id_column='IID', covariate_column='pcs',
            ancestry_id_column='IID', ancestry_column='ancestry', phenotype_columns=('y',), ancestries=('EUR',),
            num_cov=1, mask={'LoF': 'HC'}, max_maf=0.02, samples_per_cohort=0, sample_seed=42,
            role_seed=42, shared_rate=0.6, gene_selection=GeneSelectionRequest('all', 0, 42, None),
            prepared_cache_dir=root / 'cache',
        )

    def test_direct_reader_matches_plink_and_expected_values(self):
        ids = ('s5', 's0', 's3')
        wanted = self.keys[::-1] + ('1:999:A:G',)
        with patch('subprocess.run', side_effect=AssertionError('PLINK called')):
            actual, keys = pgen_extract(self.prefix, ids, wanted)
        expected, expected_keys = plink_reference(self.prefix, ids, wanted)
        self.assertEqual(keys, expected_keys)
        self.assertEqual(keys, self.keys[:-1])
        np.testing.assert_array_equal(actual, expected)
        truth = self.genotypes[:, [5, 0, 3]].T.copy()
        truth[truth == -9] = 0
        truth[:, -1] = 0  # Historical PLINK export for colon-containing symbolic ALT.
        np.testing.assert_array_equal(actual, truth)
        self.assertEqual(actual.dtype, np.int8)
        self.assertTrue(actual.flags.c_contiguous)
        for wanted in ((), ('absent',), (self.keys[-1],)):
            matrix, emitted = pgen_extract(self.prefix, ids, wanted)
            self.assertEqual(matrix.shape, (3, 0))
            self.assertEqual(emitted, ())

    def test_fractional_dosage_and_bad_alignment_are_rejected(self):
        with pgenlib.PgenWriter(bytes(self.prefix) + b'.pgen', 6, variant_ct=1, dosage_present=True) as writer:
            writer.append_dosages(np.array([0.25, 0, 1, 2, -9, 1], dtype=np.float32))
        Path(f'{self.prefix}.pvar').write_text(f'#CHROM\tPOS\tID\tREF\tALT\n1\t1\t{self.keys[0]}\tA\tG\n')
        with self.assertRaisesRegex(ValueError, 'fractional'):
            pgen_extract(self.prefix, self.ids, self.keys[:1])
        with self.assertRaisesRegex(ValueError, 'sample set mismatch'):
            pgen_extract(self.prefix, ('unknown',), self.keys[:1])
        Path(f'{self.prefix}.pvar').write_text(f'#CHROM\tPOS\tID\tREF\tALT\n1\t1\t{self.keys[0]}\tA\tC\n')
        with self.assertRaisesRegex(ValueError, 'ALT'):
            pgen_extract(self.prefix, self.ids, self.keys[:1])

    def test_prepared_bytes_match_plink_and_cache_is_reused(self):
        for rate in (0.6, 1.0):
            request = replace(self.request, shared_rate=rate)
            prepare_chromosomes(request)
            link = request.run_dir / 'prepared/EUR/chr1'
            self.assertTrue(link.is_symlink())
            cached = link.resolve()
            baseline = replace(request, run_dir=self.root / f'baseline{rate}', prepared_cache_dir=None)
            prepare_chromosomes(baseline, extractor=plink_reference)
            self.assertEqual(block_files(cached), block_files(baseline.run_dir / 'prepared/EUR/chr1'))
            with patch('rewrite.preprocessing.prepare.pgen_extract', side_effect=AssertionError('cache miss')):
                other = replace(request, run_dir=self.root / f'other{rate}')
                prepare_chromosomes(other)
            self.assertEqual((other.run_dir / 'prepared/EUR/chr1').resolve(), cached)
            self.assertTrue((cached / 'config.json').is_file())
            if rate == 1:
                self.assertTrue(all(path.stat().st_size == 0 for path in (cached / 'B/private').iterdir()))

    def test_cache_identity_tracks_inputs_and_ignores_execution_scope(self):
        config = prepared_cache_config(self.request, 1, 'EUR')
        equivalent = replace(self.request, run_dir=self.root / 'elsewhere', prepared_cache_dir=self.root / 'other-cache',
                             chromosomes=(1, 2), ancestries=('EUR', 'AFR'), plink2_bin='another-plink', mask={'LoF': ['HC']})
        self.assertEqual(config, prepared_cache_config(equivalent, 1, 'EUR'))
        for field, value in (('max_maf', 0.01), ('sample_seed', 1), ('role_seed', 1), ('shared_rate', 1),
                             ('samples_per_cohort', 2), ('phenotype_columns', ('z',)), ('num_cov', 2),
                             ('mask', {'LoF': 'LC'}), ('gene_selection', GeneSelectionRequest('random', 1, 7, None))):
            with self.subTest(field=field):
                self.assertNotEqual(config, prepared_cache_config(replace(self.request, **{field: value}), 1, 'EUR'))
        for path in (Path(f'{self.prefix}.pgen'), Path(f'{self.prefix}.pvar'), Path(f'{self.prefix}.psam'),
                     self.request.phenotype, self.request.covariate, self.root / 'panel1.tsv', self.root / 'annotation1.tsv'):
            before = prepared_cache_config(self.request, 1, 'EUR')
            stat = path.stat()
            os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000000))
            self.assertNotEqual(before, prepared_cache_config(self.request, 1, 'EUR'))
        request = replace(self.request, gene_selection=GeneSelectionRequest('file', 0, 42, self.root / 'panel1.tsv'))
        self.assertIn('mtime_ns', prepared_cache_config(request, 1, 'EUR')['gene_selection']['path'])

    def test_failed_build_is_not_published_and_real_outputs_are_preserved(self):
        def fail(out_dir):
            (out_dir / 'pos.txt').touch()
            raise RuntimeError('interrupted')
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            prepare_cached_blocks(self.request, 1, 'EUR', fail)
        self.assertFalse(list(self.request.prepared_cache_dir.rglob('config.json')))
        self.assertFalse((self.request.run_dir / 'prepared/EUR/chr1').exists())
        prepare_chromosomes(self.request)
        link = self.request.run_dir / 'prepared/EUR/chr1'
        cached = link.resolve()
        link.unlink()
        link.mkdir()
        (link / 'keep').write_text('keep')
        with self.assertRaisesRegex(ValueError, 'prepare --clear'):
            prepare_chromosomes(self.request)
        self.assertEqual((link / 'keep').read_text(), 'keep')
        self.assertTrue((cached / 'config.json').is_file())
        with self.assertRaisesRegex(ValueError, 'outside run_dir'):
            prepare_chromosomes(replace(self.request, prepared_cache_dir=self.request.run_dir / 'cache'))

    def test_concurrent_requests_build_once(self):
        entered, release = Event(), Event()
        def build(out_dir):
            entered.set()
            if not release.wait(5):
                raise TimeoutError('test did not release builder')
            (out_dir / 'pos.txt').touch()
        build = Mock(side_effect=build)
        with ThreadPoolExecutor(2) as pool:
            first = pool.submit(prepare_cached_blocks, self.request, 1, 'EUR', build)
            self.assertTrue(entered.wait(5))
            second = pool.submit(prepare_cached_blocks, replace(self.request, run_dir=self.root / 'run2'), 1, 'EUR', build)
            release.set()
            self.assertEqual(first.result().resolve(), second.result().resolve())
        self.assertEqual(build.call_count, 1)

    def test_random_and_file_gene_selection_rerun(self):
        for selection in (GeneSelectionRequest('random', 1, 42, None),
                          GeneSelectionRequest('file', 0, 42, self.root / 'panel1.tsv')):
            request = replace(self.request, gene_selection=selection)
            first = prepare_chromosomes(request)
            with patch('rewrite.preprocessing.prepare.pgen_extract', side_effect=AssertionError('cache miss')):
                self.assertEqual(first, prepare_chromosomes(request))
            self.assertTrue((request.run_dir / 'selected_genes.tsv').is_file())


if __name__ == '__main__':
    unittest.main()
