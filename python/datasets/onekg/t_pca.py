import importlib.util
import struct
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('pca_utils', Path(__file__).with_name('utils.py'))
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)


class PCAMergeTest(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.sources = [self.root / 'chr1', self.root / 'chr2']
        self.work = self.root / 'pca'
        self.work.mkdir()

    def merged(self, prefix, variants):
        Path(f'{prefix}.pgen').write_bytes(struct.pack('<3sIIB', b'\x6c\x1b\x10', variants, 50 if variants else 0, 0))
        Path(f'{prefix}.pvar').write_text('variants')
        Path(f'{prefix}.psam').write_text('samples')

    def finish(self, command, check):
        self.assertTrue(check)
        prefix = Path(command[command.index('--out') + 1])
        if '--pmerge-list' in command:
            self.assertEqual(prefix.name, 'merged.partial')
            self.assertFalse((self.work / 'merged.pgen').exists())
            self.assertIn('--merge-info-mode', command)
            self.merged(prefix, 100)
            Path(f'{prefix}.log').write_text('merge completed')
        elif '--indep-pairwise' in command:
            Path(f'{prefix}.prune.in').write_text('selected')
        else:
            Path(f'{prefix}.eigenvec').write_text('new PCs')
            Path(f'{prefix}.eigenval').write_text('new eigenvalues')

    def test_incomplete_merge_rebuilds_and_invalidates_downstream_outputs(self):
        self.merged(self.work / 'merged', 0)
        for name in ('pruned.prune.in', 'pruned.prune.out', 'pca.eigenvec', 'pca.eigenval'):
            (self.work / name).write_text('stale')
        with patch.object(utils.subprocess, 'run', side_effect=self.finish) as run:
            output = utils.create_pca(self.sources, self.work)
        self.assertEqual(run.call_count, 3)
        self.assertEqual(output.read_text(), 'new PCs')
        self.assertEqual((self.work / 'merged.log').read_text(), 'merge completed')
        with patch.object(utils.subprocess, 'run', side_effect=AssertionError('unexpected regeneration')):
            self.assertEqual(utils.create_pca(self.sources, self.work), output)

    def test_failed_merge_is_not_reused_on_retry(self):
        self.merged(self.work / 'merged', 0)
        def fail(command, check):
            self.merged(Path(command[-1]), 0)
            raise subprocess.CalledProcessError(5, command)
        with patch.object(utils.subprocess, 'run', side_effect=fail) as run:
            with self.assertRaises(subprocess.CalledProcessError):
                utils.create_pca(self.sources, self.work)
        self.assertEqual(run.call_count, 1)
        self.assertFalse((self.work / 'merged.pgen').exists())
        with patch.object(utils.subprocess, 'run', side_effect=self.finish) as run:
            utils.create_pca(self.sources, self.work)
        self.assertEqual(run.call_count, 3)
        self.assertFalse((self.work / 'merged.partial.pgen').exists())

    def test_missing_member_rebuilds_and_single_chromosome_does_not_merge(self):
        self.merged(self.work / 'merged', 100)
        (self.work / 'merged.psam').unlink()
        with patch.object(utils.subprocess, 'run', side_effect=self.finish) as run:
            utils.create_pca(self.sources, self.work)
        self.assertEqual(run.call_count, 3)
        single = self.root / 'single'
        with patch.object(utils.subprocess, 'run', side_effect=self.finish) as run:
            utils.create_pca(self.sources[:1], single)
        self.assertEqual(run.call_count, 2)
        self.assertTrue(all('--pmerge-list' not in call.args[0] for call in run.call_args_list))

    def test_real_plink_merge_preserves_pca_without_info_column(self):
        import numpy as np
        import pgenlib

        rng = np.random.default_rng(42)
        for chromosome, prefix in enumerate(self.sources, 1):
            with pgenlib.PgenWriter(bytes(prefix) + b'.pgen', 80, variant_ct=40, nonref_flags=False) as writer:
                writer.append_biallelic_batch(rng.binomial(2, 0.25, size=(40, 80)).astype(np.int8))
            Path(f'{prefix}.psam').write_text('#IID\n' + ''.join(f's{i:03}\n' for i in range(80)))
            Path(f'{prefix}.pvar').write_text(
                '##INFO=<ID=NOTE,Number=1,Type=String,Description="Unused annotation">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tINFO\n' + ''.join(
                    f'{chromosome}\t{i + 1}\t{chromosome}:{i + 1}:A:G\tA\tG\tNOTE=unused\n' for i in range(40)))
        original = {path: path.read_bytes() for prefix in self.sources
                    for path in (Path(f'{prefix}.pgen'), Path(f'{prefix}.pvar'), Path(f'{prefix}.psam'))}
        self.merged(self.work / 'merged', 0)
        repaired = utils.create_pca(self.sources, self.work, num_components=2)
        self.assertNotIn('\tINFO', (self.work / 'merged.pvar').read_text())
        baseline = self.root / 'baseline'
        baseline.mkdir()
        subprocess.run(['plink2', '--pmerge-list', str(self.work / 'pmerge_list.txt'), '--out', str(baseline / 'merged')],
                       check=True, stdout=subprocess.DEVNULL)
        expected = utils.create_pca(self.sources, baseline, num_components=2)
        np.testing.assert_allclose(np.loadtxt(repaired.with_suffix('.eigenval')), np.loadtxt(expected.with_suffix('.eigenval')))
        np.testing.assert_allclose(abs(np.loadtxt(repaired, skiprows=1, usecols=(1, 2))),
                                   abs(np.loadtxt(expected, skiprows=1, usecols=(1, 2))), atol=1e-8)
        with patch.object(utils.subprocess, 'run', side_effect=AssertionError('unexpected regeneration')):
            utils.create_pca(self.sources, self.work, num_components=2)
        self.assertEqual(original, {path: path.read_bytes() for path in original})


if __name__ == '__main__':
    unittest.main()
