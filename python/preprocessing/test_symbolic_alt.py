import tempfile
import unittest
from pathlib import Path

import numpy as np
import pgenlib

from .plink import pgen_allele_counts, pgen_extract


class SymbolicALTTest(unittest.TestCase):
    def test_dosages_counts_and_sample_order(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'chr21'
            alts = ('G', '<DEL>', '<INS:ME:ALU>', '<DEL:ME:LINE1>')
            keys = tuple(f'21:{i+1}:A:{alt}' for i, alt in enumerate(alts))
            dosages = np.array([[0, 1, 2, -9], [1, 0, -9, 2],
                                [0, 2, 1, -9], [2, -9, 0, 1]], dtype=np.int8)
            with pgenlib.PgenWriter(str(prefix).encode() + b'.pgen', 4, variant_ct=4) as writer:
                for row in dosages:
                    writer.append_biallelic(row)
            Path(f'{prefix}.psam').write_text('#IID\ns0\ns1\ns2\ns3\n')
            Path(f'{prefix}.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
                f'21\t{i+1}\t{key}\tA\t{alt}\n' for i, (key, alt) in enumerate(zip(keys, alts))))
            expected = dosages.copy()
            expected[expected == -9] = 0
            for indices in ((3, 1, 0, 2), (3, 1)):
                ids = tuple(f's{i}' for i in indices)
                matrix, emitted = pgen_extract(prefix, ids, keys)
                self.assertEqual(emitted, keys)
                np.testing.assert_array_equal(matrix, expected[:, indices].T)
                self.assertEqual(pgen_allele_counts(prefix, ids, keys),
                                 dict(zip(keys, expected[:, indices].sum(axis=1).tolist())))
            # Full ALT matching still rejects a malformed key.
            path = Path(f'{prefix}.pvar')
            path.write_text(path.read_text().replace(keys[2], '21:3:A:<INS:ME:LINE1>'))
            with self.assertRaisesRegex(ValueError, "counted allele is not the key's ALT"):
                pgen_extract(prefix, ('s0',), ('21:3:A:<INS:ME:LINE1>',))


if __name__ == '__main__':
    unittest.main()
