"""
================================================================================
  CFD Slicer Pipeline
  Author  : Julian G-A
  Contributors: John Baton
  Project : UTFR
--------------------------------------------------------------------------------
  test_compare_recruits.py
  Regression checks for comparison recruits.
  Checks RAM estimates, serial/parallel output consistency and aero-only processing.
  Do not distribute without permission.
================================================================================
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import compare
import io_utils
import recruits


class ComparisonRecruitsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.a, self.b = self.root/'a', self.root/'b'
        grid = np.arange(441, dtype=np.float32).reshape(21, 21)/441
        mask = np.zeros_like(grid, dtype=bool)
        mask[10, 10] = True
        for idx in range(1, 4):
            for directory, offset in ((self.a, 0), (self.b, 0.1)):
                io_utils.save_plane(str(directory/'NPZ'/'X'),
                    f'slice_X_{idx:03d}_{idx*10:+.1f}mm',
                    {'Cp': grid+offset, 'CpT': 2*grid+offset},
                    mask, [-0.1, 0.1, -0.1, 0.1], 0, idx/100, idx)
        resources = patch('recruits.system_resources', return_value=
                          recruits.Resources(32*recruits.GIB, 20*recruits.GIB, 8))
        resources.start()
        self.addCleanup(resources.stop)

    def run_compare(self, name, count, **kwargs):
        return compare.compare_runs(str(self.a), str(self.b), scalars=['Cp', 'CpT'],
            axes=['X'], output_dir=str(self.root/name), fig_height=2, dpi=40,
            max_recruits=count, **kwargs)

    def test_parallel_matches_serial_and_analysis_estimate(self):
        _, _, preview = compare.prepare_comparison(str(self.a), str(self.b), ['Cp', 'CpT'], ['X'], 2, 40)
        serial = self.run_compare('serial', 1, ram_budget_gb=4)
        events = []
        with ThreadPoolExecutor(max_workers=1) as background:
            parallel = background.submit(self.run_compare, 'parallel', 2,
                ram_budget_gb=4, progress_callback=events.append).result()
        self.assertEqual(serial, parallel)
        self.assertEqual(len(parallel), 6)
        self.assertEqual(events[0]['recruits'], 2)
        self.assertEqual(events[0]['per_recruit_gib'], preview/recruits.GIB)
        self.assertEqual(events[-1]['event'], 'done')
        self.assertEqual(events[-1]['completed'], 3)
        for row in parallel:
            self.assertAlmostEqual(row['mean'], 0.1, places=6)
        self.assertEqual((self.root/'serial'/'comparison_summary.csv').read_bytes(),
                         (self.root/'parallel'/'comparison_summary.csv').read_bytes())
        pngs = list((self.root/'serial'/'PNG').rglob('*.png'))
        self.assertEqual(len(pngs), 6)
        for path in pngs:
            other = self.root/'parallel'/path.relative_to(self.root/'serial')
            self.assertEqual(path.read_bytes(), other.read_bytes())

    def test_insufficient_ram_and_aero_only_bypass(self):
        with self.assertRaisesRegex(MemoryError, 'comparison recruit'):
            self.run_compare('low', 2, ram_budget_gb=0.01)
        with patch('recruits.run_tasks', side_effect=AssertionError('Must not start pool')):
            self.assertEqual(self.run_compare('aero-only', 2, images=False), [])


if __name__ == '__main__':
    unittest.main()
