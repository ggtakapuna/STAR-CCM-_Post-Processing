"""
================================================================================
  CFD Slicer Pipeline
  Author  : Julian G-A
  Contributors: John Baton
  Project : UTFR
--------------------------------------------------------------------------------
  test_recruits.py
  Regression checks for recruit allocation and the slice pipeline.
  Checks memory limits, shared-field interpolation, spawned processing and cleanup.
  Do not distribute without permission.
================================================================================
"""

import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pyvista as pv
from scipy.interpolate import griddata

import config as cfg
import interpolation
import pipeline
import recruits


class AllocationTests(unittest.TestCase):
    def test_budget_cpu_plane_limits_and_system_reserve(self):
        resources = recruits.Resources(32*recruits.GIB, 20*recruits.GIB, 4)
        allocate = lambda budget, **kw: recruits.allocate_recruits(
            budget, recruits.GIB, resources=resources, **kw)
        self.assertEqual(allocate(2).recruits, 2)
        self.assertEqual(allocate(16).recruits, 4)
        self.assertEqual(allocate(16, plane_count=1).recruits, 1)
        self.assertEqual(allocate(16, max_recruits=1).recruits, 1)
        self.assertEqual(allocate(0.25).recruits, 0)
        self.assertEqual(allocate(16, plane_count=0).recruits, 0)
        for total, available, expected in ((32, 5, 0), (8, 3, 1)):
            self.assertEqual(recruits.allocate_recruits(16, recruits.GIB,
                resources=recruits.Resources(total*recruits.GIB, available*recruits.GIB, 16)).recruits,
                expected)
        self.assertEqual(recruits.allocate_recruits(32, recruits.GIB,
            resources=recruits.Resources(64*recruits.GIB, 40*recruits.GIB, 64)).recruits, 8)

    def test_invalid_budgets_and_limits(self):
        for value in (0, -1, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                recruits.allocate_recruits(value, recruits.GIB)
        for count in (0, 9, 1.5):
            with self.assertRaises(ValueError):
                recruits.allocate_recruits(4, recruits.GIB, max_recruits=count)

    def test_estimate_scales_with_workload(self):
        bounds = [(0, 1, 0, 1)]
        baseline = recruits.estimate_recruit_bytes(bounds, 5, 1)
        self.assertGreater(recruits.estimate_recruit_bytes(bounds, 1, 1), baseline)
        self.assertGreater(recruits.estimate_recruit_bytes(bounds, 5, 6), baseline)
        self.assertGreater(recruits.estimate_recruit_bytes(bounds, 5, 1, dpi=300), baseline)
        self.assertGreater(recruits.estimate_recruit_bytes(
            bounds, 5, 1, geometry_bytes=recruits.GIB), baseline)


class InterpolationTests(unittest.TestCase):
    def test_combined_fields_match_original_scalar_interpolation(self):
        rng = np.random.default_rng(26)
        points = rng.uniform(0, 1, (400, 2))
        values = np.column_stack((points[:, 0]+points[:, 1],
                                 np.sin(points[:, 0]*4)-points[:, 1]))
        values[8, 1] = np.nan
        mask = np.zeros((31, 31), dtype=bool)
        mask[:, 15] = True
        result = interpolation.interpolate_per_region(points, values, mask, 0, 0, 30)
        for field in range(2):
            expected = interpolation.interpolate_per_region(points, values[:, field], mask, 0, 0, 30)
            np.testing.assert_allclose(result[:, :, field], expected, equal_nan=True)
        self.assertTrue(np.isnan(result[mask]).all())
        mask[:] = False
        rows, cols = np.indices(mask.shape)
        expected = griddata(points, values, np.column_stack((cols.ravel()/30, rows.ravel()/30)))
        result = interpolation.interpolate_per_region(points, values, mask, 0, 0, 30)
        np.testing.assert_allclose(result, expected.reshape(result.shape),
                                   rtol=1e-6, atol=1e-7, equal_nan=True)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bounds = {a: (-0.2, 0.2, -0.2, 0.2) for a in range(3)}
        scalars = copy.deepcopy({k: cfg.SCALARS[k] for k in ('Cp', 'CpT')})
        scalars['Cp']['array'] = 'pressure_test'
        scalars['CpT']['array'] = 'total_test'
        config_patch = patch.multiple(cfg, CAR_BOUNDS_2D=self.bounds, SCALARS=scalars, DILATION_PX=2)
        config_patch.start()
        self.addCleanup(config_patch.stop)
        pv.MultiBlock([pv.Sphere(radius=0.07)]).save(self.root/'geometry.vtm')
        blocks = pv.MultiBlock()
        for i, x in enumerate((-0.05, 0.05, 0.15)):
            plane = pv.Plane(center=(x, 0, 0), direction=(1, 0, 0),
                             i_size=0.38, j_size=0.38, i_resolution=8, j_resolution=8)
            points = plane.cell_centers().points
            plane.cell_data['pressure_test'] = points[:, 1]+points[:, 2]
            plane.cell_data['total_test'] = 2*points[:, 1]-points[:, 2]
            blocks[f'X_Planes_{i}'] = plane
        blocks.save(self.root/'data.vtm')
        resources = patch('recruits.system_resources', return_value=
                          recruits.Resources(32*recruits.GIB, 20*recruits.GIB, 4))
        resources.start()
        self.addCleanup(resources.stop)

    def run_pipeline(self, output, count=2, **kwargs):
        return pipeline.run_pipeline(scalars=['Cp', 'CpT'], axes=['X'],
            geom_case=str(self.root/'geometry.vtm'), data_case=str(self.root/'data.vtm'),
            output_dir=str(output), resolution_mm=20, fig_height=2, dpi=40,
            max_recruits=count, **kwargs)

    def test_windows_spawn_outputs_match_serial_and_settings_propagate(self):
        profile = pipeline.inspect_recruit_memory(str(self.root/'geometry.vtm'),
            str(self.root/'data.vtm'), ['Cp', 'CpT'], ['X'])
        preview = recruits.estimate_recruit_bytes([self.bounds[0]], 20, 2,
                                                 fig_height=2, dpi=40, **profile)
        original_resolution = cfg.RESOLUTION_MM
        events = []
        serial = self.run_pipeline(self.root/'serial', count=1, ram_budget_gb=4)
        with ThreadPoolExecutor(max_workers=1) as background:
            parallel = background.submit(self.run_pipeline, self.root/'parallel',
                ram_budget_gb=4, progress_callback=events.append).result()
        self.assertEqual(serial['recruits'], 1)
        self.assertEqual(parallel['recruits'], 2)
        self.assertEqual(parallel['completed'], 3)
        self.assertEqual(events[0]['per_recruit_gib'], preview/recruits.GIB)
        self.assertEqual(cfg.RESOLUTION_MM, original_resolution)
        self.assertEqual([e['event'] for e in events], ['allocation', 'plane', 'plane', 'plane', 'done'])
        files = list((self.root/'serial'/'NPZ'/'X').glob('*.npz'))
        self.assertEqual(len(files), 3)
        indices = []
        for path in files:
            other = self.root/'parallel'/path.relative_to(self.root/'serial')
            with np.load(path) as a, np.load(other) as b:
                self.assertEqual(set(a.files), set(b.files))
                for key in a.files:
                    np.testing.assert_array_equal(a[key], b[key])
                self.assertEqual(float(b['resolution_mm']), 20)
                self.assertEqual(b['mask'].shape, (21, 21))
                indices.append(int(b['plane_index']))
                self.assertTrue(np.isfinite(b['Cp']).any())
        self.assertEqual(sorted(indices), [1, 2, 3])
        pngs = list((self.root/'parallel'/'PNG').rglob('*.png'))
        self.assertEqual(len(pngs), 6)
        for path in pngs:
            serial_path = self.root/'serial'/path.relative_to(self.root/'parallel')
            self.assertEqual(path.read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
            self.assertEqual(path.read_bytes(), serial_path.read_bytes())

    def test_insufficient_budget_and_recruit_failure_clean_up_staging(self):
        original_temp = tempfile.TemporaryDirectory
        staged = []
        def track_temp(*args, **kwargs):
            result = original_temp(*args, **kwargs)
            staged.append(Path(result.name))
            return result
        with patch('pipeline.tempfile.TemporaryDirectory', side_effect=track_temp):
            with self.assertRaisesRegex(MemoryError, 'one recruit'):
                self.run_pipeline(self.root/'too-small', ram_budget_gb=0.01)
            with patch.dict(cfg.SCALARS['Cp'], cmap='invalid-colormap-for-test'):
                with self.assertRaisesRegex(RuntimeError, 'Recruit failed on X#'):
                    self.run_pipeline(self.root/'error', ram_budget_gb=4)
        self.assertTrue(staged)
        self.assertTrue(all(not path.exists() for path in staged))

    def test_staged_stack_matches_original_plane_selection(self):
        stack = pv.read(self.root/'data.vtm').combine()
        blocks = pv.MultiBlock({'X_Planes': stack})
        cc = stack.cell_centers()
        with tempfile.TemporaryDirectory() as staging:
            jobs = pipeline._extract_plane_jobs(blocks, ['Cp', 'CpT'], ['X'], staging)
            self.assertEqual(len(jobs), 3)
            for job in jobs:
                pts = cc.points
                selected = ((np.abs(pts[:, 0]-job['plane_val']) < 1e-4) &
                    (pts[:, 1] >= -0.2) & (pts[:, 1] <= 0.2) &
                    (pts[:, 2] >= -0.2) & (pts[:, 2] <= 0.2))
                self.assertNotIn('pts_all', job)
                with np.load(job['path']) as data:
                    np.testing.assert_array_equal(data['points'], pts[selected][:, [1, 2]])
                    np.testing.assert_array_equal(data['values'][:, 0], cc['pressure_test'][selected])


if __name__ == '__main__':
    unittest.main()
