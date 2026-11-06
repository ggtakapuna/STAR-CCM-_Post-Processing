"""
================================================================================
  CFD Slicer Pipeline
  Author  : Julian G-A
  Contributors: John Baton
  Project : UTFR
--------------------------------------------------------------------------------
  Checks Y reflection of displayed grids, geometry overlays and surface meshes.
================================================================================
"""
import unittest
from unittest.mock import patch

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt

import config as cfg
import plotting
import surface_render


class DisplayOrientationTests(unittest.TestCase):
    def test_reflects_each_slice_axis_and_coordinates(self):
        grid = np.arange(6).reshape(2, 3)
        extent = np.array([1., 4., 2., 5.])
        with patch.object(cfg, 'REFLECT_Y_IMAGES', True):
            for axis, expected, bounds, value in (
                    (0, grid[:, ::-1], [-4, -1, 2, 5], 30),
                    (1, grid, [1, 4, 2, 5], -30),
                    (2, grid[::-1], [1, 4, -5, -2], 30)):
                actual_extent, actual_value, (actual,) = plotting._orient_slice_for_display(
                    axis, extent, 30, grid)
                np.testing.assert_array_equal(actual, expected)
                np.testing.assert_array_equal(actual_extent, bounds)
                self.assertEqual(actual_value, value)
        np.testing.assert_array_equal(grid, np.arange(6).reshape(2, 3))
        np.testing.assert_array_equal(extent, [1, 4, 2, 5])

    def test_slice_and_delta_geometry_stay_aligned(self):
        grid = np.arange(6, dtype=float).reshape(2, 3)
        mask = np.zeros_like(grid, dtype=bool)
        mask[0, 0] = True
        other_mask = np.zeros_like(mask)
        other_mask[1, 2] = True
        with patch.object(cfg, 'REFLECT_Y_IMAGES', True):
            for axis in (0, 1, 2):
                args = (axis, 1, 30, *cfg.PLOT_LABELS[axis], [1, 4, 2, 5], 'Cp')
                single = plotting.make_slice_plot(grid, mask, *args)
                delta = plotting.make_diff_plot(grid, grid+1, mask, other_mask, *args)
                try:
                    _, _, (g, m, mb) = plotting._orient_slice_for_display(
                        axis, [1, 4, 2, 5], 30, grid, mask, other_mask)
                    np.testing.assert_array_equal(single.axes[0].images[0].get_array(), g)
                    np.testing.assert_array_equal(single.axes[0].images[1].get_array(),
                                                  plotting._solid_overlay(m))
                    np.testing.assert_array_equal(delta.axes[0].images[1].get_array(),
                                                  plotting._diff_geometry_overlay(m, mb))
                    if axis == 1:
                        self.assertIn('-30.0 mm', single.axes[0].get_title())
                finally:
                    plt.close(single)
                    plt.close(delta)

    def test_surface_reflects_y_and_keeps_scalar_values_and_original_mesh(self):
        surf = pv.Plane(center=(1, 2, 3), i_resolution=2, j_resolution=2)
        surf.cell_data['Cp'] = np.arange(surf.n_cells, dtype=float)
        edges = surf.extract_all_edges()
        parts = [dict(name='asymmetric', surf=surf, edges=edges)]
        original_points = surf.points.copy()
        with patch.object(cfg, 'REFLECT_Y_IMAGES', True):
            result = surface_render._parts_for_display(parts)
        for key in ('surf', 'edges'):
            expected = parts[0][key].points.copy()
            expected[:, 1] *= -1
            np.testing.assert_allclose(result[0][key].points, expected)
        np.testing.assert_array_equal(result[0]['surf']['Cp'], surf['Cp'])
        np.testing.assert_array_equal(surf.points, original_points)

    def test_reflection_can_be_disabled(self):
        grid = np.arange(6).reshape(2, 3)
        with patch.object(cfg, 'REFLECT_Y_IMAGES', False):
            extent, val, (actual,) = plotting._orient_slice_for_display(0, [1, 4, 2, 5], 30, grid)
            self.assertIs(actual, grid)
            self.assertEqual(val, 30)
            self.assertEqual(extent, [1, 4, 2, 5])
            parts = [dict(surf=None, edges=None)]
            self.assertIs(surface_render._parts_for_display(parts), parts)


if __name__ == '__main__':
    unittest.main()
