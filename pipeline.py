"""
================================================================================
  CFD Slicer Pipeline
  Author  : Julian G-A
  Project : UTFR
--------------------------------------------------------------------------------
  Slice, interpolate and compare STAR-CCM+ EnSight Gold exports.
  Do not distribute without permission.
================================================================================
"""

import os
import sys
import time
import copy
import gc
import multiprocessing
import tempfile
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import numpy as np
import pyvista as pv

import config as cfg
import io_utils
import geometry as geom_mod
import interpolation as interp_mod
import plotting
import recruits


# ── timer / logger ───────────────────────────────────────────────
_T0 = time.perf_counter()

def log(msg, t_ref=None):
    elapsed = time.perf_counter() - _T0
    suffix  = f"  (+{(time.perf_counter()-t_ref)*1000:.0f}ms)" if t_ref else ""
    print(f"  [{elapsed:7.2f}s] {msg}{suffix}")
    sys.stdout.flush()


# ── helpers ───────────────────────────────────────────────────────
def _extract_plane_jobs(cfd_mesh, scalars_to_run, axes_to_run, staging_dir, config=cfg):
    """
    Return list of plane job dicts from the CFD mesh.
    Jobs carry paths to per-plane arrays instead of the entire CFD dataset.
    With staging_dir=None, inspect sizes without writing temporary arrays.
    """
    cfg = config
    jobs = []
    axis_counts = {0: 0, 1: 0, 2: 0}
    for block_idx in range(cfd_mesh.n_blocks):
        b = cfd_mesh[block_idx]
        if b is None or b.n_cells == 0:
            continue

        cc   = b.cell_centers()
        pts  = cc.points

        # ── Auto-detect the constant (slicing) axis ───────────────────────────
        # The slicing axis is the one with fewest unique coordinate values.
        # We also check the block name as a hint, but the data takes precedence.
        n_unique = [len(np.unique(np.round(pts[:, a], 5))) for a in range(3)]
        detected_ax = int(np.argmin(n_unique))

        # Block name hint — read the axis letter from the leading token so a
        # name like 'XY_Planes' or 'Z_Planes' resolves unambiguously instead of
        # matching whichever axis letter appears first anywhere in the string.
        block_name = cfd_mesh.get_block_name(block_idx).strip().upper()
        lead_token = block_name.split('_')[0]
        hint_ax = None
        for a, label in [(0, 'X'), (1, 'Y'), (2, 'Z')]:
            if lead_token == label or lead_token.startswith(label + 'PLANE'):
                hint_ax = a
                break
        if hint_ax is None:
            # Fall back to a substring check for less structured names.
            for a, label in [(0, 'X'), (1, 'Y'), (2, 'Z')]:
                if label in block_name:
                    hint_ax = a
                    break

        # Trust the block-name hint when present; the constant-axis heuristic
        # is fooled by slab aspect ratio. A stack of Y-slices is thin in Y but
        # spans the full car in X and Z, and if the car is shorter vertically
        # than it is long, Z can have fewer unique values than Y — so argmin
        # picks the wrong (spanning) axis and np.unique then returns one
        # "plane" per cell-center coordinate (hundreds of thousands of them).
        # The block name (e.g. 'Y_Planes') is authoritative; fall back to the
        # data heuristic only when the name carries no axis letter.
        if hint_ax is not None:
            ax = hint_ax
            if hint_ax != detected_ax:
                log(f"  [INFO] Block {block_idx} ('{cfd_mesh.get_block_name(block_idx).strip()}'): "
                    f"using name-hint axis {cfg.AXIS_LABEL[hint_ax]} "
                    f"(data heuristic guessed {cfg.AXIS_LABEL[detected_ax]}; "
                    f"unique counts X={n_unique[0]:,} Y={n_unique[1]:,} Z={n_unique[2]:,})")
        else:
            ax = detected_ax
            log(f"  [INFO] Block {block_idx} ('{cfd_mesh.get_block_name(block_idx).strip()}'): "
                f"no axis hint in name, detected axis {cfg.AXIS_LABEL[detected_ax]} "
                f"(unique counts X={n_unique[0]:,} Y={n_unique[1]:,} Z={n_unique[2]:,})")

        ha, va = cfg.PLOT_AXES[ax]

        if cfg.AXIS_LABEL[ax] not in axes_to_run:
            continue

        # Collect arrays for all requested scalars.
        # First try exact match, then fall back to partial/substring match
        # so the pipeline works regardless of how the .case file was exported.
        scalar_vals = {}
        for skey in scalars_to_run:
            arr_name = cfg.SCALARS[skey]['array']
            if arr_name in cc.array_names:
                # Exact match
                scalar_vals[skey] = cc[arr_name]
            else:
                # Partial match — find any array whose name contains arr_name
                # or whose name is contained within arr_name
                matches = [a for a in cc.array_names
                           if arr_name.lower() in a.lower()
                           or a.lower() in arr_name.lower()]
                # Deduplicate (pyvista sometimes returns arrays twice)
                matches = list(dict.fromkeys(matches))
                if matches:
                    scalar_vals[skey] = cc[matches[0]]
                    if matches[0] != arr_name:
                        log(f"  [{skey}] array '{arr_name}' not found — "
                            f"using '{matches[0]}'")
                else:
                    log(f"  [{skey}] WARNING: no array matching '{arr_name}' "
                        f"in block {block_idx} — skipping")
        if not scalar_vals:
            continue

        unique = np.unique(np.round(pts[:, ax], 5))

        # Filter planes to car region using CAR_BOUNDS_2D
        # For each slice axis, the slicing coordinate range is
        # the extent of the OTHER two axes' bounds
        # e.g. X-slice: X range derived from Y-slice h_min/h_max isn't right —
        # use a dedicated per-axis range based on geometry bounds
        car_extent = {
            0: (cfg.CAR_BOUNDS_2D[1][0], cfg.CAR_BOUNDS_2D[1][1]),  # X: from Y-slice h range
            1: (cfg.CAR_BOUNDS_2D[0][0], cfg.CAR_BOUNDS_2D[0][1]),  # Y: from X-slice h range
            2: (cfg.CAR_BOUNDS_2D[0][2], cfg.CAR_BOUNDS_2D[0][3]),  # Z: from X-slice v range
        }
        lo, hi = car_extent[ax]
        unique_filtered = unique[(unique >= lo) & (unique <= hi)]

        if len(unique_filtered) == 0:
            log(f"  [WARN] Block {block_idx} axis {cfg.AXIS_LABEL[ax]}: "
                f"all {len(unique)} planes outside car bounds [{lo:.3f}, {hi:.3f}]m — "
                f"data range [{unique.min():.3f}, {unique.max():.3f}]m. "
                f"Check CAR_BOUNDS_2D in Settings.")
            continue
        unique = unique_filtered

        log(f"  Block {block_idx} (axis {cfg.AXIS_LABEL[ax]}): "
            f"{len(unique)} planes  scalars={list(scalar_vals.keys())}")

        # Sort once instead of scanning the complete block for every plane.
        order = np.argsort(pts[:, ax], kind='stable')
        coordinates = pts[order, ax]
        h_min, h_max, v_min, v_max = cfg.CAR_BOUNDS_2D[ax]
        for plane_val in sorted(unique):
            start = np.searchsorted(coordinates, plane_val-1e-4, side='left')
            stop = np.searchsorted(coordinates, plane_val+1e-4, side='right')
            indices = order[start:stop]
            selected = pts[indices]
            keep = ((np.abs(selected[:, ax]-plane_val) < 1e-4) &
                    (selected[:, ha] >= h_min) & (selected[:, ha] <= h_max) &
                    (selected[:, va] >= v_min) & (selected[:, va] <= v_max))
            indices = np.sort(indices[keep])
            axis_counts[ax] += 1
            plane_idx = axis_counts[ax]
            path = os.path.join(staging_dir, f'plane_{len(jobs):06d}.npz') if staging_dir else None
            if staging_dir is not None:
                np.savez(path, points=pts[indices][:, [ha, va]],
                         values=np.column_stack([v[indices] for v in scalar_vals.values()]))
            jobs.append(dict(
                axis=ax, ha=ha, va=va,
                plane_idx=plane_idx,
                plane_val=float(plane_val),
                path=path, scalar_keys=list(scalar_vals), point_count=len(indices),
            ))
    return jobs


def inspect_recruit_memory(geom_case, data_case, scalars, axes):
    """Read actual input sizes for the same estimate used during execution."""
    geometry = io_utils.load_geometry(geom_case)
    geometry_bytes = sum(b.actual_memory_size * 1024 for b in geometry)
    mesh = io_utils.load_cfd(data_case)
    jobs = _extract_plane_jobs(mesh, scalars, axes, None)
    if not jobs:
        raise ValueError('No matching planes found for the selected scalars and axes.')
    return dict(point_count=max(j['point_count'] for j in jobs), geometry_bytes=geometry_bytes)


_RECRUIT_GEOMETRY = None


def _init_recruit(settings, geometry_path):
    global _RECRUIT_GEOMETRY
    for name, value in settings.items():
        setattr(cfg, name, value)
    _RECRUIT_GEOMETRY = list(pv.read(geometry_path))


def _process_plane(job, out_dir, fig_height, dpi):
    """Process one complete plane and return only small progress metadata."""
    ax, plane_idx, plane_val = job['axis'], job['plane_idx'], job['plane_val']
    ax_label = cfg.AXIS_LABEL[ax]
    val_mm = plane_val * 1000.0
    description = f'{ax_label}#{plane_idx:03d} {val_mm:+8.1f}mm'
    timings = dict(mask=0.0, interpolate=0.0, save=0.0, plot=0.0)
    if job['point_count'] < 4:
        return dict(description=description, skipped=True, timings=timings)
    with np.load(job['path'], allow_pickle=False) as data:
        points, values = data['points'], data['values']
    h_min, h_max, v_min, v_max = cfg.CAR_BOUNDS_2D[ax]
    px = 1000.0 / cfg.RESOLUTION_MM
    stem = f'slice_{ax_label}_{plane_idx:03d}_{val_mm:+.1f}mm'
    t = time.perf_counter()
    mask = geom_mod.build_mask(_RECRUIT_GEOMETRY, ax, plane_val,
                              h_min, h_max, v_min, v_max, px, plane_pts_2d=points)
    timings['mask'] = time.perf_counter()-t
    t = time.perf_counter()
    grids = interp_mod.interpolate_per_region(points, values, mask, h_min, v_min, px)
    scalar_grids = {key: grids[:, :, i] for i, key in enumerate(job['scalar_keys'])}
    timings['interpolate'] = time.perf_counter()-t
    t = time.perf_counter()
    extent = np.array([h_min, h_max, v_min, v_max], dtype=np.float32)
    npz_path = io_utils.save_plane(os.path.join(out_dir, 'NPZ', ax_label), stem,
                                  scalar_grids, mask, extent, ax, plane_val, plane_idx)
    timings['save'] = time.perf_counter()-t
    t = time.perf_counter()
    import matplotlib.pyplot as plt
    h_label, v_label = cfg.PLOT_LABELS[ax]
    for key, grid in scalar_grids.items():
        fig = plotting.make_slice_plot(grid, mask, ax, plane_idx, val_mm,
                                      h_label, v_label, extent, key,
                                      fig_height=fig_height, dpi=dpi)
        try:
            png_dir = os.path.join(out_dir, 'PNG', key, ax_label)
            os.makedirs(png_dir, exist_ok=True)
            fig.savefig(os.path.join(png_dir, stem+'.png'),
                        bbox_inches='tight', facecolor=cfg.BG_COLOR, dpi=dpi)
        finally:
            plt.close(fig)
    timings['plot'] = time.perf_counter()-t
    return dict(description=description, skipped=False, timings=timings,
                npz_kb=os.path.getsize(npz_path)/1024)


def _run_jobs(jobs, count, settings, geometry_path, out_dir, fig_height, dpi):
    """Keep at most one pending plane per recruit."""
    with ProcessPoolExecutor(max_workers=count,
                             mp_context=multiprocessing.get_context('spawn'),
                             initializer=_init_recruit,
                             initargs=(settings, geometry_path)) as pool:
        remaining = iter(jobs)
        pending = {}
        def submit_one():
            job = next(remaining, None)
            if job is not None:
                pending[pool.submit(_process_plane, job, out_dir, fig_height, dpi)] = job
        for _ in range(count):
            submit_one()
        try:
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    job = pending.pop(future)
                    try:
                        result = future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            f"Recruit failed on {cfg.AXIS_LABEL[job['axis']]}"
                            f"#{job['plane_idx']:03d}: {exc}") from exc
                    os.remove(job['path'])
                    yield result
                    submit_one()
        finally:
            for future in pending:
                future.cancel()


def run_pipeline(scalars=None, axes=None, geom_case=None, data_case=None,
                 output_dir=None, resolution_mm=None, fig_height=7.0, dpi=None,
                 ram_budget_gb=None, max_recruits=recruits.MAX_RECRUITS,
                 progress_callback=None):
    """Process planes using RAM/CPU-limited spawned recruits.

    RAM is an estimated budget in GiB, not an enforced memory cap. Callbacks
    receive allocation, plane and done dictionaries on the calling thread.
    Use a main guard when calling from another script on Windows.
    """
    global _T0
    _T0 = time.perf_counter()
    settings = {name: copy.deepcopy(value) for name, value in vars(cfg).items() if name.isupper()}
    scalars_to_run = list(dict.fromkeys(scalars if scalars is not None else cfg.SCALARS))
    axes_to_run = list(dict.fromkeys(axes if axes is not None else ['X', 'Y', 'Z']))
    if not scalars_to_run or any(s not in settings['SCALARS'] for s in scalars_to_run):
        raise ValueError('Select at least one configured scalar.')
    if not axes_to_run or any(a not in ('X', 'Y', 'Z') for a in axes_to_run):
        raise ValueError('Select at least one of X, Y, Z.')
    res_mm = cfg.RESOLUTION_MM if resolution_mm is None else resolution_mm
    dpi = cfg.DPI if dpi is None else dpi
    settings['RESOLUTION_MM'], settings['DPI'] = res_mm, dpi
    out_dir = output_dir if output_dir is not None else cfg.OUTPUT_DIR
    if not out_dir:
        raise ValueError('Choose an output directory.')
    out_dir = os.path.abspath(out_dir)
    bounds = [settings['CAR_BOUNDS_2D'][a] for a in range(3)
              if settings['AXIS_LABEL'][a] in axes_to_run]
    estimate = recruits.estimate_recruit_bytes(bounds, res_mm, len(scalars_to_run),
                                              fig_height=fig_height, dpi=dpi)
    budget = recruits.default_budget_gib() if ram_budget_gb is None else ram_budget_gb
    recruits.allocate_recruits(budget, estimate, max_recruits=max_recruits)
    os.makedirs(out_dir, exist_ok=True)
    def emit(event):
        if progress_callback is not None:
            progress_callback(event)
    log('='*70)
    log(f'CFD Slicer | Scalars: {scalars_to_run} | Axes: {axes_to_run}')
    log(f'Resolution: {res_mm} mm/px | Output: {out_dir}')
    log('STEP 1/4 — Loading geometry and CFD data ...')
    geom_blocks = io_utils.load_geometry(geom_case)
    geometry_bytes = sum(b.actual_memory_size*1024 for b in geom_blocks)
    cfd_mesh = io_utils.load_cfd(data_case)
    totals = dict(mask=0.0, interpolate=0.0, save=0.0, plot=0.0)
    completed = skipped = 0
    with tempfile.TemporaryDirectory(prefix='cfd-recruits-') as staging:
        geometry_path = os.path.join(staging, 'geometry.vtm')
        pv.MultiBlock(geom_blocks).save(geometry_path)
        log('STEP 2/4 — Staging plane data ...')
        from types import SimpleNamespace
        jobs = _extract_plane_jobs(cfd_mesh, scalars_to_run, axes_to_run, staging,
                                   config=SimpleNamespace(**settings))
        del cfd_mesh, geom_blocks
        gc.collect()
        total = len(jobs)
        if not total:
            raise ValueError('No matching planes found. Check axes, scalar arrays and car bounds.')
        point_count = max(j['point_count'] for j in jobs)
        estimate = recruits.estimate_recruit_bytes(
            bounds, res_mm, len(scalars_to_run), point_count=point_count,
            geometry_bytes=geometry_bytes, fig_height=fig_height, dpi=dpi)
        allocation = recruits.allocate_recruits(budget, estimate,
                                               max_recruits=max_recruits, plane_count=total)
        log(f'Recruits: {allocation.recruits} | RAM budget: {budget:.2f} GiB | '
            f'usable now: {allocation.usable_gib:.2f} GiB | '
            f'estimate/recruit: {allocation.per_recruit_gib:.2f} GiB | '
            f'system reserve: {allocation.reserve_gib:.2f} GiB')
        emit(dict(event='allocation', recruits=allocation.recruits, total=total,
                  budget_gib=budget, per_recruit_gib=allocation.per_recruit_gib,
                  point_count=point_count, geometry_bytes=geometry_bytes))
        if allocation.recruits == 0:
            raise MemoryError(
                f'Not enough RAM budget or available memory for one recruit '
                f'(estimated {allocation.per_recruit_gib:.2f} GiB). '
                'Increase the budget, free memory, or use a coarser resolution.')
        log(f'STEP 3+4/4 — Processing {total} planes with {allocation.recruits} recruits ...')
        for processed, result in enumerate(_run_jobs(
                jobs, allocation.recruits, settings, geometry_path, out_dir, fig_height, dpi), start=1):
            skipped += int(result['skipped'])
            completed += int(not result['skipped'])
            for key in totals:
                totals[key] += result['timings'][key]
            status = 'SKIP (<4 points)' if result['skipped'] else f"{result['npz_kb']:.0f} KB"
            log(f"[{processed}/{total}] {result['description']}  {status}")
            emit(dict(event='plane', processed=processed, total=total, completed=completed,
                      skipped=skipped, recruits=allocation.recruits))
    wall_seconds = time.perf_counter()-_T0
    summary = dict(event='done', total=total, completed=completed, skipped=skipped,
                   recruits=allocation.recruits, wall_seconds=wall_seconds)
    log(f'DONE — {completed} completed, {skipped} skipped | {out_dir}')
    log('Cumulative recruit time (overlaps when recruits run in parallel):')
    for key, seconds in totals.items():
        log(f'  {key:12s}: {seconds:.1f}s')
    log(f'  Wall clock  : {wall_seconds:.1f}s')
    emit(summary)
    return summary


