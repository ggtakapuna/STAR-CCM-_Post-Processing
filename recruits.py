"""RAM estimates and CPU limits for independent plane-processing recruits.

The budget is advisory: SciPy, VTK and plotting allocations vary with the mesh.
It excludes the parent process; available RAM is sampled again after staging.
"""

from dataclasses import dataclass
import math
import os
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED

import psutil

GIB = 1024 ** 3
MAX_RECRUITS = 8


def _init_settings(settings):
    import config
    for name, value in settings.items():
        setattr(config, name, value)


def run_tasks(function, jobs, count, settings):
    """Bound queued tasks; propagate failures before accepting further work."""
    with ProcessPoolExecutor(max_workers=count,
                             mp_context=multiprocessing.get_context('spawn'),
                             initializer=_init_settings, initargs=(settings,)) as pool:
        remaining = iter(enumerate(jobs))
        pending = {}
        def submit():
            item = next(remaining, None)
            if item is not None:
                index, job = item
                pending[pool.submit(function, job)] = index
        for _ in range(count):
            submit()
        try:
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    index = pending.pop(future)
                    yield index, future.result()
                    submit()
        finally:
            for future in pending:
                future.cancel()


@dataclass(frozen=True)
class Resources:
    total: int
    available: int
    cpus: int


@dataclass(frozen=True)
class Allocation:
    recruits: int
    budget_gib: float
    usable_gib: float
    per_recruit_gib: float
    reserve_gib: float


def system_resources():
    memory = psutil.virtual_memory()
    try:
        cpus = len(psutil.Process().cpu_affinity())
    except (AttributeError, NotImplementedError, psutil.Error):
        cpus = os.cpu_count() or 1
    return Resources(memory.total, memory.available, max(1, cpus))


def reserve_bytes(resources):
    return max(2 * GIB, int(resources.total * 0.15))


def default_budget_gib(resources=None):
    resources = resources or system_resources()
    return max(0.25, min(8.0, (resources.available - reserve_bytes(resources)) / GIB))


def estimate_recruit_bytes(bounds, resolution_mm, scalar_count, *,
                           point_count=100_000, geometry_bytes=0,
                           fig_height=7.0, dpi=150):
    """Conservative estimate for grids, triangulation, mesh and image buffers.

    Before loading, assume 100k source points and no known geometry size.
    At launch the largest staged plane and loaded geometry replace these values.
    """
    if not math.isfinite(resolution_mm) or resolution_mm <= 0:
        raise ValueError("Resolution must be a positive finite number.")
    if not math.isfinite(fig_height) or fig_height <= 0 or dpi <= 0:
        raise ValueError("Figure height and DPI must be positive finite numbers.")
    grid_peak = image_peak = 0
    for h_min, h_max, v_min, v_max in bounds:
        if not all(math.isfinite(v) for v in (h_min, h_max, v_min, v_max)):
            raise ValueError("Plot bounds must be finite.")
        if h_max <= h_min or v_max <= v_min:
            raise ValueError("Plot bounds must have positive width and height.")
        width = int((h_max - h_min) * 1000 / resolution_mm) + 1
        height = int((v_max - v_min) * 1000 / resolution_mm) + 1
        grid_peak = max(grid_peak, width * height * (128 + 24 * scalar_count))
        fig_width = min(24.0, max(5.0, fig_height * (h_max-h_min) / (v_max-v_min) + 1.8))
        image_peak = max(image_peak, fig_width * fig_height * dpi ** 2 * 24)
    source = point_count * (256 + 24 * scalar_count)
    return int(1.3 * (0.35 * GIB + grid_peak + image_peak + source + 3 * geometry_bytes))


def allocate_recruits(budget_gib, per_recruit_bytes, *, resources=None,
                      max_recruits=MAX_RECRUITS, plane_count=None):
    resources = resources or system_resources()
    if not math.isfinite(budget_gib) or budget_gib <= 0:
        raise ValueError("RAM budget must be a positive finite number.")
    if not isinstance(max_recruits, int) or not 1 <= max_recruits <= MAX_RECRUITS:
        raise ValueError(f"Maximum recruits must be between 1 and {MAX_RECRUITS}.")
    if per_recruit_bytes <= 0:
        raise ValueError("Per-recruit memory must be positive.")
    reserve = reserve_bytes(resources)
    usable = min(budget_gib * GIB, max(0, resources.available - reserve))
    count = min(max_recruits, resources.cpus, int(usable // per_recruit_bytes))
    if plane_count is not None:
        count = min(count, max(0, plane_count))
    return Allocation(count, budget_gib, usable / GIB,
                      per_recruit_bytes / GIB, reserve / GIB)
