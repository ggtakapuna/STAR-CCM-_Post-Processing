# UT27 CFD Post-Processing Pipeline
**Author:** Julian G-A  
**Contributors:** John Baton  
**Team:** University of Toronto Formula Racing (UTFR)  
**Vehicle:** UT27 Formula SAE

A desktop application and Python pipeline for post-processing STAR-CCM+ aerodynamic simulation data exported in EnSight Gold format.

---

## Features

### Slice Pipeline
Interpolates volumetric CFD scalar fields onto regular 2D grids across X, Y, and Z planes. Each plane is masked against the geometry to correctly handle solid regions, multi-element wings, and hollow surfaces using a CFD-point-veto flood fill. Outputs PNG images and NPZ data files.

**Scalars supported:** Cp, CpT, Vi, Vj, Vk (and any custom array via auto-detection)

### Recruits (parallel processing)
The **Run** and **Compare** tabs include **Recruits & RAM Budget** controls.
Each recruit is an independent CPU process handling one plane, or one matching
plane pair for comparison, across all selected scalars. Up to eight recruits run
at once, limited by RAM, the logical CPUs available to this process, and the
number of planes. Uncheck **Use multiple recruits** to process one at a time.

The RAM slider maximum is currently available RAM minus the greater of 2 GiB or
15% of total RAM, rounded down to 0.25 GiB. It refreshes every two seconds while
idle and reduces the selection if available RAM falls. The CPU limit uses process
CPU affinity, with the system logical CPU count as fallback, capped at eight;
it does not measure current CPU utilization.

Click **Analyze RAM from input files** to calculate memory per recruit from the
actual inputs. Run reads the geometry size and largest selected plane; Compare
reads saved NPZ grid headers. The GUI and execution use the same calculation.
Until analysis finishes, the GUI asks for analysis instead of showing a generic
memory estimate. Changing inputs invalidates the analysis; changing Run resolution
updates its estimate using the analyzed input sizes. Execution reports the actual
recruit count and memory estimate in the GUI and Log tab.

RAM is an estimated processing budget, **not a hard memory cap**. The parent
process also needs memory to load EnSight files. Allocation is checked again
after loading; if even one recruit does not fit, the run stops with an error.

Run stages per-plane data in the system temporary directory, requiring adequate
temporary disk space. Recruits read their own plane data and each load one copy of
the geometry. Pending jobs are bounded by the recruit count; temporary files are
cleaned up after success or failure. Fields share interpolation triangulations
within each fluid region. Existing PNG and NPZ formats are retained.

Compare recruits produce delta PNGs and statistics for matching plane pairs.
The parent writes the summary CSV in a consistent order and then generates the
optional aero report. Aero-only CLI comparisons use no recruits.

### Comparison Pipeline
Computes plane-by-plane deltas between two simulation runs. Delta images use a three-colour geometry overlay:
- **Gray** — geometry present in both runs
- **Magenta** — new geometry (Run B only)
- **Brown** — removed geometry (Run A only)

Optionally parses STAR-CCM+ `.log` files to generate an aerodynamic balance delta report.

### Surface Rendering
Maps surface CFD scalar data (Cp, Cf) directly onto the geometry mesh via KD-tree cell matching. Produces:
- Full-car renders across 8 views (bottom, top, iso front/rear, side L/R, front, rear)
- Per-component zoomed renders with auto-fit framing
- Front wing + EPs and Rear wing + EPs merged into single parts

### Interactive 3D Export
Exports a self-contained `.html` file using Plotly. Opens in any browser — no software installation required. Supports rotation, zoom, and per-part toggling via the legend.

### 3D File Export
- **VTP / VTK** — preserves scalar data, opens in ParaView
- **STL / PLY / OBJ** — geometry only, opens in Fusion360, Meshmixer, Windows 3D Viewer

### Visualization Viewer
A standalone Tkinter PNG browser (`viz_viewer.py`) for scrolling through pipeline and compare output folders. Pairs parent, new, and delta slices on (axis, plane index) rather than filename, so mismatched stems still line up. View one set at a time or all three side by side, switch scalar/axis, and navigate with arrow keys, mouse wheel, or PgUp/PgDn. Requires Pillow.

---

## Colormaps

### Cp
Physics-informed diverging scale. Grey is pinned exactly at Cp = 0 regardless of the data range.

| Range | Colour |
|---|---|
| Most negative (suction) | Deep blue |
| Mid suction | Cyan → green |
| Zero | Neutral grey |
| Positive (pressure) | Orange → deep red |

### Cf
Sequential scale with a hard sentinel colour at zero.

| Value | Colour |
|---|---|
| Cf = 0 (separation / stagnation) | Pink (hard step, not blended) |
| Low Cf | Light blue |
| High Cf (attached flow) | Deep blue |

Range: 0 – 0.03

---

## Installation

```bash
pip install pyvista numpy scipy matplotlib plotly pillow scikit-image psutil
```

Tkinter is included with most Python distributions. If missing:
```bash
# Ubuntu/Debian
sudo apt install python3-tk
```

---

## File Structure

```
cfd_slicer/
├── config.py          # All user settings — paths, scalars, colormaps, bounds
├── geometry.py        # Geometry masking (CFD-point-veto flood fill)
├── interpolation.py   # Per-region griddata interpolation
├── plotting.py        # Slice and delta plot generation
├── io_utils.py        # Case file loading, NPZ save/load, plane discovery
├── pipeline.py        # Slice pipeline orchestration
├── compare.py         # Run comparison and delta rendering
├── aero_report.py     # Log file parser and aero balance report
├── surface_render.py  # Surface rendering, part views, HTML and 3D export
├── cli.py             # Command-line interface
├── gui.py             # Desktop GUI (Tkinter)
└── viz_viewer.py      # Standalone PNG browser for output folders
```

---

## Output Structure

```
output/
├── PNG/
│   ├── Cp/
│   │   ├── X/   slice_X_001_+100.0mm.png
│   │   ├── Y/
│   │   └── Z/
│   └── CpT/
├── NPZ/
│   ├── X/       slice_X_001_+100.0mm.npz
│   ├── Y/
│   └── Z/
```

Comparison output (separate directory):
```
compare_RunB_vs_RunA/
├── PNG/
│   ├── Cp/
│   │   ├── X/   diff_Cp_X_001_+100.0mm.png
│   │   ├── Y/
│   │   └── Z/
│   └── CpT/
└── comparison_summary.csv
```

Surface renders:
```
output/
├── PNG/
│   ├── surface_Cp_bottom.png
│   ├── surface_Cp_iso_front.png
│   ├── surface_Cp_Front_wing_bottom.png   ← per-component, auto-zoomed
│   └── ...
├── NPZ/
│   └── ...
└── surface_Cp.html                         ← interactive 3D
```

---

## Usage

### GUI
```bash
python gui.py
```

### CLI
```bash
# Run slice pipeline
python cli.py run --scalars Cp CpT --axes X Y Z --res 2

# Automatically choose up to eight recruits within an estimated 4 GiB budget
python cli.py run --ram-gb 4 --recruits 8

# Process one plane at a time
python cli.py run --ram-gb 4 --recruits 1

# Compare two runs
python cli.py compare DIR_A DIR_B --label-a "Baseline" --label-b "Winglet"

# Compare plane pairs in parallel with the same RAM/CPU allocation rules
python cli.py compare DIR_A DIR_B --ram-gb 4 --recruits 8

# Inspect a saved plane
python cli.py info path/to/plane.npz

# List all saved planes
python cli.py list output/NPZ --axis Y
```

### Visualization Viewer
```bash
python viz_viewer.py
```
Point each field at an output directory containing a `PNG/` subfolder, then scroll through paired slices.

---

## Configuration

Images use a Y reflection across the X-Z plane (`Y -> -Y`) by default. This
reflects slice/compare scalar fields and geometry together, with matching axis
coordinates and Y-plane title positions. Surface PNGs and interactive HTML views
use the same orientation. Set `REFLECT_Y_IMAGES = False` in `config.py` to restore
the original display orientation. Saved NPZ data, slice filenames, comparison
statistics and mesh file exports retain the simulation coordinates. Regenerate
PNGs/HTML to apply the display change to existing results.

All settings live in `config.py`. The most commonly changed values:

```python
GEOMETRY_CASE  = r"path/to/Geometry.case"
DATA_CASE      = r"path/to/Raw_Data.case"
OUTPUT_DIR     = r"path/to/output"

RESOLUTION_MM  = 1       # grid resolution (1mm recommended)
DPI            = 150     # output image DPI
DILATION_PX    = 3       # geometry mask dilation (px)

CAR_BOUNDS_2D  = {
    0: (-0.90,  0.90, -0.10, 1.30),  # X-slice: h=Y, v=Z
    1: (-1.10,  2.20, -0.10, 1.30),  # Y-slice: h=X, v=Z
    2: (-1.10,  2.20, -0.90, 0.90),  # Z-slice: h=X, v=Y
}
```

New scalar arrays can be added in the Settings tab or via the **Detect from CFD file** button, which reads array names directly from the loaded `.case` file.

---

## Data Format

EnSight Gold (`.case` + binary data files) as exported by STAR-CCM+.

The pipeline expects three types of case files:
| File | Used for |
|---|---|
| `Geometry.case` | Surface geometry mesh (one block per part) |
| `Raw_Data_XXXXX.case` | Volumetric slice data (X/Y/Z planes) |
| `Surface_XXXXX.case` | Surface scalar data (Cp, Cf on mesh) |

---

## Dependencies

| Package | Purpose |
|---|---|
| `pyvista` | Mesh I/O, surface extraction, offscreen rendering |
| `numpy` | Grid operations, masking, NPZ I/O |
| `scipy` | KD-tree matching, griddata interpolation |
| `scikit-image` | Line rasterization for geometry masking |
| `matplotlib` | Plot generation, colormaps |
| `plotly` | Interactive 3D HTML export |
| `pillow` | Image loading for the visualization viewer |
| `psutil` | Available RAM and CPU limits for recruit allocation |
| `tkinter` | Desktop GUI |

Python 3.10+ required.

Python callers can pass `ram_budget_gb`, `max_recruits` and `progress_callback`
to `run_pipeline()` and `compare_runs()`. Callbacks receive allocation, plane and
done event dictionaries on the calling thread; GUI callers must marshal them to
the Tk thread. Place calls inside an `if __name__ == '__main__':` guard because
recruits use spawned processes, including on Windows.

Run the regression checks with `python -m unittest discover -s tests -v`.

---

*Do not distribute without permission.*  
*Julian G-A *
