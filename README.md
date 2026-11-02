# UT26 CFD Post-Processing Pipeline
**Author:** Julian G-A  
**Team:** University of Toronto Formula Racing (UTFR)  
**Vehicle:** UT26 Formula SAE

A desktop application and Python pipeline for post-processing STAR-CCM+ aerodynamic simulation data exported in EnSight Gold format.

---

## Features

### Slice Pipeline
Interpolates volumetric CFD scalar fields onto regular 2D grids across X, Y, and Z planes. Each plane is masked against the geometry to correctly handle solid regions, multi-element wings, and hollow surfaces using a CFD-point-veto flood fill. Outputs PNG images and NPZ data files.

**Scalars supported:** Cp, CpT, Vi, Vj, Vk (and any custom array via auto-detection)

### Recruits (parallel plane processing)
The **Run** tab includes a **Recruits & RAM Budget** slider. Each recruit is an
independent CPU process that masks, interpolates and exports one plane at a time.
The estimated count updates with the RAM budget, resolution, selected fields,
axes, figure size and available system memory. Up to eight recruits run at once,
also limited by CPU availability and the number of planes. Uncheck **Use multiple
recruits** to process one plane at a time.

The slider maximum follows currently available RAM minus the system reserve,
rounded down to 0.25 GiB, and refreshes every two seconds while idle. If available
RAM drops below the selected budget, the selection is reduced to match. The CPU
limit is the number of logical CPUs this process is allowed to use, capped at
eight; it does not measure current CPU utilization.

Use **Analyze RAM from input files** in either Run or Compare to calculate RAM per
recruit from the actual inputs. Run reads the geometry size and largest selected
plane; Compare reads saved NPZ grid headers. The GUI and execution use the same
calculation. Until analysis finishes, the GUI asks for analysis instead of showing
a generic memory estimate. Changing inputs invalidates the analyzed estimate;
changing resolution updates the Run estimate from the analyzed input sizes.

The pipeline recalculates the actual count after loading and staging the data,
using the largest plane and geometry size. It reserves the greater of 2 GiB or
15% of total RAM for the system. The budget covers estimated recruit memory;
it is **not a hard memory cap**, and the parent process needs additional memory
to load the EnSight files. If even one recruit does not fit, the run stops with
an actionable error instead of exceeding the estimate deliberately.

Per-plane data is staged in the system temporary directory, so adequate temporary
disk space is required. Recruits read only their plane data and each load one copy
of the geometry. Pending jobs are bounded by the recruit count; temporary files
are cleaned up on success or failure. Fields share interpolation triangulations
within each fluid region. Existing PNG and NPZ formats are retained. The Log tab
reports the actual recruit count, completed/skipped planes and wall-clock time.

Compare has the same recruit controls. Each recruit compares one matching plane
pair across the selected scalars. The parent writes the ordered summary CSV;
the optional aero report runs after plane processing. CLI comparisons accept
`--ram-gb` and `--recruits` too. Aero-only comparisons (`--no-images`) use no recruits.

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

# Let an estimated 4 GiB budget select up to 8 recruits
python cli.py run --ram-gb 4 --recruits 8

# Process one plane at a time
python cli.py run --ram-gb 4 --recruits 1

# Compare two runs
python cli.py compare DIR_A DIR_B --label-a "Baseline" --label-b "Winglet"

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

When calling `run_pipeline()` from another Python script, place the call inside
an `if __name__ == '__main__':` guard because recruits use spawned processes.
`ram_budget_gb`, `max_recruits`, and an optional `progress_callback` are available
through the Python API. Callbacks run on the calling thread.

---

*Do not distribute without permission.*  
*Julian G-A *
