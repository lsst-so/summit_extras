# Package Overview

`lsst.summit.extras` is a collection of analysis tools, interactive utilities,
and visualization code used for Rubin Observatory summit operations and
commissioning work. It sits in the `lsst.summit` namespace alongside
`summit_utils` (which provides the lower-level data access and utility
primitives that this package builds on).

## Module Map

### Analysis Modules

| Module | Key Classes / Functions | Purpose |
|--------|------------------------|---------|
| `focusAnalysis.py` | `SpectralFocusAnalyzer`, `NonSpectralFocusAnalyzer`, `FitResult` | Analyze focus sweeps — fits Gaussians across spectral slices or PSF sizes, then fits a parabola to find best focus. |
| `assessQFM.py` | `AssessQFM` | Compare QuickFrameMeasurementTask results against a shipped baseline dataset (`data/qfm_baseline_assessment.parq`). Uses multiprocessing for parallel evaluation. |
| `fastStarTrackerAnalysis.py` | `Source`, `NanSource`, `findFastStarTrackerImageSources()`, `plotSourceMovement()` | Detect and characterize sources in fast star tracker camera images. Provides source centroid, flux, cutout extraction, and motion tracking plots. |
| `ringssSeeing.py` | `SeeingConditions` | Retrieve seeing measurements from the RINGSS turbulence profiler via the EFD. Returns a dataclass with seeing, ground layer, free atmosphere, and FWHM values. |
| `slewTimingAuxTel.py` | `plotExposureTiming()` | Plot exposure and slew timing for AuxTel by querying ATPtg/ATMCS command sequences from the EFD. |
| `slewTimingSimonyi.py` | `plotExposureTiming()` | Same as above but for the Simonyi Survey Telescope (MTPtg/MTMount). |
| `logUtils.py` | `LogBrowser` | Browse butler task logs to identify failure modes. Groups errors by message, with special-case handling for known variable-text errors. |

### Interactive / Notebook Modules

| Module | Key Classes / Functions | Purpose |
|--------|------------------------|---------|
| `animation.py` | `Animator` | Create video animations (MP4/GIF) from a list of butler dataIds. Supports custom remapping, centroid overlays, and ffmpeg-based encoding. |
| `imageSorter.py` | `ImageSorter` | Interactive image tagging UI for Jupyter notebooks. Presents images for annotation with single-key tags (Q=bad star, F=focus, D=donut, etc.) and persists results to pickle files. |
| `annotations.py` | `Annotations` | Read-only interface to annotation data written by `ImageSorter`. Query by dataId to get tags, notes, and check whether images have been examined. |
| `monitoring.py` | `Monitor` | Real-time AuxTel image monitoring via bestEffortIsr + Firefly display. **Largely superseded by RubinTV** but still importable. |

### Utility Modules

| Module | Key Classes / Functions | Purpose |
|--------|------------------------|---------|
| `headerFunctions.py` | `compareHeaders()`, `keyValuesSetFromFiles()`, `buildHashAndHeaderDicts()` | FITS header inspection: compare headers between files, extract specific key values across file sets, and cache header libraries as pickle files. |

### Plotting Subpackage (`plotting/`)

The `plotting/` directory is a namespace subpackage (no `__init__.py`).
Each module is imported directly:

```python
from lsst.summit.extras.plotting.psfPlotting import makeFocalPlanePlot
```

See [plotting.md](plotting.md) for detailed documentation.

## CLI Entry Points (`bin.src/`)

These scripts are installed to `$PATH` by the SCons build:

| Script | Usage | Purpose |
|--------|-------|---------|
| `fitsheaderCompare.py` | `fitsheaderCompare.py file1.fits file2.fits` | Compare FITS headers between two files, highlighting differences. |
| `imageSorter.py` | `imageSorter.py "glob/pattern/*.fits" output.pkl` | Run the interactive image annotation tool on a set of FITS files. |
| `scrapeHeaders.py` | `scrapeHeaders.py "glob/*.fits" [-k KEY] [-j KEY1 KEY2]` | Extract FITS header values. Supports per-file output, recursive walks, and library caching. |

## Data Flow

Most modules in this package follow a common pattern:

1. **Input**: Butler data access (exposures, catalogs, visit summaries) and/or
   EFD telemetry queries.
2. **Processing**: Domain-specific analysis (fitting, measurement, filtering).
3. **Output**: Matplotlib figures, text summaries, or persisted annotation
   files.

The package does not write back to the butler or modify telescope state — it is
purely analytical and visual.

## Dependency Layers

```
  ┌──────────────────────────────┐
  │     summit_extras            │  ← This package
  │  (analysis, plotting, tools) │
  └──────────┬───────────────────┘
             │ imports
  ┌──────────▼───────────────────┐
  │     summit_utils             │  ← Utility primitives (bestEffortIsr,
  │  (data access, EFD, etc.)    │     EFD client wrappers, camera utils)
  └──────────┬───────────────────┘
             │ imports
  ┌──────────▼───────────────────┐
  │     LSST Science Pipelines   │  ← afw, daf_butler, ip_isr,
  │  (core DM stack)             │     meas_algorithms, pipe_tasks, etc.
  └──────────┬───────────────────┘
             │ imports
  ┌──────────▼───────────────────┐
  │  Third-party (numpy, scipy,  │
  │  matplotlib, astropy, etc.)  │
  └──────────────────────────────┘
```
