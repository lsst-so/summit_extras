# Plotting Subsystem

The `lsst.summit.extras.plotting` subpackage provides specialized
visualization tools for observatory data. It is a **namespace package** with
no `__init__.py` — each module is imported directly.

## Modules

### focusSweep.py

**Purpose**: Collect focus sweep data and fit/plot focus parabolas.

**Key functions**:
- `collectSweepData(dayObs, seqNums, ...)` — Queries the consolidated
  database and EFD to build a table of hexapod positions, PSF sizes (FWHM),
  and ellipticity components (e1, e2) for a set of exposures in a focus
  sweep sequence.
- `fitSweepParabola(table, column, ...)` — Fits a parabola to the collected
  sweep data to find the best-focus position.
- `plotSweepParabola(table, fitResult, ...)` — Renders the sweep data and
  fitted parabola as a matplotlib figure.

**Data sources**: Consolidated database (via `consdb`), EFD hexapod telemetry.

### fwhmFocalPlane.py

**Purpose**: Visualize FWHM across the focal plane for a single visit.

**Key functions**:
- `getFwhmValues(butler, dataId)` — Extracts per-detector FWHM values from
  the `visitSummary` table. Converts PSF sigma to FWHM in arcseconds using
  the pixel scale.
- `makeFocalPlaneFWHMPlot(butler, dataId, ...)` — Renders a focal plane map
  colored by FWHM. Uses camera geometry to position detectors and applies a
  colormap to the FWHM values.

**Data sources**: Butler `visitSummary` tables, camera geometry.

### psfPlotting.py

**Purpose**: Visualize PSF properties and source catalogs on the focal plane
and in sky coordinates.

**Key functions**:
- `makeFocalPlanePlot(...)` — Plot per-detector values on the focal plane
  layout with optional quiver (vector field) overlays.
- `makeEquatorialPlot(...)` — Plot source properties in RA/Dec coordinates.
- `makeAzElPlot(...)` — Plot source properties in Azimuth/Elevation.
- `makeTableFromSourceCatalog(...)` — Convert a butler source catalog into a
  flat table with extended metadata (detector ID, type, position).
- Helper functions for compass roses, coordinate grids, and axis formatting.

**Data sources**: Butler source catalogs, camera geometry, WCS solutions.

### zernikePredictedFwhm.py

**Purpose**: Visualize how individual Zernike modes and degrees of freedom
(DOFs) contribute to predicted FWHM.

**Key functions**:
- `makeZernikePredictedFWHMPlot(zernikes, ...)` — Bar chart of per-Zernike-mode
  FWHM contributions.
- `makeDofPredictedFWHMPlot(dofs, ...)` — Bar chart grouping DOFs by category
  (M2 hexapod, camera hexapod, M1M3 bending modes, M2 bending modes) with
  multi-column formatting for large DOF sets.

**Data sources**: AOS (Active Optics System) pipeline outputs — Zernike
coefficient arrays and DOF vectors, typically from `lsst.ts.wep`.

**AOS support**: The `makeZernikePredictedFWHMPlot` function supports both
paired (intra/extra-focal pair) and unpaired (single donut) AOS modes. When
unpaired, certain layout adjustments apply.

## Common Patterns

All plotting modules share these patterns:

1. **matplotlib-based** — Every function produces matplotlib figures. Most
   accept an optional `fig`/`ax` parameter for embedding in larger layouts.
2. **Butler as data source** — Most functions take a butler instance and a
   dataId as their primary inputs.
3. **Camera geometry** — Focal plane plots use `lsst.afw.cameraGeom` to map
   detector IDs to physical positions.
4. **No side effects** — Plots are returned or displayed; nothing is written
   back to the butler or external storage.
