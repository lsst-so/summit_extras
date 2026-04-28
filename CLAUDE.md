# summit_extras - Agent Guide

`lsst.summit.extras` is a **library** of analysis tools, plotting utilities,
and interactive helpers for Vera Rubin Observatory summit operations. It
provides focus analysis, slew timing, seeing measurements, image annotation,
FITS header inspection, and focal-plane visualizations used by observers and
commissioning scientists at the summit and in Jupyter notebooks.

## This is a library, not an application

Unlike `rubintv_production`, this package **is** imported by other code and
notebooks. Its public surface matters:
- Renaming or removing a public symbol is a breaking change for downstream
  notebook users. Prefer deprecation when altering widely-used APIs.
- The `__init__.py` re-exports from several modules via `from .X import *`.
  Keep this in sync when adding new public classes or functions.
- The `plotting/` subpackage has no `__init__.py` — modules are imported
  directly (e.g. `from lsst.summit.extras.plotting.psfPlotting import ...`).

## Architecture Documentation

Detailed architecture docs live in `architecture/` and MUST be kept up to date
whenever architectural changes are made:

- [Package Overview](architecture/overview.md) - module map, what lives where,
  and how the pieces fit together
- [Plotting Subsystem](architecture/plotting.md) - the plotting subpackage,
  its modules, and their data sources
- [Testing Guide](architecture/testing.md) - how to run tests, what's covered,
  and environment requirements

**When modifying code that changes any of the following, update the
corresponding architecture doc(s) in the same change:**
- Adding or removing a module or public class (`overview.md`)
- Changing the plotting subpackage layout (`plotting.md`)
- Changing test infrastructure or dependencies (`testing.md`)

## Quick Orientation

```
python/lsst/summit/extras/          # Main Python package
  __init__.py                        # Wildcard re-exports (animation, focus, headers, etc.)
  animation.py                       # Animator - video/GIF creation from dataIds
  annotations.py                     # Annotations - read image tags from ImageSorter
  assessQFM.py                       # AssessQFM - compare QuickFrameMeasurement vs baseline
  fastStarTrackerAnalysis.py         # Fast star tracker image source detection & plotting
  focusAnalysis.py                   # SpectralFocusAnalyzer, NonSpectralFocusAnalyzer
  headerFunctions.py                 # FITS header comparison, extraction, library caching
  imageSorter.py                     # ImageSorter - interactive image tagging/annotation
  logUtils.py                        # LogBrowser - failure zoology from butler task logs
  monitoring.py                      # Monitor - AuxTel real-time display (largely superseded)
  ringssSeeing.py                    # SeeingConditions - RINGSS seeing data via EFD
  slewTimingAuxTel.py                # plotExposureTiming() for AuxTel
  slewTimingSimonyi.py               # plotExposureTiming() for Main Telescope
  plotting/                          # Plotting subpackage (no __init__.py)
    focusSweep.py                    # Focus sweep data collection & parabola fitting
    fwhmFocalPlane.py                # FWHM focal plane map from visit summaries
    psfPlotting.py                   # PSF/source catalog visualization on focal plane
    zernikePredictedFwhm.py          # Zernike mode & DOF FWHM prediction plots

bin.src/                             # CLI entry points (installed to bin/ by SCons)
  fitsheaderCompare.py               # Compare FITS headers between two files
  imageSorter.py                     # Interactive image annotation from glob pattern
  scrapeHeaders.py                   # Extract FITS header values from files

data/                                # Shipped data files
  qfm_baseline_assessment.parq       # Baseline QuickFrameMeasurement results

tests/                               # Unit tests (pytest via lsst.utils.tests)
  test_animation.py                  # Animator (skips without ffmpeg/LATISS)
  test_annotations.py                # Annotations with mock data
  test_focusAnalysis.py              # Focus analyzers (requires LATISS butler)
  test_logUtils.py                   # LogBrowser import smoke test
```

## Key Concepts

- **Butler**: LSST's data access layer for reading/writing datasets.
- **dayObs / seqNum**: Integer observing date (YYYYMMDD) and sequence number
  identifying an exposure within a night.
- **EFD**: Engineering Facilities Database — time-series telemetry from the
  telescope, accessed via `lsst_efd_client`.
- **QuickFrameMeasurement (QFM)**: A fast pipeline task that measures basic
  image quality metrics.
- **RINGSS**: A turbulence profiler instrument that measures seeing conditions.
- **ISR**: Instrument Signature Removal — the first processing step that
  removes detector artifacts.

## Development

### Naming

- camelCase for all variables, functions, methods, and attributes.
- PascalCase for classes.
- No snake_case except when required by external APIs.
- All function/method names must contain a verb (including private ones),
  e.g. ``getTrackingKey`` not ``trackingKey``. Exception: ``fromX``
  class methods for constructors do not need a verb.
- Prefer longer, descriptive variable names over short, abbreviated ones.
  Established abbreviations are fine (``dets`` for ``detectors``, ``expId``
  for exposure identifier, etc.).

### Formatting

- black (line-length 110), isort (black profile)

### Type Annotations

- Use built-in types (``int``, ``str``, ``float``, ``dict``, ``list``,
  ``tuple``, etc.).
- Never import ``Dict``, ``List``, ``Tuple``, ``Optional``, or ``Union``
  from ``typing``.
- Use ``| None`` instead of ``Optional[...]``.

### Docstrings

- Use numpydoc format.
- Include types for every parameter and for the return value (if not
  ``None``).
- Always name the return value unless the return type is ``None`` (omit
  the Returns section in that case).
- If a parameter is ``| None``, describe its type as ``<type>, optional``.
- Argument order in docstrings must match the function signature, and
  types must be correct.
- No docstrings for class ``__init__``; document the class instead.

Example:

```python
def myFunction(param1: int, param2: str | None = None) -> bool:
    """This function does something.

    Parameters
    ----------
    param1 : `int`
        The first parameter.
    param2 : `str`, optional
        The second parameter.

    Returns
    -------
    result : `bool`
        The result of the function.
    """
    return param1 > 0 and param2 != "hello"
```

### Environment Setup

**ALWAYS run code from this package, run its tests, or import any of its
modules through the DM Stack environment.** The package lives in the
``lsst.*`` namespace and imports many sibling packages from the LSST stack
(``lsst.daf.butler``, ``lsst.afw``, ``lsst.summit.utils``, etc.).
Running ``python``, ``pytest`` or ``python -c "import lsst.summit..."``
without the stack sourced will fail with ``ModuleNotFoundError`` every
time. Do not waste a turn discovering this — source the stack first.

**Every time** you need to run, test, or import package code, prefix the
command so that both lines below are sourced first, in this order:

```bash
source ~/stack.sh && . ~/setup_packages.sh && <your command>
```

For example:

```bash
source ~/stack.sh && . ~/setup_packages.sh && pytest tests/test_annotations.py
source ~/stack.sh && . ~/setup_packages.sh && python -c "from lsst.summit.extras.focusAnalysis import SpectralFocusAnalyzer; print('ok')"
```

Do **not** try ``PYTHONPATH=python python ...`` or any other workaround —
sourcing the stack is the only supported way and is required for every
single invocation, because the shell state from previous Bash tool calls
does not persist.

### Linting & Tooling

- **Linting**: flake8 (max-line-length 110, max-doc-length 79)
- **Type checking**: mypy (Python 3.11 target, pydantic plugin)
- **Pre-commit hooks**: trailing whitespace, YAML, isort, black, flake8
- **Build system**: LSST SCons + pyproject.toml (Python-only, no C++)
- **License**: GPLv3

### Running Linters

```bash
source ~/stack.sh && . ~/setup_packages.sh && python -m flake8 python/
source ~/stack.sh && . ~/setup_packages.sh && python -m mypy python/ tests/
source ~/stack.sh && . ~/setup_packages.sh && python -m black --check python/
source ~/stack.sh && . ~/setup_packages.sh && python -m isort --check python/
```

### EUPS Dependencies

This package depends on these LSST stack packages (declared in
`ups/summit_extras.table`): afw, atmospec, base, daf_butler, geom, ip_isr,
meas_algorithms, pipe_tasks, utils, pex_exceptions, summit_utils,
display_matplotlib.
