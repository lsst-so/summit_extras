# Testing Guide

## Test Framework

Tests use `pytest` via the `lsst.utils.tests` framework, which provides
memory leak detection and standard LSST test infrastructure.

## Running Tests

**The DM Stack environment must be sourced before running any tests.** See the
main CLAUDE.md for details.

### Run all tests

```bash
source ~/stack.sh && . ~/setup_packages.sh && pytest tests/
```

### Run a single test file

```bash
source ~/stack.sh && . ~/setup_packages.sh && pytest tests/test_annotations.py -v
```

### Run with output

```bash
source ~/stack.sh && . ~/setup_packages.sh && pytest tests/ -v -s
```

## Test Files and Coverage

| Test File | Module Under Test | What It Tests | Requirements |
|-----------|-------------------|---------------|-------------|
| `test_annotations.py` | `annotations.py` | Tag/note retrieval, existence checks, dataId lookup. Uses mock pickle data created in setUp. | None (self-contained) |
| `test_animation.py` | `animation.py` | Animator creation and basic operation. | Skips if `ffmpeg` is not installed or LATISS butler repo is not available. |
| `test_focusAnalysis.py` | `focusAnalysis.py` | SpectralFocusAnalyzer and NonSpectralFocusAnalyzer against real data. | Requires LATISS butler access with specific dayObs/seqNum data. Skips in CI. |
| `test_logUtils.py` | `logUtils.py` | Import smoke test for LogBrowser. | None (import-only). |

## Test Environment Notes

- **Butler access**: Several tests require access to a real butler repository
  with LATISS data. These tests are decorated with `@unittest.skipUnless` and
  will skip gracefully in environments without data access.
- **ffmpeg**: The animation test requires `ffmpeg` to be on PATH.
- **CI**: GitHub Actions runs linting and formatting checks only (see
  `.github/workflows/`). The unit tests that require butler data are intended
  to be run locally or on systems with data access (summit, USDF).

## Coverage Gaps

The following modules have no dedicated tests:
- `assessQFM.py` — requires baseline parquet data and butler
- `fastStarTrackerAnalysis.py` — requires star tracker images
- `headerFunctions.py` — requires FITS files
- `imageSorter.py` — interactive, requires display
- `monitoring.py` — requires live butler repo and Firefly
- `ringssSeeing.py` — requires EFD access
- `slewTimingAuxTel.py` / `slewTimingSimonyi.py` — require EFD access
- `plotting/*` — all plotting modules (require butler + specific data)

Most of these are difficult to test in isolation because they depend on real
observatory data or external services.

## Linting (CI)

The GitHub Actions workflows run these checks on every PR:

- **lint.yaml** — flake8 via shared `lsst/rubin_workflows`
- **formatting.yaml** — black + isort via shared `lsst/rubin_workflows`
- **do_not_merge.yaml** — blocks PRs with "DO NOT MERGE" in commits
- **rebase_checker.yaml** — validates clean rebase state

### Running linters locally

```bash
source ~/stack.sh && . ~/setup_packages.sh && python -m flake8 python/
source ~/stack.sh && . ~/setup_packages.sh && python -m black --check python/
source ~/stack.sh && . ~/setup_packages.sh && python -m isort --check python/
source ~/stack.sh && . ~/setup_packages.sh && python -m mypy python/ tests/
```

### Auto-formatting

```bash
source ~/stack.sh && . ~/setup_packages.sh && python -m black python/ tests/
source ~/stack.sh && . ~/setup_packages.sh && python -m isort python/ tests/
```
