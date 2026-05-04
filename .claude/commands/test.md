Run the test suite for summit_extras. Source the DM Stack environment first — this is mandatory.

Run:

```
source ~/stack.sh && . ~/setup_packages.sh && pytest tests/ -v --tb=short 2>&1
```

Report:
1. How many tests passed, failed, skipped
2. For any failures, show the traceback and diagnose the root cause
3. For skipped tests, note why they were skipped (missing data, missing ffmpeg, etc.)

If the user provided arguments like a specific test file or -k filter, pass those through to pytest instead of running the full suite.
