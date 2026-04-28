Run mypy type checking on summit_extras. Source the DM Stack environment first.

Run:

```
source ~/stack.sh && . ~/setup_packages.sh && python -m mypy python/ tests/ 2>&1
```

Report:
1. Total number of errors
2. Group errors by category (missing imports, type mismatches, etc.)
3. For any new or surprising errors, suggest fixes

Note: mypy.ini has extensive ignore rules for LSST and third-party packages. Errors from those packages can usually be addressed by adding ignore rules to mypy.ini rather than fixing the upstream code.
