Run all linters and formatters for summit_extras. Source the DM Stack environment first.

Run these in sequence, reporting results for each:

1. **black** (formatting check):
   ```
   source ~/stack.sh && . ~/setup_packages.sh && python -m black --check python/ tests/ 2>&1
   ```

2. **isort** (import ordering):
   ```
   source ~/stack.sh && . ~/setup_packages.sh && python -m isort --check python/ tests/ 2>&1
   ```

3. **flake8** (linting):
   ```
   source ~/stack.sh && . ~/setup_packages.sh && python -m flake8 python/ 2>&1
   ```

If the user says "fix" or "autoformat", run black and isort WITHOUT --check to auto-fix, then re-run flake8 to report remaining issues.

Summarize: how many files need formatting, how many lint errors, and the most common error codes.
