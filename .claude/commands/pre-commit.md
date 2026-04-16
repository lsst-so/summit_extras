Run the full pre-commit check suite, simulating what would run before a git commit.

```
source ~/stack.sh && . ~/setup_packages.sh && pre-commit run --all-files 2>&1
```

If pre-commit is not installed, fall back to running the individual hooks manually:

1. Trailing whitespace check
2. YAML validation
3. isort --check
4. black --check
5. flake8

Report all failures and offer to auto-fix formatting issues (black, isort, trailing whitespace).
