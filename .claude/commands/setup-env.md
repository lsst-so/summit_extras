Verify the DM Stack environment is working for this package. Run:

```
source ~/stack.sh && . ~/setup_packages.sh && python -c "from lsst.summit.extras import Animator; print('summit_extras imports OK')"
```

If it succeeds, report the result. If it fails, diagnose the error — common issues:
- `~/stack.sh` or `~/setup_packages.sh` missing or broken
- EUPS package not set up (check `eups list summit_extras`)
- Missing dependency (check the import traceback)

Then run a quick lint check:

```
source ~/stack.sh && . ~/setup_packages.sh && python -m flake8 python/ --count --statistics
```

Report a summary of environment health and any lint issues found.
