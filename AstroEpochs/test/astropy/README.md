# Astropy and ERFA reference values

The scripts here produced the truth for two test files. The test suite does not run them and does
not need Python.

| Script | Writes | Used by |
|---|---|---|
| `make_reference.py` | `reference.csv`: 1793 epochs, each with Astropy's two-part JD and ISOT string in every scale | `test_correctness_astropy_benchmark.jl` |
| `make_erfa_parity.py` | printed Julia literals: pyerfa's output for each ported ERFA routine at its edge cases | `test_correctness_erfa_parity.jl` |

Rerun them when AstroEpochs' conversions change on purpose, or to benchmark against a newer Astropy.
They need Python with the versions below installed, and run from the AstroEpochs folder:

```
python test/astropy/make_reference.py
python test/astropy/make_erfa_parity.py
```

`make_reference.py` turns off Astropy's automatic IERS download, so the leap-second table is the
one bundled with the installed `astropy-iers-data`, and a rerun is reproducible. Versions used:
Astropy 8.0.1, pyerfa 2.0.1.5 (ERFA 2.0.1), astropy-iers-data 0.2026.9.28.
