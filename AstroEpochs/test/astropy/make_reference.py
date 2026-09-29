"""Reference values for AstroEpochs from Astropy. Not run by the test suite.

Writes reference.csv next to this file: for each case, an input epoch in one scale and format, and
what Astropy gives for it in every scale, as a two-part Julian date and as an ISOT string at
millisecond precision. test_correctness_astropy_benchmark.jl reads the file, so the suite never
runs Python.

Run with the Astropy harness (C:\\Users\\steve\\Dev\\TestHarnesses\\astropy):
    uv run --project C:\\Users\\steve\\Dev\\TestHarnesses\\astropy python test/astropy/make_reference.py

Cases, by `kind`:
  random      uniform epochs 1972-2100, as a JD in each scale and as a UTC ISOT string
  leap        either side of every leap second, in UTC and in TAI
  leap_inst   instants inside a leap second (UTC 23:59:60.x); AstroEpochs cannot represent these
              yet, so the benchmark reports them as a known gap
  rounding    epochs whose ISOT string rounds up across a second, minute, hour or day
  special     J2000, the 1977 TCG/TCB epoch, MJD 0, and the ends of the range
  pre1972     UTC before 1972-01-01, where Astropy applies the drifting pre-1972 offsets and
              AstroEpochs does not; a known gap by decision
"""
import csv
import math
import random
from pathlib import Path

import astropy
import erfa
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False          # the bundled leap-second table; no network, reproducible

SCALES = ("tai", "tt", "utc", "tdb", "tcb", "tcg")
HERE = Path(__file__).parent
random.seed(20260929)

cases = []


def add(kind, scale, fmt, value, value2=None):
    cases.append((kind, scale, fmt, value, value2))


# ── random ──────────────────────────────────────────────────────────────────
jd_lo, jd_hi = 2441317.5, 2488069.5                  # 1972-01-01 to 2100-01-01
for scale in ("tai", "tt", "tdb", "tcb", "tcg"):
    for _ in range(200):
        jd1 = float(random.randint(int(jd_lo), int(jd_hi)))
        jd2 = random.uniform(-0.5, 0.5)
        add("random", scale, "jd", jd1, jd2)
for _ in range(200):
    t = Time(random.uniform(jd_lo, jd_hi), format="jd", scale="utc")
    add("random", "utc", "isot", t.isot)

# ── leap seconds ────────────────────────────────────────────────────────────
leap = erfa.leap_seconds.get()
for row in leap[1:]:                                   # each date on which TAI-UTC stepped
    y, m = int(row["year"]), int(row["month"])
    if y < 1972 or (y, m) == (1972, 1):
        continue
    new = Time(f"{y:04d}-{m:02d}-01T00:00:00", scale="utc")
    prev = (new - 1).utc.isot[:10]                     # the day that ends with the leap second
    for hms in ("12:00:00.000", "23:59:59.000", "23:59:59.999"):
        add("leap", "utc", "isot", f"{prev}T{hms}")
    for hms in ("60.000", "60.500", "60.999"):
        add("leap_inst", "utc", "isot", f"{prev}T23:59:{hms}")
    for hms in ("00:00:00.000", "00:00:00.001"):
        add("leap", "utc", "isot", f"{new.isot[:10]}T{hms}")
    # The same boundary approached in TAI, which is continuous across it.
    tai = new.tai
    for ds in (-1.5, -1.0 - 1e-6, -1.0 + 1e-6, -0.5, -1e-6, 0.0, 1e-6, 1e-3):
        t = tai + ds / 86400.0
        kind = "leap_inst" if -1.0 < ds < 0.0 else "leap"
        add(kind, "tai", "jd", t.jd1, t.jd2)

# ── rounding ────────────────────────────────────────────────────────────────
for s in ("2024-06-01T23:59:59.9996", "2024-06-01T23:59:59.9994", "2024-06-01T12:59:59.9995",
          "2024-06-01T00:00:59.9999", "2023-12-31T23:59:59.99951", "2016-12-31T23:59:59.9996"):
    add("rounding", "utc", "isot", s)
for scale in ("tt", "tai", "tdb"):
    for frac_s in (-0.0004, -0.0006, 0.0004):          # seconds either side of midnight
        add("rounding", scale, "jd", 2460463.0, -0.5 + frac_s / 86400.0 + (1.0 if frac_s < 0 else 0.0))

# ── special epochs ──────────────────────────────────────────────────────────
add("special", "tt", "jd", 2451545.0, 0.0)                              # J2000
add("special", "tai", "isot", "1977-01-01T00:00:00.000")               # TCG/TCB epoch
add("special", "tt", "isot", "1977-01-01T00:00:32.184")
add("special", "tcg", "isot", "1977-01-01T00:00:32.184")
add("special", "tcb", "isot", "1977-01-01T00:00:32.184")
add("special", "tdb", "isot", "1977-01-01T00:00:32.184")
add("special", "tt", "jd", 2400000.0, 0.5)                              # MJD 0
add("special", "utc", "isot", "1972-01-01T00:00:00.000")               # first integer offset
for y in (1600, 1700, 1900, 2200, 2500):
    add("special", "tt", "isot", f"{y}-03-01T06:00:00.000")
    add("special", "tdb", "isot", f"{y}-09-15T18:30:00.000")

# ── pre-1972 UTC ────────────────────────────────────────────────────────────
for s in ("1960-01-01T00:00:00.000", "1965-06-15T12:00:00.000", "1971-12-31T23:59:59.000"):
    add("pre1972", "utc", "isot", s)

# ── evaluate ────────────────────────────────────────────────────────────────
rows = []
for kind, scale, fmt, value, value2 in cases:
    t = Time(value, value2, format="jd", scale=scale) if fmt == "jd" else \
        Time(value, format="isot", scale=scale, precision=3)
    t.precision = 3
    row = {"kind": kind, "in_scale": scale, "in_format": fmt,
           "in_value": repr(value) if fmt == "jd" else value,
           "in_value2": repr(value2) if fmt == "jd" else ""}
    for s in SCALES:
        ts = getattr(t, s)
        row[f"{s}_jd1"] = repr(float(ts.jd1))
        row[f"{s}_jd2"] = repr(float(ts.jd2))
        row[f"{s}_isot"] = ts.isot
    rows.append(row)

out = HERE / "reference.csv"
with out.open("w", newline="") as fh:
    fh.write(f"# Astropy {astropy.__version__}, pyerfa {erfa.__version__} "
             f"(ERFA {erfa.version.erfa_version}); generated by make_reference.py\n")
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"wrote {len(rows)} cases to {out}")
