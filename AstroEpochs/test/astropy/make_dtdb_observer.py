"""pyerfa's TDB - TT with observer terms, at random dates and observers. Not run by the test suite.

Writes dtdb_observer.csv, read by test_correctness_arithmetic_precision.jl. The core conversions
evaluate dtdb at the geocentre, so nothing else in the suite holds the observer terms.

Run from the AstroEpochs folder, with pyerfa 2.0.1.5 installed:
    python test/astropy/make_dtdb_observer.py
"""
import csv
import math
import random

import erfa

random.seed(20261102)
with open("test/astropy/dtdb_observer.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["date1", "date2", "ut", "elong_rad", "u_km", "v_km", "tdb_minus_tt_s"])
    for _ in range(100):
        date1 = float(math.floor(random.uniform(2378496.5, 2524593.5)))
        date2, ut = random.random(), random.random()
        elong = random.uniform(-math.pi, math.pi)
        u, v = random.uniform(0.0, 6400.0), random.uniform(-6400.0, 6400.0)
        w.writerow([repr(x) for x in (date1, date2, ut, elong, u, v, float(erfa.dtdb(date1, date2, ut, elong, u, v)))])
