"""pyerfa's outputs for each ERFA routine AstroEpochs ports. Not run by the test suite.

Prints the Julia literals in test_correctness_erfa_parity.jl. The inputs are chosen for the edges
each routine has: half-days (where rounding direction matters), both orders of a two-part date,
leap-second days, and times that round across a day.

Run from the AstroEpochs folder, with pyerfa 2.0.1.5 installed:
    python test/astropy/make_erfa_parity.py
"""
import erfa
from astropy.utils import iers

iers.conf.auto_download = False
r = repr

print(f"# pyerfa {erfa.__version__} (ERFA {erfa.version.erfa_version})")

dates = [(2451545.0, 0.0), (2400000.5, 50123.2), (2450123.5, 0.2), (0.25, 2459143.0),
         (2378496.5, 0.3), (2488069.5, 0.75), (2305447.5, -0.1), (2560000.0, 0.49999999)]
print("const _DATES = [" + ", ".join(f"({r(a)}, {r(b)})" for a, b in dates) + "]")
for name in ("taitt", "tttai", "tttcg", "tcgtt", "tdbtcb", "tcbtdb"):
    vals = [tuple(float(v) for v in getattr(erfa, name)(a, b)) for a, b in dates]
    print(f"    :{name} => [" + ", ".join(f"({r(x)}, {r(y)})" for x, y in vals) + "],")
print("const _DTDB = [" + ", ".join(r(float(erfa.dtdb(a, b, 0.0, 0.0, 0.0, 0.0))) for a, b in dates) + "]")

cal = [(1583, 1, 1), (1858, 11, 17), (1900, 2, 28), (2000, 2, 29), (2016, 12, 31),
       (2400, 12, 31), (-4000, 3, 1)]
print("const _CAL2JD = [" + ", ".join(
    f"({y}, {m}, {d}) => ({r(float(erfa.cal2jd(y, m, d)[0]))}, {r(float(erfa.cal2jd(y, m, d)[1]))})"
    for y, m, d in cal) + "]")
jds = [(2451545.0, 0.0), (2451544.5, 0.0), (2451544.5, 0.5), (2451545.5, -0.5),
       (2400000.5, 50123.2), (2457754.0, 0.5 - 1e-10), (2459143.0, -0.25), (0.5, 2451544.0),
       (2299160.5, 0.0), (-68569.5, 0.0)]
print("const _JD2CAL = [" + ", ".join(
    f"({r(a)}, {r(b)}) => ({', '.join(r(v.item() if hasattr(v, 'item') else v) for v in erfa.jd2cal(a, b))})"
    for a, b in jds) + "]")

d2dtf = [("UTC", 2457753.5, 0.99999), ("UTC", 2457753.5, 0.999999999), ("UTC", 2457753.5, 0.9999942),
         ("TT", 2460462.5, 0.9999999954), ("UTC", 2460462.5, 0.9999999954), ("TAI", 2451544.5, 0.5),
         ("UTC", 2441499.5, 0.99999999)]
rows = []
for sc, a, b in d2dtf:
    iy, im, id_, ihmsf = erfa.d2dtf(sc, 3, a, b)
    h, m, s, f = (int(x) for x in ihmsf)
    rows.append(f'(:{sc.lower()}, {r(a)}, {r(b)}) => ({int(iy)}, {int(im)}, {int(id_)}, {h}, {m}, {s}, {f})')
print("const _D2DTF = [" + ", ".join(rows) + "]")

dtf2d = [("UTC", 2016, 12, 31, 23, 59, 60.5), ("UTC", 2016, 12, 31, 12, 0, 0.0),
         ("UTC", 2017, 1, 1, 0, 0, 0.0), ("TT", 2016, 12, 31, 12, 0, 0.0), ("TAI", 2024, 2, 29, 23, 59, 59.999)]
print("const _DTF2D = [" + ", ".join(
    f"(:{sc.lower()}, {y}, {m}, {d}, {h}, {mi}, {r(s)}) => ({r(float(erfa.dtf2d(sc, y, m, d, h, mi, s)[0]))}, "
    f"{r(float(erfa.dtf2d(sc, y, m, d, h, mi, s)[1]))})" for sc, y, m, d, h, mi, s in dtf2d) + "]")

utc = [(2457753.5, 0.5), (2457753.5, 0.99999), (2457754.5, 0.0), (0.25, 2457753.5), (2460462.5, 0.3)]
print("const _UTCTAI = [" + ", ".join(
    f"({r(a)}, {r(b)}) => ({r(float(erfa.utctai(a, b)[0]))}, {r(float(erfa.utctai(a, b)[1]))})" for a, b in utc) + "]")
tai = [(2457754.5, 0.000428), (2457754.5, 0.000417), (2457754.5, 0.000440), (2460462.5, 0.3)]
print("const _TAIUTC = [" + ", ".join(
    f"({r(a)}, {r(b)}) => ({r(float(erfa.taiutc(a, b)[0]))}, {r(float(erfa.taiutc(a, b)[1]))})" for a, b in tai) + "]")
