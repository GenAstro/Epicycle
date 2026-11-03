# Changelog

## 0.4.0

### Changes that can break a script

| In 0.3 | In 0.4 |
|:-------|:-------|
| `translate(earth, moon, jd_tdb)`, and the same for `translate_state` | `translate(earth, moon, t)`, where `t` is an AstroEpochs `Time` in any scale. A bare number is an error that shows the fix, `Time(jd, TDB(), JD())`. |
| `body_axes_rotation(model, naifid, jd_tdb)` and `iau2015_orientation(body, jd_tdb)` | Both take a `Time`. A custom orientation model's method takes one too, and reads `t.tdb.jd`. |
| `CelestialBody(name = "MyAsteroid", mu = 1.0, equatorial_radius = 1.0)` took Earth's NAIF ID and flattening | `name`, `mu`, `equatorial_radius` and `naifid` are required. `flattening` defaults to 0. |
| A body's `mu` or radius could be set to a negative, zero or nonfinite value on some routes | Every route checks the value. An integer value is stored as `Float64`. |
| `list_downloaded_spice_kernels()` and `list_cached_spice_kernels()` printed | `cached_spice_kernels()` returns the names of files on disk, and `loaded_spice_kernels()` the names of kernels in the session, of every kernel type. |
| `IauPolynomialOrientation` | `IAU2015`. The old name is removed. |
| `AstroUniverse.EARTH_DEFAULTS` | Removed. Read `earth`. |
| `println(body)` printed seven lines | `println(body)`, and a body inside a vector or a string, print one line. The block with units stays at the REPL. |

### Body constants

The built-in values now come from named sources, and results in some orbits move slightly.

- Earth keeps WGS84: μ 398600.4418 km³/s², radius 6378.137 km, flattening 1/298.257223563. The flattening was rounded to 0.00335281 before.
- Other bodies take μ from NAIF's `gm_de440.tpc`, using the system value for a planet with moons, since from Earth a planet and its moons act as one mass.
- Radii and flattening of other bodies come from NAIF's `pck00011.tpc`.

### New

- Named orientation models: `IAU2015`, `IAU1991` for Mars, `LunarPA` and `LunarME`.
- Earth orientation data has two sources for each series, the IERS data centre and the US Naval Observatory. A failed download leaves the cached copy in place.

### Fixed

- `using AstroUniverse` no longer reaches the network. Kernels are fetched at the first call that needs them, and if that fails, the error names each file, its URL and the folder to place it in.
- `eop_load` read a fixed-width file of the other theory's series into the wrong nutation fields, with no message. It now throws an error naming the mismatch.
- `eop_load` now checks a file's values as a download is checked, and an unreadable cached table is replaced.
- `eop(theory, ::Time)` was one second off in UT1 during a leap second.
- A correction carrying a derivative failed at a plain date.
- `translate` from several threads could abort the process. SPICE calls now share one lock.
- IAU2015 derivatives and `BigFloat` failed for seven of the nine bodies.
- Processes fetching one kernel at the same time failed. Each attempt now has its own temporary file.
- A truncated cached kernel loaded silently, and later calls failed with a SPICE handle error. A cached file of the wrong size is now fetched again.
- `unload_spice_kernel` could hang when a file of the same name had been loaded from another folder.
- `translate` with a NaN or a dual-valued time now gives a clear error.
- Texture paths were fixed at precompile time, and are now computed when the package loads.

### Removed

- The unused `Pkg` dependency.
