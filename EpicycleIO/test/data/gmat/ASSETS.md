# Data files

Sample ephemerides that ship with GMAT R2026a (`data/vehicle/ephem/` in the release), read by
`test_ephemeris_oem.jl` and `test_ephemeris_stk.jl` as files another tool wrote. Each is trimmed
to its header, its metadata, its first 20 data lines and its last 2, so the span, the declared
start, stop and usable times, and the final step are those of the original. Line endings are LF.

| File | Original in the release | Data lines kept | What it covers |
|---|---|---|---|
| `DemoLEO-1.oem` | `ccsds/DemoLEO-1.oem` | 22 of 577 | OEM 1.0 as GMAT writes it; a low Earth orbit at 5 min steps |
| `DemoLuna-1.oem` | `ccsds/DemoLuna-1.oem` | 22 of 2161 | OEM 1.0; an Earth-centred orbit at lunar distance |
| `SampleOEMEphem.oem` | `ccsds/SampleOEMEphem.oem` | 22 of 2881 | OEM 1.0 from another originator: a header comment, mixed-case `Lagrange`, microsecond epochs |
| `SampleSTKEphem.e` | `stk/SampleSTKEphem.e` | 22 of 2610 | `stk.v.10.0`: three-digit exponents, no interpolation keywords. `NumberOfEphemerisPoints` is changed from 2610 to 22 to match |

GMAT is released by NASA under the Apache License 2.0.
