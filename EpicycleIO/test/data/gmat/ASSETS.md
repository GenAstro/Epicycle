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

## Moon-centred reference files

`gmat_luna.oem`, `gmat_luna.e` and `gmat_luna.bsp` are one low lunar orbit, two hours of
point-mass Moon gravity at 60 s steps, written by GMAT R2026a in all three formats from
`EpicycleIO/test/truth/ephemeris_truth_luna.script`. The Moon's gravitational parameter was set
to Epicycle's, 4902.8 km³/s², so `test_ephemeris_api.jl` checks each file against the exact
two-body orbit. Sizes and SHA-256 as stored in the repository (text files with LF line endings;
GMAT on Windows writes CRLF, which git normalizes, and the readers accept either):

| File | Bytes | SHA-256 |
|---|---|---|
| `gmat_luna.oem` | 20824 | `8ff9ff5f15123d69518be4bf8be969c1efaa15669aa430eff6478e9d187327be` |
| `gmat_luna.e` | 20349 | `24fc0a9ccb986648c94a5aed7e9aabbde6d40e3a66be3f5465a8f24c8f30309d` |
| `gmat_luna.bsp` | 14336 | `5b98cb52333151e94738e9de48c1687a9d54264195bf2e919ba9932d1692a110` |

GMAT is released by NASA under the Apache License 2.0.
