# Changelog

## 0.7.0

### Changes that can break a script

| In 0.6 | In 0.7 |
|:-------|:-------|
| `Time(value, :tt, :jd)` | `Time(value, TT(), JD())`. Symbols are no longer accepted, and the error shows the tag form. |
| `t.scale == :tt` | `t.scale == TT()`. `t.scale` and `t.format` return tags. |
| `==` was false for the same date in different formats | `==` compares the scale and the date. The format plays no part, and `hash` agrees. |
| `println(t)` printed four lines | `println(t)`, and a time inside a vector or a string, print one line. The block stays at the REPL. |
| `t.jd` printed in the block as `2.451545e6` | The block prints a Julian date in fixed notation, `2451545.000000000`. |
| `tai_minus_utc(days)` and `tdb_minus_tt(jd1, jd2)`, both unexported | Now private. Use `tai_minus_utc(t)` and `tdb_minus_tt(t)`. |
| Rebuilding a time from `t.mjd` or `t.isot` to change its display | `Time(t, MJD())` and `Time(t, ISOT())` hold the date exactly. |

### New

- `Time(t, scale)`, `Time(t, format)`, `Time{T}(t)` and `convert`, to change a time you hold.
- `Time(year, month, day, hour, minute, second, scale, format)`, from calendar fields.
- `isot(t; digits)`, a timestamp with 0 to 12 decimal places of seconds.
- `tai_minus_utc(t)` and `tdb_minus_tt(t[, position])`, the offsets between scales.
- `<`, `isless`, `sort`, `minimum`, `maximum`, and `isapprox` with a tolerance in seconds.
- `t .+ offsets`, since a time is one value in a broadcast.
- `Time` from any string type.

### Fixed

- A whole-number Modified Julian Date threw `InexactError`.
- A truncated or altered leap-second list was accepted and stored, and gave a wrong offset until it expired. A list must now match its hash and hold every change through 2017.
- A depot that could not be written stopped the first UTC conversion.
- Processes starting together broke each other's download of the leap-second list.
- Precompiling a package that converts a UTC time could reach for the network.
- Subtracting an unsigned interval wrapped: `t - 0x01` moved the time 255 days forward.
- The derivative of `tdb_minus_tt(t, position)` was NaN at exactly UTC midnight.
- `tai_minus_utc(t)` returned the next value in the last tens of microseconds of a leap second.
- A time with no finite date threw from UTC conversion and from printing. It now passes through, prints, and has no calendar date to write.
