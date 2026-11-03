# Changelog

## 0.7.0

### Changes that can break a script

| In 0.6 | In 0.7 |
|:-------|:-------|
| `Time(value, :tt, :jd)` | `Time(value, TT(), JD())`. Symbols are no longer accepted, and the error shows the tag form. |
| `t.scale == :tt` | `t.scale == TT()`. `t.scale` and `t.format` return tags. |
| `==` was false for the same date in different formats | `==` compares the scale and stored date parts. Display format is ignored, and equal times have equal hashes. |
| `println(t)` printed four lines | `println(t)`, and a time inside a vector or a string, print one line. The block stays at the REPL. |
| `t.jd` printed in the block as `2.451545e6` | The block prints a Julian date in fixed notation, `2451545.000000000`. |
| `tai_minus_utc(days)` and `tdb_minus_tt(jd1, jd2)`, both unexported | Now private. Use `tai_minus_utc(t)` and `tdb_minus_tt(t)`. |
| Rebuilding a time from `t.mjd` or `t.isot` to change its display | `Time(t, MJD())` and `Time(t, ISOT())` preserve the stored date parts. |

### New

- `Time(t, scale)`, `Time(t, format)`, `Time{T}(t)`, and `convert`, to change the scale, display format, or storage type of an existing time.
- `Time(year, month, day, hour, minute, second, scale, format)`, from calendar fields.
- `isot(t; digits)`, a timestamp with 0–12 decimal places in the seconds field. Extra displayed digits do not guarantee extra accuracy.
- `tai_minus_utc(t)` and `tdb_minus_tt(t[, position])`, the offsets between scales.
- `<`, `isless`, `sort`, `minimum`, `maximum`, and `isapprox` with a tolerance in seconds.
- `t .+ offsets`, since a time is one value in a broadcast.
- `Time` from any string type.

### Fixed

- A whole-number Modified Julian Date threw `InexactError`.
- A truncated or altered leap-second list could be accepted and stored. Validation now checks its hash and expected historical entries through 2017. A truncated file fails validation; comments and lines after the hash are not covered.
- An unusable depot stopped the first UTC conversion. The list now goes to a temporary folder for the session, with a warning.
- Processes starting together broke each other's download of the leap-second list.
- Precompiling a package that converts a UTC time could reach for the network.
- Subtracting an unsigned interval wrapped: `t - 0x01` moved the time 255 days forward.
- The derivative of `tdb_minus_tt(t, position)` was NaN at exactly UTC midnight.
- `tai_minus_utc(t)` returned the next value in the last tens of microseconds of a leap second.
- A nonfinite date threw during UTC conversion and display. It now passes through those operations, but calendar timestamp formatting still throws.
