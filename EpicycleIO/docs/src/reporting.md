```@meta
CurrentModule = EpicycleIO
```

# Reporting

Use [`report`](@ref) to write named data columns to a text file. Reports are
useful for reviewing simulation results, sharing selected quantities, and
importing data into other analysis tools. Each keyword supplies a column name
and its values; the file contains a header followed by one row per sample.

## Write a report

This example writes three samples of elapsed time and altitude:

```julia
using EpicycleIO

elapsed_s = [0.0, 60.0, 120.0]
altitude_km = [400.0, 401.5, 403.0]

path = report("flight.txt";
              elapsed_s = elapsed_s,
              altitude_km = altitude_km)
```

The file contains the following values, with columns padded for readability:

```text
elapsed_s                 altitude_km
0.0                       400.0
60.0                      401.5
120.0                     403.0
```

Columns appear in the order supplied. The keyword names become the headings
exactly as written. Include units in those names where useful: `altitude_km`
is more informative than `altitude`. The writer does not convert units or add
metadata such as the coordinate system or time scale.

`report` returns the output path. Relative paths are resolved from Julia's
current working directory. Create the destination directory before writing.
An existing file is overwritten; each call writes a complete report rather
than appending rows.

## Include vector data

A column whose samples are vectors expands into one output column per
component. Supply a vector of sample vectors, with one sample for each row:

```julia
position_km = [[7000.0, 0.0, 0.0],
               [6985.0, 450.0, 0.0],
               [6940.0, 898.0, 0.0]]

report("position.txt";
       elapsed_s = elapsed_s,
       position_km = position_km)
```

The headings are `elapsed_s`, `position_km_1`, `position_km_2`, and
`position_km_3`. Component numbers follow the order within each sample. Use
vectors with the same number of components throughout the series.

To choose individual component names, supply them as separate columns:

```julia
report("position.txt";
       elapsed_s = elapsed_s,
       x_km = [r[1] for r in position_km],
       y_km = [r[2] for r in position_km],
       z_km = [r[3] for r in position_km])
```

The report does not record the position's reference frame. Identify it in the
headings or accompanying documentation when sharing the file.

## Report epochs

Convert Epicycle `Time` values to numbers or single-line strings before passing
them to `report`. For an existing sequence of epochs, `sample_epochs`, choose
the representation needed by the reader:

```julia
# sample_epochs contains Time values expressed in UTC.
report("epochs.txt";
       epoch_utc_isot = [t.isot for t in sample_epochs])

report("julian_dates.txt";
       epoch_utc_jd = [t.jd for t in sample_epochs])
```

Use `t.mjd` for modified Julian dates. These examples assume the epochs are
already in UTC; a column name does not change their time scale. Label the
heading with the actual scale and representation so readers can interpret the
values correctly.

Passing `Time` objects directly raises an `ArgumentError` because their printed
representation spans multiple lines.

## Choose a delimiter

By default, `report` separates columns with two spaces and pads each field to
a minimum width of 24 characters. Set `delim` to use another separator:

```julia
report("flight.tsv";
       delim = "\t",
       elapsed_s = elapsed_s,
       altitude_km = altitude_km)
```

Padding is retained when you change the delimiter. A comma delimiter is also
accepted; the writer does not quote or escape fields, so a value containing the
delimiter is refused, and with the default delimiter so is a value containing any
whitespace. Configure the receiving tool to trim padding if needed. The filename
extension does not select the delimiter.

## Formatting and data requirements

Integers are written as integers. Other real numbers are written as
floating-point values in the shortest form that reads back to the same value, so
a Julian date keeps its full precision. Other values are written using their
string representation.

Before writing a report:

- Supply at least one named column.
- Give every column the same number of samples.
- Align the samples yourself: values at the same index become one row. The
  writer does not match epochs, interpolate data, or resample series.
- Convert values whose printed representation spans multiple lines to numbers
  or single-line strings.

Missing columns, unequal column lengths, multiline values, and values containing
the delimiter raise an `ArgumentError`. Column names and lengths are checked before the file is opened,
but values are formatted while writing. A formatting error can therefore leave
a partially written file or replace an existing file with incomplete output.

For trajectory files that carry epochs, coordinate systems, and other metadata
in a supported exchange format, see [Ephemeris files](ephemeris.md).
