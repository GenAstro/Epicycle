```@meta
CurrentModule = EpicycleIO
```

# Reporting

`report` writes named columns to a delimited text file: the same arrays that would be plotted,
in a form that can be read, diffed, or loaded by something else.

Like the plotting functions, `report` operates on arrays and has no knowledge of where they came
from, so the examples below build their own data and run as written.

```julia
using EpicycleIO

t = collect(0.0:0.25:6.0)
r = [[7000.0 + 10i, 100.0i, 1300.0] for i in t]
v = [[0.0, 7.35 - 0.001i, 1.0] for i in t]

report("flight.txt"; time = t, position = r, velocity = v)
```

## Column Headings

The columns are named by the keywords given, rather than by whatever the quantity is called
internally. A report is usually read by someone who did not perform the run, and `time` tells
that reader more than `epoch` does.

A vector-valued column expands into one column per component:

```
time        position_1    position_2    position_3    velocity_1  …
```

Columns must be the same length. Columns from two different recordings have their own sample
times and cannot share a report, so a length mismatch raises an error naming the columns and
their lengths rather than writing a file that silently misaligns.

Columns carry no units of their own. They are written as given, so a keyword should name what
the number is when the reader will need to know.

```julia
using EpicycleIO

t   = collect(0.0:0.25:6.0)
alt = 400.0 .+ 25.0 .* sin.(t)

report("flight.txt"; hours_from_epoch = t, altitude_km = alt)
```

## Epochs Are Rendered by the Caller

`report` does not interpret epochs, and this is deliberate. A bare Julian date does not record
whether it is TT or UTC, and a text writer is the wrong place to make that decision. A `Time`
passed directly raises an error naming the conversions available.

```julia
using Epicycle
using EpicycleIO

t = [Time(2458849.5 + 0.01k, 0.0, :tdb, :jd) for k in 0:5]
alt = 400.0 .+ 25.0 .* collect(0.0:5.0)

# Choose the scale and format, and the file then records which was meant
report("flight.txt"; time = [x.jd   for x in t], altitude = alt)
report("flight.txt"; time = [x.isot for x in t], altitude = alt)
```

The same rule catches anything else that would break a row. A value whose printed form spans
several lines would turn one row into several without any error, so `report` refuses it and says
what to convert.

## The Same Data, Reported and Plotted

Both take the arrays from one `history` call. There the resemblance stops, and it should. A file
has column headings that must be called something, so `report` names them with keywords. A plot
has series, which Plotly names with `name`. The two are alike in taking the same data, not in
their syntax.

```@raw html
<!-- doc-fragment -->
```
```julia
t, r = history(Calc(epoch, sat), Calc(position_vector, sat, EarthMJ2000Eq))

report("flight.txt"; time = [x.jd for x in t], position = r)
xyplot("Position", t, r; name = ["x", "y", "z"])
```
