# EpicycleIO examples

These run. They were written against the interface before it was implemented, to be read and
argued with while changing it was still cheap, and the package caught up with them.

    julia --project=<environment> Ex_06_Porkchop.jl

The data is made up. `synthetic.jl` produces arrays shaped like real results from sines and
random numbers: no propagation, no `history`, no spacecraft. That was the claim being tested and
it held — the plotting side never knew where the numbers came from, so the arrays could have been
replaced by `history(Calc(...), ...)` calls with nothing else in these files changing.

| File | Shows | The question it is asking |
|---|---|---|
| `Ex_01_Basics.jl` | a time history, attributes, two spacecraft | is `xyplot(t, alt)` what you want to type? |
| `Ex_02_Components.jl` | a position column drawing as three lines | is `position[1]` the right automatic name? is the nesting rule too clever? |
| `Ex_03_Panels.jl` | named panels, replace versus add, layout | is a leading string the right way to name a panel? |
| `Ex_04_PolarPass.jl` | an az/el sky plot | should the compass convention be the default for polar? |
| `Ex_05_Residuals.jl` | residuals inside a 3σ band | is `band(x, lo, hi)` right, or centre-and-half-width? |
| `Ex_06_Porkchop.jl` | contours over a grid | is `Z` indexed `[arrival, departure]` the orientation you expect? |
| `Ex_07_GroundTrack.jl` | a ground track and stations | explicit `rad2deg`, or an equivalent of `thetaunit`? |
| `Ex_08_LiveUpdate.jl` | a solver publishing per iteration | one point per call, or a `push!` that does not resend? |
| `Ex_09_Report.jl` | `report`, and epochs rendered by the caller | is requiring the conversion right, or too blunt? |
| `Ex_10_EscapeHatch.jl` | a raw PlotlyBase trace | should the trace constructors be re-exported? |

## Reading order

`Ex_01` and `Ex_02` carry most of the interface. `Ex_04` and `Ex_05` are where the domain shows
up. `Ex_10` is the argument that we do not have to enumerate every plot type.

## What the examples already changed

Writing them surfaced one thing the spec had missed: the wrappers need bang forms.
`Ex_04` adds a second pass and a horizon mask to one sky plot, which needs `scatterpolar!`, and
the spec only ever defined `xyplot!`. Every wrapper now has one.
