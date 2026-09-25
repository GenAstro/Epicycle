# Manual tests

The automated suite never opens a browser. It checks that the JSON we send is the JSON we
meant to send, and it cannot check whether Plotly agrees, whether the picture is legible, or
whether the figure means what its caption claims.

Every bug found in EpicycleIO so far was found by looking at a screen, and every one of them
passed the whole suite. That is not a gap to close — it is the boundary. These procedures cover
the other side of it.

**What only a person can judge**

- whether Plotly drew anything at all
- whether what it drew is the plot that was asked for
- whether the data means what the title says
- whether it is readable at the size it actually appears
- whether interaction — zoom, live update, closing a window — behaves

---

## Setup

```julia
julia --project=<environment>
using EpicycleIO
```

A browser tab opens on the first plot. The examples live in the package, at
`joinpath(pkgdir(EpicycleIO), "examples")`.

Between checks, `clear_all!()` clears the board — panels vanish from the browser too.

---

## Two-minute smoke test

Run after any change. If all four pass, nothing structural is broken.

```julia
include(".../examples/Ex_01_Basics.jl")
```

1. **Something is drawn.** A line, cyan, in a dark panel.
2. **The axes are labelled.** Not bare numbers — `Ex_03` sets titles and they must appear.
3. **The legend names the series.** "Sat A" and "Sat B", not "trace 0".
4. **Re-running leaves one panel, not two.** Run the same include twice.

---

## Regression checks

Each of these is a bug that shipped. They are listed by what wrong looks like, because in every
case wrong looked plausible.

### R1 — `panel!` actually applies

```julia
include(".../examples/Ex_03_Panels.jl")
```

Look at the **Decay** panel.

- **Right:** y axis is logarithmic — gridlines bunched toward the top, not evenly spaced. x axis
  runs 6 on the left down to 0 on the right. Both axes have titles.
- **Wrong:** evenly spaced y gridlines, x running 0 to 6, no axis titles.

*What happened:* layout keywords were published flat, as `{"yaxis_type":"log"}`. plotly.js does
not know that spelling and ignores it silently, so every `panel!` call in every example did
nothing while the plots still looked fine.

### R2 — polar radial units

```julia
include(".../examples/Ex_04_PolarPass.jl")
```

Look at **Sky view**.

- **Right:** each pass rises from the rim, arcs toward the centre, and returns to the rim. The
  10 degree mask is a ring near the outer edge. 0° is at the top and 90° is to the **right**.
- **Wrong:** the passes hug the outer edge as thin arcs and never approach the centre. Or 90° is
  on the **left**, which means the azimuth is mirrored.

*What happened:* `thetaunit` converts the angle; Plotly has no equivalent for the radial axis.
Elevation in radians against a 0–90 degree range put every point within 1.6% of the rim.

The mirrored case is worth staring at. A mirrored compass plot is not obviously broken — it is a
plausible-looking picture of a pass that went the other way.

### R3 — the figure means what the caption says

```julia
include(".../examples/Ex_10_EscapeHatch.jl")
```

- **Monte Carlo dispersion:** one box, centred near 501 km, whiskers spanning roughly 480–525.
  It is a spread of outcomes over 500 runs.
- **Residual distribution:** two violins. "range" is a narrow shape centred on zero; "range-rate"
  is visibly **fatter and offset upward**. That difference is the point — it is what a wrong
  noise model looks like.
- **Wrong:** two violins the same shape, or a box spanning 440–560 with a flat top and bottom
  (that is the distribution of a sine wave, not a dispersion).

### R4 — box and violin draw at all

Same panels as R3.

- **Right:** a filled box, and two filled violin shapes.
- **Wrong:** correct axes and legend, but empty plot area.

*What happened:* an explicit category column was added to fix a tick label and Plotly stopped
drawing entirely. The payload was valid; only the browser knew.

---

## Rendering

### M1 — every example draws

Run all ten. Every panel must contain marks.

```julia
for f in readdir(".../examples", join=true)
    endswith(f, ".jl") && startswith(basename(f), "Ex_") && include(f)
end
```

Check each panel has visible data. An empty panel with correct axes is the failure this whole
document exists for.

### M2 — the picture is the right kind

| Panel | Should be |
|---|---|
| Altitude, Position | lines |
| Sky view | polar, passes arcing inward |
| Range residuals | grey band with points scattered inside it |
| Porkchop | contours with labelled levels and a colour bar |
| Ground track | coastlines with a cyan track and three orange stations |
| Convergence | markers descending toward a floor |

### M4 — the 3D view

```julia
include(".../examples/Demo_PropagateAndView.jl")
```

Coast, burn, coast, so there are two propagation segments and one burn of 0.20 km/s.

- **Right:** a globe with the trajectory drawn around it in **two colours**, one per segment,
  and **one point where the colours meet, labelled `maneuver ΔV=0.2 km/s`**. Day and night
  sides of the body are distinguishable. Play runs the spacecraft along the path; the clock
  time under the view advances.
- **Wrong:** one colour where there should be two, which means segment extraction lost the
  burn. Or the globe untextured, which means the imagery layers failed to load.

The ΔV on that label is worth reading rather than glancing at. It is not the number that was
commanded — it is differenced back out of the recorded velocities either side of the burn, so
if it reads 0.2 the extraction found the right pair of arcs. A wrong number here means the
trajectory is being drawn from the wrong samples, which the picture alone would not show.

!!! note "One point per burn"
    Three burns must give three points. They previously gave one: the packet id was built from
    the maneuver's name, `solve_trajectory!` names every burn "maneuver", and Cesium merges
    packets that share an id. If a multi-burn run ever shows fewer points than burns, that is
    where to look.

**The check worth staring at: press play and watch the orbit against the stars.** The camera is
locked in the inertial frame, so the orbit should hold still and the body should turn beneath
it. If the orbit drifts or appears to precess, the ICRF camera lock did not survive — that is
the setting most likely to have been disturbed in the port, and the one that looks plausible
while being wrong.

Then check the panel bar: play toggles to pause, reset returns to the start, slower and faster
change the rate and the readout agrees, stars toggles the sky.

### M5 — a view and plots on one page

Same script. The semi-major axis plot should show a step at the same moment the trajectory
changes colour. That pairing is the reason the globe and the plots share a page, and it is the
one thing neither could show alone.

Cesium should not load at all on a dashboard with no 3D panel — run `Ex_01_Basics.jl` in a fresh
session and check the network tab stays free of `cesium.js`.

### M3 — legibility

With all panels in the grid at once:

- axis numbers readable without leaning in
- series distinguishable by colour
- titles not overlapping the plot area — the polar 0° label is the one that collides
- legends not covering data

This is a judgement call, and it is the one the automated tests can never make. Principle D1
says user experience is the product; a correct plot that cannot be read fails anyway.

---

## Interaction

### I1 — live update, and zoom survives it

```julia
clear_all!()
for k in 1:200
    xyplot!("Live", [Float64(k)], [sin(k/10)]; mode = "markers")
    sleep(0.05)
end
```

- points appear progressively, not all at the end
- **while it runs**, zoom into the panel with the scroll wheel. The view must stay where you put
  it as new points arrive. If it snaps back on every update, `Plotly.react` has been replaced by
  a full redraw somewhere.
- the REPL stays responsive throughout — publishing must not block the loop

### I2 — a closed window comes back

1. Close the browser tab.
2. `xyplot(1:10, rand(10))`
3. A new tab must open.

*What happened:* the flag suppressing repeat launches was a one-way latch, so once a tab had
existed the package believed one was there forever.

### I3 — grid, tabs, and scoped windows

- the `grid`/`tabs` button top right switches between all panels at once and one at a time
- `open_dashboard("Sky view")` opens a tab with only that panel
- a panel shown in two tabs updates in both
- resizing the window reflows the grid and the plots resize with it

### I4 — the server stops cleanly

```julia
EpicycleIO.close_dashboard()
```

The tab should report "disconnected — is Julia still running?" rather than hanging or erroring.
Then `xyplot(1:10, rand(10))` must bring everything back.

---

## Reporting

### P1 — a report is readable by someone who did not write it

```julia
include(".../examples/Ex_09_Report.jl")
```

Open `flight.txt`. Headings are your keywords, `position` has become `position_1`, `position_2`,
`position_3`, and the columns line up.

### P2 — an epoch is refused with an explanation

Passing `Time` values to `report` must fail and the message must name `t.jd`, `t.mjd` and
`t.isot`. This one is checked automatically too, but the wording is worth reading as a user
would, because that is the point of it.

---

## Recording a run

Date, Julia version, PlotlyBase version, browser. Then the check ids that passed, and for
anything that failed, a screenshot — a screenshot has found four bugs here and a description
has found none.

```
2026-09-01  julia 1.12.4  PlotlyBase 0.8.23  Edge
smoke 1-4 ok · R1-R4 ok · M1-M3 ok · I1-I4 ok · P1-P2 ok
```
