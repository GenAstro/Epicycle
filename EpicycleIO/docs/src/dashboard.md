```@meta
CurrentModule = EpicycleIO
```

# The dashboard

The EpicycleIO dashboard brings interactive plots and animated spacecraft
trajectories into one browser page. Use it to compare results, inspect a
propagation, or watch a solver converge while your Julia session runs.

Your first plotting call starts a local server and opens the dashboard in your
default browser. Subsequent calls update the same set of panels. No account or
API key is required.

## Create your first dashboard

Create two panels to compare altitude and speed:

```julia
using EpicycleIO

hours = collect(0.0:0.1:6.0)
altitude_km = 400.0 .+ 20.0 .* sin.(hours)
speed_km_s = 7.67 .- 0.01 .* sin.(hours)

xyplot("Altitude", hours, altitude_km; mode = "lines", name = "Nominal")
panel!("Altitude";
       xaxis_title = "Hours from epoch", yaxis_title = "Altitude (km)")

xyplot("Speed", hours, speed_km_s; mode = "lines", name = "Nominal")
panel!("Speed";
       xaxis_title = "Hours from epoch", yaxis_title = "Speed (km/s)")
```

Both panels appear in the same dashboard. Hover over a plot to inspect values
and use the plot's toolbar to explore the data. Animated trajectory panels
created with `orbitview` can appear alongside data plots; see
[3D plotting](plotting_3d.md) for trajectory setup and playback controls.

## Arrange panels in the browser

Use the **grid / tabs** button at the top of the dashboard to switch layouts:

- **Grid** shows panels together for comparison. The number of columns adjusts
  to the number of panels.
- **Tabs** shows one panel at a time. Select a panel's name to display it.

To give a group of panels its own browser tab, pass their names to
`open_dashboard`:

```julia
open_dashboard("Altitude", "Speed")  # Show these two panels.
open_dashboard("Altitude")           # Give one panel the full view.
open_dashboard()                     # Show all panels.
```

Each browser tab is a view of the same panels in Julia. Updating a panel updates
every tab that displays it. A tab opened with selected panel names continues to
show only those panels as you create others.

## Add data and reuse panels

The leading string in a plotting call identifies the destination panel. Omit it
to use the default panel, `"Plot"`. The `name` keyword labels a data series, also
called a trace, in the panel's legend.

The replacing and adding forms let you control what happens to existing data:

| Call | Effect |
|:--|:--|
| `xyplot("Altitude", x, y)` | Replace the panel's data and reset its interactive view. |
| `xyplot!("Altitude", x, y)` | Add a trace while keeping existing traces and the current zoom. |
| `panel!("Altitude"; ...)` | Update layout settings while keeping previously supplied settings. |

For example, add a second altitude history to the first example:

```julia
xyplot!("Altitude", hours, altitude_km .+ 30;
        mode = "lines", name = "Alternative", line_dash = "dash")
panel!("Altitude"; showlegend = true)
```

The panel now contains two traces, `Nominal` and `Alternative`. Reusing panel
names lets you rerun a script without creating duplicate panels. When rerunning
a sequence of adding calls, first replace or clear the panel to avoid accumulating
copies of the traces.

The other data plotting functions use the same replacing and adding convention.
See [Data plotting](data_plotting.md) for plot types and accepted data shapes.

## Configure titles and axes

Use `panel!` for titles, axis labels, limits, and legends. Keyword names follow
Plotly layout properties, with underscores for nested settings:

```julia
panel!("Altitude";
       title = "Altitude comparison",
       xaxis_title = "Hours from epoch",
       yaxis_title = "Altitude (km)",
       xaxis_range = [0, 6],
       yaxis_range = [350, 500],
       legend_orientation = "h")
```

Layout settings survive data replacement. A later `panel!` call changes only the
properties you supply. For example, `xaxis_range = [6, 0]` reverses the horizontal
axis, and `yaxis_autorange = true` restores automatic vertical limits.

You can configure a panel before adding data. If other panels already exist,
`panel!` warns when the supplied name is new. Check the spelling if you intended
to update an existing panel.

## Update results during a run

Call plotting functions from a propagation or solver callback to display results
as they become available. Publishing is asynchronous, and rapid updates are
combined; the browser may skip intermediate states.

Updates reach the browser only when Julia has a chance to run other tasks, which
happens at any `sleep`, I/O, or `yield()`. A tight loop that does none of these
holds back every update until it finishes; call `yield()` after plotting in such
a loop.

For a growing history, keep the accumulated values and replace the plotted
series at a suitable interval:

```julia
iterations = Int[]
costs = Float64[]

for k in 1:100
    push!(iterations, k)
    push!(costs, 10.0 / k)  # Replace with the current solver cost.
    if k % 10 == 0
        xyplot("Convergence", iterations, costs; mode = "lines", name = "Cost")
    end
end
panel!("Convergence"; xaxis_title = "Iteration", yaxis_title = "Cost")
```

Each `xyplot!` call creates a new trace; it does not append samples to an
existing trace. Use it to compare separate series rather than to add every
sample of a long run. Replacing the data resets the view to fit the new results.

## Add captions to a demonstration

Use `EpicycleIO.narrate` to explain results above the panels during a tutorial,
presentation, or live run. A caption can include a heading, explanatory text,
labeled values, a status message, and a code sample:

```julia
EpicycleIO.narrate("Compare the nominal and alternative altitude histories.";
    title = "Altitude comparison",
    subtitle = "Six-hour example",
    facts = ["Duration" => "6 h", "Altitude offset" => "30 km"],
    status = "Ready to inspect",
    code = "open_dashboard(\"Altitude\", \"Speed\")",
    panels = ["Altitude", "Speed"])

open_dashboard("Altitude", "Speed")
```

Each call replaces the current caption. Text is plain text, with line breaks
preserved. Format numbers and units before passing them in `facts`. The `code`
field displays source text; it does not execute it.

The full dashboard shows the latest caption. Supplying `panels` also keeps a
caption for browser tabs showing exactly that set of panels, in any order.
Those tabs retain their caption when later captions describe other panel sets.

Remove all captions with:

```julia
EpicycleIO.clear_narration()
```

## Clear results and manage the session

| Function | Effect |
|:--|:--|
| `clear!("Altitude")` | Empty the panel, retaining its layout and place on the dashboard. |
| `clear!()` | Empty the default panel, `"Plot"`. |
| `clear_all!()` | Remove every panel and its layout. |
| `EpicycleIO.close_dashboard()` | Stop the local server, retaining the panels in Julia. |

Clearing a panel that does not exist has no effect. After stopping the server,
the next plotting call starts it again.

Closing a browser tab leaves the panels in Julia. If no dashboard viewer is
detected, a later plotting call can open another tab. Panels belong to the
current Julia session; ending that session ends live updates.

### Control automatic browser opening

For an automated run, disable automatic browser opening before plotting:

```julia
EpicycleIO.auto_open!(false)
```

Plotting and publishing remain enabled. Open the dashboard explicitly with
`open_dashboard()`, or restore automatic opening with
`EpicycleIO.auto_open!(true)`.

To obtain the address without opening a browser, use:

```julia
EpicycleIO.dashboard_url()
EpicycleIO.dashboard_url("Altitude", "Speed")
```

These calls start the local server if needed. The address uses `127.0.0.1` and
is intended for a browser on the same machine as Julia.

## Browser requirements and troubleshooting

The dashboard loads Plotly.js from a public content delivery network. It loads
CesiumJS when a trajectory view is first displayed, and trajectory views request
Earth imagery from an external map service. An internet connection is needed
when these resources are not cached by the browser.

| Symptom | What to check |
|:--|:--|
| The browser does not open automatically | Check whether automatic opening is disabled. Call `open_dashboard()` or obtain the address with `EpicycleIO.dashboard_url()` and open it manually. |
| A panel remains blank | Check that the browser can load the required external libraries and, for trajectory views, map imagery. |
| The dashboard reports that it is disconnected | Check that Julia is still running and the local server has not been stopped. A new plotting call restarts a stopped server. |
| A panel is missing from a selected view | Check the panel name and the names passed to `open_dashboard`. Use `open_dashboard()` to view all panels. |
| Earlier data disappear after a plotting call | Use the adding form, such as `xyplot!`, to retain existing traces. |
| Traces multiply when a script is rerun | Start the panel with a replacing call or use `clear!` before adding traces. |
| Axis labels or legend settings do not change | Apply layout settings with `panel!` and check the panel name. |

If EpicycleIO cannot launch the browser, it issues a warning containing the
local dashboard URL. Open that address manually on the machine running Julia.

See [Data plotting](data_plotting.md), [3D plotting](plotting_3d.md), and the
[API reference](api.md) for more detail.
