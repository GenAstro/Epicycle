```@meta
CurrentModule = EpicycleIO
```

# API

Every name `EpicycleIO` exports. For what an attribute does, see
[Plotly's reference](https://plotly.com/javascript/reference/) — attribute names here are
Plotly's own.

To see this list from the REPL rather than reading it here:

```julia
names(EpicycleIO)
```

```@autodocs
Modules = [EpicycleIO]
Order   = [:function]
Filter  = t -> nameof(t) in names(EpicycleIO)
```

```@autodocs
Modules = [EpicycleIO]
Order   = [:type]
Filter  = t -> nameof(t) in names(EpicycleIO)
```

## Qualified names

These are not exported; call them with the package name, as in
`EpicycleIO.close_dashboard()`.

```@docs
EpicycleIO.close_dashboard
EpicycleIO.dashboard_url
EpicycleIO.auto_open!
EpicycleIO.narrate
EpicycleIO.clear_narration
EpicycleIO.Diagnostic
```

```@index
Pages = ["api.md"]
```
