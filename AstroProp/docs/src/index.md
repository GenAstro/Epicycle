```@meta
CurrentModule = AstroProp
```

# AstroProp

AstroProp provides force models, orbital propagators, and stopping conditions
for modelling spacecraft motion. AstroProp provides interfaces to the extensive 
numerical integration libraries in Julia's OrdinaryDiffEq.jl.  AstroProp is tested against the General Mission Analysis Tool (GMAT).

## Quick Start

The example below shows how to propagate a spacecraft to various stopping conditions:

```julia
using Epicycle

# Spacecraft
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
)sa

# Forces + integrator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
stop = StopAt(sat, PosDotVel(), 0.0; direction=+1)
propagate(prop, sat, stop)
println(get_state(sat, Keplerian()))

# Stop when |r| reaches 7000 km 
propagate(prop, sat, StopAt(sat, PosMag(), 7000.0))
println(get_state(sat, SphericalRADEC()))       

# Propagate backwards to node
stop = StopAt(sat, PosZ(), 0.0)
sol = propagate(prop, sat, stop; direction= :backward)

# Propagate multiple spacecraft with multiple stopping conditions
sc1 = Spacecraft(); sc2 = Spacecraft() 
stop_sc1_node = StopAt(sc1, PosZ(), 0.0)
stop_sc2_periapsis = StopAt(sc2, PosDotVel(), 0.0; direction=+1)
propagate(prop, [sc1,sc2], stop_sc1_node, stop_sc2_periapsis)

```

## Core Functions

```@docs
OrbitPropagator
IntegratorConfig
propagate
StopAt
```

## API Reference

```@autodocs
Modules = [AstroProp]
Order = [:type, :function, :macro, :constant]
Public = true
```
# Index

```@index
```