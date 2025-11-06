# Component Cheat Sheets

## A Julia Language Cheat Sheet

A Julia Language [Cheat Sheet](https://cheatsheet.juliadocs.org/).

## Time Systems (AstroEpochs)

AstroEpochs provides comprehensive time system handling for astrodynamics applications. It supports multiple time scales (UTC, TAI, TT, TDB) and formats (Julian Date, Modified Julian Date, ISO 8601) with high-precision conversions between them.

```julia
using AstroEpochs

# Create from Julian Date
t1 = Time(2451545.0, TT(), JD())

# Create from Modified Julian Date  
t2 = Time(51544.5, UTC(), MJD())

# Create from ISO string
t3 = Time("2000-01-01T12:00:00.000", TAI(), ISOT())

# Access different representations
t1.jd        # Julian Date value
t1.mjd       # Modified Julian Date value  
t1.isot      # ISO 8601 string

# Convert between scales (creates new Time object)
t_utc = t1.utc
t_tdb = t1.tdb
```

See the full [Reference Material](https://genastro.github.io/Epicycle/AstroEpochs/dev/) for more details.

## Orbital States (AstroStates)

AstroStates handles spacecraft orbital state representations and conversions between state representations. It supports 10 representations include Cartesian, Keplerian, B-plane and other representations.

```julia
using AstroStates

# Define a Cartesian state
cart = CartesianState([7000.0, 0.0, 100.0, 0.0, 7.5, 2.5])

# Convert to Keplerian then back to Cartesian
mu = 398600.4418 
kep   = KeplerianState(cart, mu)     
cart2 = CartesianState(kep, mu)     

# Display some state elements
kep.sma
kep.raan

# Generate a vector containing the state struct data
to_vector(kep)

# See a list of all supported representations
subtypes(AbstractOrbitState)
```

See the full [Reference Material](https://genastro.github.io/Epicycle/AstroStates/dev/) for more details.

## Spacecraft Modeling (AstroModels)

AstroModels provides the spacecraft model defining properties such as time, state, and mass.

```julia
using AstroModels, AstroStates, AstroEpochs

# Method 1: Using a CartesianState struct
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = 1000.0
)

# Method 2: Direct construction with OrbitState
sc2 = Spacecraft(
    state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03], Cartesian()),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = 1000.0
)
``` 

See the full [Reference Material](https://genastro.github.io/Epicycle/AstroModels/dev/) for more details. 

## Celestial Bodies (AstroUniverse)

AstroUniverse provides access to celestial body data including gravitational parameters, physical properties, and NAIF identification codes. It includes predefined bodies and supports creation of custom celestial objects for specialized applications and performs translations between coordinate origins. 

```julia
using AstroUniverse

# Access predefined celestial bodies
earth.mu
venus.naifid

# Create a custom body
phobos = CelestialBody(
    name = "Phobos",
    naifid = 401,                    # NAIF ID for Phobos
    mu = 7.0875e-4,       # km³/s² (gravitational parameter)
    equatorial_radius = 11.1,                   # km (mean radius)
)

```

See the full [Reference Material](https://genastro.github.io/Epicycle/AstroUniverse/dev/) for more details.
## Maneuvers (AstroMan)

AstroMan provides maneuver modeling capabilities for trajectory modifications. It supports impulsive maneuvers with various coordinate frame options and specific impulse specifications for realistic propulsion modeling.

```julia
using Epicycle
m = ImpulsiveManeuver(axes=Inertial(), 
                      Isp=300.0, 
                      element1=0.01, 
                      element2=0.0, 
                      element3=0.0)

sc = Spacecraft()
maneuver(sc, m)
```

See the full [Reference Material](https://genastro.github.io/Epicycle/AstroMan/dev/) for more details.

## Calculations Framework (AstroFun)

AstroFun provides a unified calculation framework for extracting and setting orbital parameters, celestial body properties, and maneuver characteristics. It offers a consistent interface for accessing computed quantities across the Epicycle ecosystem.

```julia
using AstroFun, AstroStates, AstroModels

# Create a spacecraft with orbital state
sc = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]), 
                time = Time("2024-01-01T12:00:00", UTC(), ISOT()), 
                mass = 1000.0)

# Get semi-major axis from current state
sma_calc = OrbitCalc(sc, SMA())
a = get_calc(sma_calc)           
set_calc!(sma_calc, 10000.0)  

# Set target incoming asymptote (rp = 6900, C3 = 14.0)
hyp = OrbitCalc(sc, IncomingAsymptote())
set_calc!(hyp, [6900.0, 14.0, 0.0, 0.0, 0.0, 0.0])  
    
# Set and get Earth's mu
mu_calc = BodyCalc(earth, GravParam())
μ = get_calc(mu_calc)            
set_calc!(mu_calc, 3.986e5)      

# Set and get maneuver elements
toi = ImpulsiveManeuver()
dvvec_calc = ManeuverCalc(toi, sc, DeltaVVector())
Δv = get_calc(dvvec_calc)   
set_calc!(dvvec_calc, [0.2, 0.3, 0.4])
```

See the full [Reference Material](https://genastro.github.io/Epicycle/AstroFun/dev/) for more details.


