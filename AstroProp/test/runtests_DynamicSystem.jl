
using Test
using AstroStates
using AstroBase
using AstroEpochs
using AstroUniverse
using AstroProp

# === Dummy state types ===
struct Cr <: AbstractParam end
struct TotalMass <: AbstractState end
struct TT_Time <: AbstractParam end
struct CustomUserVar1 <: AbstractParam end

# === Spacecraft instances ===
sat1 = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT())
)

sat2 = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT())
)

# === Force model setup ===
pm_grav = PointMassGravity(earth, ())
forces = ForceModel(pm_grav)

# === STM configuration ===
sat_stms = STMConfig([
    [PosVel(), TT_Time(), CustomUserVar1()],
    [PosVel(), TotalMass(), Cr()]
])

# === Dynamics system with STM configuration ===
dynsys = DynSys(
    forces = forces, 
    spacecraft = [sat1, sat2],
    stm_cfg = sat_stms
)

# === Verification ===
@assert dynsys.stm_cfg == sat_stms
@assert length(dynsys.spacecraft) == 2

println("✅ DynSys STM registration test passed.")
