# Graph walk/cache probe
# Run after loading prototype:
# include("FrameTransformPrototype.jl")

using StaticArrays

println("\n=== Graph Walk Cache Probe ===")

# Start from empty caches
clear_transform_caches!()
println("Initial PRECOMPUTED_PATH_CACHE size: ", length(PRECOMPUTED_PATH_CACHE))
println("Initial SIM_PLAN_CACHE size: ", length(SIM_PLAN_CACHE))

# Build a transform pair that is typically not precomputed in normal smoke tests
# PEF -> CIRS should route through: PEF -> ITRF -> TIRS -> CIRS

time = Time("2004-04-06T12:00:00.000", UTC(), ISOT())
pef_frame = CoordinateFrame(earth, PEF())
cirs_frame = CoordinateFrame(earth, CIRS())

coord_pef = Coordinate(
    SVector{3,Float64}(7000.0, 100.0, 50.0),
    SVector{3,Float64}(0.0, 7.5, 1.0),
    time,
    pef_frame,
)

before_l2 = length(PRECOMPUTED_PATH_CACHE)
before_l1 = length(SIM_PLAN_CACHE)

coord_cirs = transform(coord_pef, cirs_frame)

after_l2 = length(PRECOMPUTED_PATH_CACHE)
after_l1 = length(SIM_PLAN_CACHE)

meta = edge_transform_metadata_from_context(coord_pef, cirs_frame, nothing)
entry = get_precomputed_path(PEF, CIRS, meta)

println("\nPRECOMPUTED_PATH_CACHE: ", before_l2, " -> ", after_l2)
println("SIM_PLAN_CACHE: ", before_l1, " -> ", after_l1)
println("Resolved path PEF -> CIRS: ", join(string.(entry.path), " -> "))

println("\nResult in CIRS:")
println(coord_cirs)

println("\n=== Probe Complete ===")
