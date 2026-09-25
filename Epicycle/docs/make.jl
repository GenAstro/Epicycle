# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using Epicycle
using Documenter

# The umbrella's pages document names owned by the library packages, so Documenter needs
# those modules to attribute their docstrings. With only [Epicycle] listed, every `@docs`
# entry for a force model or geometry reported "no docs found" and rendered nothing.
using AstroModels, AstroProp, AstroStates, AstroEpochs, AstroFrames, AstroUniverse
using AstroCallbacks, AstroManeuvers, AstroSolve, AstroRoutines, EpicycleBase, EpicycleIO

DocMeta.setdocmeta!(Epicycle, :DocTestSetup, :(using Epicycle); recursive=true)

# One copy of an example: the script in examples/, rendered into a page. EPICYCLE_DOC_EXAMPLES says
# which of them this build runs; see docs/examples.jl.
include("examples.jl")

const EXAMPLE_NAV = example_pages(
    ["Propagation" => ["Propagation and stopping conditions" => "Ex_PropagationBasics",
                       "Impulsive maneuvers" => "Ex_ImpulsiveManeuver",
                       "The spacecraft history" => "Ex_SpacecraftHistory",
                       "Propagation about the Moon" => "Ex_Propagation_Moon",
                       "Propagation about Mars" => "Ex_Propagation_Mars"],
     "Targeting" => ["Targeting a single maneuver" => "Ex_SimpleTarget",
                     "Hohmann transfer" => "Ex_HohmannTransfer",
                     "GEO transfer, as an event sequence" => "Ex_GeoTransfer",
                     "GEO transfer, in flight order" => "Ex_GeoTransferTargetBlock",
                     "Station keeping" => "Ex_StationKeeping"],
     "Optimal control" => ["The brachistochrone" => "Ex_Brachistochrone",
                           "The Hull problem" => "Ex_HullProblem",
                           "Linear tangent steering" => "Ex_LinearTangentSteering",
                           "Obstacle avoidance" => "Ex_ObstacleAvoidance",
                           "Soft lunar landing" => "Ex_MoonLanding",
                           "The Goddard rocket" => "Ex_GoddardRocket",
                           "Low thrust orbit raising" => "Ex_LowThrustOrbitRaising",
                           "Orbit raising with a custom force model" => "Ex_ForceModelCustomForceModel",
                           "Earth-Moon Lyapunov transfer" => "Ex_LyapunovTransfer",
                           "Reference tracking with a path bound" => "Ex_ReferenceTracking",
                           "Parameter identification" => "Ex_ParameterIdentification"],
     "Low thrust and interplanetary" =>
         ["Earth to Mars, Sims-Flanagan" => "Ex_MarsTransferSimsFlanagan",
          "Earth to Apophis rendezvous, Sims-Flanagan" => "Ex_ApophisRendezvousSimsFlanagan",
          "Two transcriptions in one problem" => "Ex_MixedTranscription",
          "Earth-Earth-Venus gravity assist" => "Ex_GravityAssistMGA",
          "Gravity assist with an Epicycle maneuver" => "Ex_GravityAssistMGAManeuver"],
     "Orbit determination" => ["Batch least squares" => "Ex_OrbitDetermination",
                               "Extended Kalman filter" => "Ex_ExtendedKalmanFilter",
                               "Stepping the filter" => "Ex_SteppedKalmanFilter"],
     "Plotting and visualization" => ["Plotting results" => "Ex_PlottingResults",
                                      "Viewing a trajectory in 3D" => "Ex_TrajectoryView"]],
    joinpath(@__DIR__, "..", "examples"),
    joinpath(@__DIR__, "src", "examples"))

makedocs(;
    modules=[Epicycle, EpicycleBase, AstroStates, AstroEpochs, AstroUniverse, AstroFrames,
             AstroModels, AstroManeuvers, AstroRoutines, AstroCallbacks, AstroProp,
             AstroSolve, EpicycleIO],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="Epicycle.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/",
        edit_link="main",
        assets=String[],
        sidebar_sitename=false,
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Packages" => "packages.md",
        "Examples" => EXAMPLE_NAV,
    ],
    # Each package's doctests run in its own build, with its own DocTestSetup. The
    # umbrella lists every module so it can attribute their docstrings to its pages, and
    # re-running their doctests here would execute each one twice under the wrong setup.
    doctest=false,
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="Epicycle",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
