```@meta
CurrentModule = AstroSolve
```

# How to Think in AstroSolve

AstroSolve uses the same model for parameter optimization, optimal control, and orbit estimation. A
problem contains a sequence of discrete and interval events, variables the solver may change, and
constraints that must hold. An optimization problem may include an objective to minimize or
maximize. An estimation problem includes measurement models. This page introduces these concepts
before the later guides apply them to complete problems.

## The trajectory graph

AstroSolve models optimization and estimation problems as sequences of events arranged in a
directed acyclic graph (DAG). Nodes are discrete events, such as an impulsive maneuver. Edges are
interval events, such as propagation or a finite maneuver. The direction of each connection defines
the event sequence. A graph may be a simple sequence, or it may branch and merge.

The graph separates mission structure from the method used to solve each interval. One interval
may be an ordinary propagation, another a collocation phase, and another a shooting phase. The
events connecting them may apply maneuvers, introduce variables, or evaluate constraints.

Parameter optimization, optimal control, and estimation use this structure differently, but share
the same formulation wherever practical.

## A transfer as a graph

A transfer with a departure maneuver, a coast to an intermediate correction, a thrust phase, and
an arrival condition shows the sequence structure:

```text
 departure              correction                              arrival
    event                   event                                  event
      ●──────── coast ─────────●───────────── thrust phase ──────────●
   vary Δv             vary correction                       constrain orbit
                                                               minimize cost
```

The departure and correction maneuvers are discrete events. The coast and thrust phase are interval
events. The maneuver components are variables, the arrival orbit is constrained, and the thrust phase may
carry an objective such as propellant use or final mass. The same graph can mix propagation and
optimal control because the connections describe the trajectory independently of the method used
on each interval.

The rest of this page describes the elements in a sequence: events, phases, variables, constraints,
and objectives. The later guides show how to use them to solve problems.

## Sequences, events, intervals, and phases

Optimization and estimation problems are modeled as sequences of discrete and interval events. A
discrete event acts at one point in the trajectory, such as an impulsive maneuver. An interval event
has a finite duration, such as a propagation that stops at periapsis or a finite maneuver.

A phase is an interval event represented by a transcription. It contains the dynamics, time span,
state, control, and any static parameters needed for that part of the trajectory. A transcription
may also introduce its own variables and constraints; Sims-Flanagan is one example.

## Variables, constraints, and objectives

`Vary`, `Constraint`, and `Objective` turn a trajectory or estimation model into a solver problem.
`Vary` identifies values the solver may change, `Constraint` sets conditions those values must
satisfy, and `Objective` supplies the scalar value to minimize or maximize. Variables and
constraints may belong to an event, a phase, or a link between phases. Bounds define the permitted
values; scales keep unlike quantities numerically comparable. In estimation, a varied quantity has
a prior covariance rather than optimization bounds.

## State flow and continuity

Each event receives the state left by the preceding event. A maneuver changes that state, and a
propagation advances it to the next event. A phase enforces state flow internally through its
transcription. Continuity between phases is separate and must be imposed on their `Link`; it is not
inferred from their order in the sequence. The standard continuity constraint matches boundary
state and time. A custom link constraint represents a deliberate discontinuity such as staging or
a flyby.

## Transcriptions

A transcription turns a continuous phase into the finite variables and equations passed to the
solver. Collocation represents state and control at mesh points and constrains the resulting
trajectory to satisfy the dynamics. Shooting propagates arcs and constrains their endpoints or
match points. The choice changes the numerical problem, not the mission definition: variables,
constraints, objectives, and links use the same vocabulary across transcriptions. This also allows
different transcription methods to appear in one sequence.

## Derivatives

AstroSolve needs derivatives of objectives, constraints, and dynamics. It uses derivatives declared
with `@partial` where provided and automatic differentiation where supported. Event-sequence
parameter optimization can use finite differences. Functions differentiated automatically must not
assume `Float64` inputs. `check_partials` compares declared derivatives with finite differences
before a solve.

## References

- Hughes, S. P. (2026), "A Transcription-Agnostic Formulation for Optimal Control and Estimation
  in Astrodynamics," AAS/AIAA Astrodynamics Specialist Conference, Vancouver, British Columbia,
  July 2026.
- Williams, J., Falck, R., and Beekman, I. (2018), "Application of Modern Fortran to Spacecraft
  Trajectory Design and Optimization," AIAA/AAS Space Flight Mechanics Meeting.
  [Online](https://ntrs.nasa.gov/api/citations/20180000413/downloads/20180000413.pdf)
