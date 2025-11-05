# Introduction

Welcome to Epicycle, an application for space mission analysis, trajectory optimization, and navigation. 

Epicycle provides an integrated application for astrodynamics analysis with focus on breadth and extensible interfaces. The system covers orbital state representations, coordinate transformations, trajectory propagation, and optimization.

The system handles Cartesian, Keplerian, and Modified Equinoctial orbital elements. Trajectory propagation uses Julia's differential equation solvers. Optimization connects SNOW algorithms with IPOPT.

Current implementation emphasizes application integration over individual model depth. Interfaces are designed for systematic expansion of the model library.