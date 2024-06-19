---
title: 'Thrase.jl: A high-order accurate SBP-SAT finite difference code on unstructured meshes for SEAS problems'
tags:
  - fault dynamics
  - friction
  - geophysics
  - Julia
authors:
  - name: Brittany Erickson
    orcid: 0000-0001-9457-8572
    affiliation: 1
  - name: Dewi Yokelson
    orcid: 0000-0003-1453-5906
    affiliation: 1
  - name: Alexandre Chen
    orcid:
    affiliation: 1
affiliations:
 - name: University of Oregon
   index: 1
date: 30 March 2024
bibliography: paper.bib
---

(The whole paper should be 250-1000 words)

# Summary
(A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.)

Earthquake cycle simulations are important as they provide more understanding of the behavior along seismic faults. These simulations can help predict not only when an earthquake might occur, but the intensity as well. However, earthquake cycles are complex to model, as they depend on many different dynamic parameters that change over time. Translating these problems so that they can be solved computationally enables much faster time to solution but requires using numerical methods to discretize and model them accurately. Thrase.jl provides such a mathematical framework for numerically determining how these parameters interact and simulating many different earthquake cycle problems. We can study the effects of different boundary conditions, friction parameters, and slip rates and how they affect material displacement over time. There are three example benchmark problems, and the code can be modified and extended to solve variations or new problems entirely.
[@Erickson19]

# Statement of need
(A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.)

(TODO: include an example figure, maybe from BP1)

Thrase.jl provides a high-order accurate SBP-SAT finite difference code on unstructured meshes for SEAS (Sequences of Earthquakes and Aseismic Slip) problems. It is written in Julia and contains both serial and GPU-enabled parallel capabilities. Thrase.jl solves the wave equation (Poisson) in second order form by tracking the evolution of boundary and interface displacements [@Erickson19]. This allows for a non-stiff character formulation and use of explicit time-stepping methods such as Runge-Kutta by rewriting the nonlinear friction laws with the method in [Kozdon2012]. This contrasts with previous work by [@Duru2019], [@Virta2014], which are limited by the stiff character formulation. 

# Acknowledgements
(Any funding acknowledgements?)

This study is supported by...

# References
