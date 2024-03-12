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

Earthquake cycle simulations are important as they provide more understanding of the behavior along seismic faults. We can study the effects of different boundary conditions, friction parameters, and slip rates and how they affect material displacement over time. 
[@Erickson19]

# Statement of need
(A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.)

Thrase.jl provides a high-order accurate SBP-SAT finite difference code on unstructured meshes for SEAS (Sequences of Earthquakes and Aseismic Slip) problems. It is written in Julia and contains both a serial and GPU-enabled parallel (TODO specific sections enabled?) version. Thrase.jl solves the wave equation (Poisson) in second order form by tracking the evolution of boundary and interface displacements [@Erickson19]. This allows for a non-stiff character formulation and use of explicit time-stepping methods such as Runge-Kutta by rewriting the nonlinear friction laws with the method in [Kozdon2012]. This contrasts with previous work by [@Duru2019], [@Virta2014], which are limited by the stiff character formulation. 

# Acknowledgements
(Any funding acknowledgements?)

This study is supported by...

# References
