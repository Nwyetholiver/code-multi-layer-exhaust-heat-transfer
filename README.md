# Multi-Layer Exhaust Heat Transfer Analysis

## Objective

Develop a transient heat transfer model for a multi-layer cylindrical exhaust system and validate the thermal behavior through numerical simulation and real-world experimentation.

## Overview

This project models radial heat conduction through a three-layer cylindrical wall representing an automotive exhaust system. The system includes steel, insulation, and aluminum layers. A transient finite difference method (FDM) was implemented to analyze temperature distribution over time and evaluate thermal performance.

To strengthen the analysis, an experimental setup using a 3D printed muffler prototype was used to observe real temperature trends and compare them with simulation results.

## Engineering Problem

Automotive exhaust systems operate under high thermal loads. Managing heat transfer is critical for:

* Preventing material failure
* Improving efficiency
* Reducing external surface temperatures for safety

This project investigates how layered materials affect heat dissipation and insulation performance.

## Tools & Methods

* MATLAB (Finite Difference Method - Explicit)
* SolidWorks (geometry modeling)
* 3D Printing (prototype muffler)
* IR Thermometer (temperature measurement)
* Heat Transfer Theory (conduction + convection)

## Model Description

* 1D radial heat conduction in cylindrical coordinates
* No internal heat generation
* Multi-layer system:

  * Steel (inner layer)
  * Insulation (middle layer)
  * Aluminum (outer layer)

## Governing Equation

The transient heat conduction equation in cylindrical coordinates:

∂T/∂t = α (1/r ∂/∂r (r ∂T/∂r))

Boundary Conditions:

* Inner surface: constant high temperature (exhaust gas)
* Outer surface: convection to ambient air

## Numerical Method

* Explicit finite difference scheme
* Stability controlled using Fourier Number (Fo)

Two cases analyzed:

* Stable case (Fo < 0.5)
* Unstable case (Fo > 0.5)

## Results

### Key Observations

* Temperature drops significantly across the insulation layer
* Steeper gradients indicate higher heat flux regions
* System approaches steady-state as time increases
* Instability occurs when time step is too large

### Validation

* Transient results converge toward steady-state solution
* Experimental trends align with simulation behavior
* Differences attributed to convection assumptions and material uncertainty

### Sensitivity Analysis

* Smaller Δr improves spatial accuracy
* Smaller Δt improves stability and smoothness
* Coarser grids lead to numerical error

## Experimental Setup

* 3D printed muffler prototype
* Attached to idling vehicle exhaust
* Temperature measured over time at outer surface
* Compared against simulation predictions

## Key Takeaways

* Multi-layer insulation drastically reduces outer surface temperature
* Explicit FDM requires strict stability criteria
* Numerical modeling can effectively predict real thermal behavior
* Combining simulation with experimentation strengthens engineering validation

## Future Improvements

* Include radiation heat transfer
* Use implicit methods for better stability
* Improve convection modeling with CFD
* Use higher accuracy sensors (thermocouples)

## Repository Structure

/code → MATLAB scripts
/data → experimental temperature data
/images → plots, diagrams, setup photos
/cad → SolidWorks models
/report → full project report

## Author

Noah Oliver
Mechanical Engineering Student
