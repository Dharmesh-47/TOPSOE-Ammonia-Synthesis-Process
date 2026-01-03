# TOPSOE Ammonia Process – Process Modeling & Kinetic Analysis

A comprehensive process modeling and analysis of the Haldor Topsoe ammonia synthesis process, integrating separation design, thermodynamic equilibrium analysis, and modern reaction kinetics.

This repository contains MATLAB-based simulations and datasets used to analyze key unit operations (letdown vessel, distillation column, compressor) and to study ammonia synthesis efficiency using equilibrium thermodynamics and the Extended Stoltze kinetic model.

<pre>
TOPSOE-Ammonia-Process/
│── Group-1 final_report_Topsoe.pdf     # Final project report
│── Group 1-new.xlsx                    # Process & operating data
│── mc_cabe.m                           # McCabe–Thiele distillation analysis
│── Mccabe_new.m.txt                    # Updated distillation calculations
│── paper_application.m                 # Kinetics & efficiency analysis
│── README.md                           # Documentation
</pre>

## Project Scope

- Process understanding of the Haldor Topsoe ammonia synthesis loop
- Separation system analysis using flash calculations and distillation design
- Equilibrium-based ammonia yield estimation
- Efficiency analysis incorporating modern kinetic insights
- Validation of trends reported in recent ammonia synthesis literature

## Tools & Methods

- MATLAB for numerical modeling and visualization
- Chemical thermodynamics (equilibrium, fugacity concepts)
- Mass and energy balance formulations
- McCabe–Thiele graphical distillation design
- Literature-backed kinetic modeling (Extended Stoltze formulation)

## Features

### Separation & Distillation Analysis

- Flash and letdown vessel modeling under adiabatic conditions.
- Component-wise material and energy balances.
- McCabe–Thiele method to determine:
  - Minimum number of stages
  - Operating reflux ratio
  - Optimal feed stage location
- MATLAB implementation for reproducible column design.

### Equilibrium Yield Estimation

- Calculation of equilibrium ammonia mole fraction using thermodynamic relations.
- Temperature and pressure dependence of equilibrium yield.
- Comparison between theoretical equilibrium and actual outlet composition.

### Ammonia Synthesis Efficiency Analysis

- Definition of efficiency as:
  
  η = y_out / y_eq

- Parametric studies over:
  - Temperature
  - Pressure
  - H₂:N₂ feed ratio
- Identification of kinetically optimal operating windows.

### Extended Stoltze Kinetic Model Application

- Incorporation of temperature- and hydrogen-dependent active-site density.
- Analysis of hydrogen inhibition effects at low temperatures.
- Correlation between catalyst site density and ammonia synthesis efficiency.
- Validation of trends against reported experimental literature.

### Visualization & Results

- Efficiency vs. temperature plots.
- Efficiency vs. H₂:N₂ ratio plots.
- Efficiency vs. predicted site density plots.
- Parity plots comparing equilibrium and outlet ammonia compositions.
- MATLAB figure files retained for reproducibility.

## Repository Notes

- MATLAB scripts are written for clarity and direct interpretability.
- Excel files contain both raw and processed data used in simulations.
- The final report provides full theoretical background, assumptions, derivations, and discussion of results.

## References

- Nadiri, S. et al. (2024). *Ammonia synthesis rate over a wide operating range: From experiments to validated kinetic models*. ChemCatChem.
- Stoltze, P., & Nørskov, J. K. (1985). Bridging the pressure gap in ammonia synthesis kinetics.
- McCabe, W. L., Smith, J. C., & Harriott, P. *Unit Operations of Chemical Engineering*.
- Haldor Topsoe technical literature and process documentation.

---

*This project was developed as part of an Applied Process Engineering course and extended to include modern kinetic modeling and literature-backed analysis of industrial ammonia synthesis.*
