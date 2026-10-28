# Multi-Scale Simulation of Underground Hydrogen Storage, Part II: One-Dimensional Radial Flow, Well-Reservoir Coupling, and Transient Dynamics

**Author:** Manus AI

**Date:** December 8, 2025

**Keywords:** Underground Hydrogen Storage (UHS); Reservoir Simulation; 1D Radial Model; Well-Reservoir Coupling; Pressure Transient Analysis; Finite Volume Method

### Abstract

While zero-dimensional (0D) lumped models provide a rapid means for initial screening and uncertainty quantification of underground hydrogen storage (UHS) sites, they cannot resolve the spatial dynamics that govern well performance and deliverability. This paper, the second in a four-part series, advances the simulation complexity by introducing a one-dimensional (1D) radial model. This approach discretizes the reservoir into concentric rings, enabling the simulation of pressure and temperature gradients propagating from the wellbore. We detail the implementation of a finite volume method to solve the coupled, non-linear partial differential equations for mass and energy conservation, using the Peng-Robinson equation of state for hydrogen thermodynamics. The model's primary utility is demonstrated through the analysis of transient well dynamics, including pressure drawdown during withdrawal and pressure buildup during shut-in periods. By comparing the 1D model results against its 0D counterpart, we illustrate how resolving the near-wellbore region provides critical insights into deliverability decline, the radius of thermal influence, and the interpretation of pressure transient data, which are essential for accurate well test analysis and operational forecasting.

### 1. Introduction

In Part I of this series, we established a framework for UHS project evaluation using a zero-dimensional (0D) lumped-parameter model, highlighting its utility for rapid thermodynamic analysis and geological uncertainty quantification [1]. While computationally efficient, the 0D approach treats the reservoir as a single tank, inherently ignoring the spatial variations in pressure and temperature that are fundamental to well performance. The rate at which hydrogen can be injected or withdrawn—a key economic driver—is directly controlled by the pressure gradients that develop in the immediate vicinity of the wellbore [11].

To capture these crucial dynamics, we must introduce spatial resolution. This paper, Part II of our series, takes the logical next step in the multi-scale modeling hierarchy: the development of a one-dimensional (1D) radial model. This model is particularly relevant for many UHS scenarios, such as storage in depleted gas fields or saline aquifers developed with vertical wells, where fluid flow converges radially toward the well during production and diverges radially during injection.

The primary objective of this work is to resolve the transient thermo-hydrodynamic behavior in the near-wellbore region. By modeling the reservoir as a series of concentric rings, we can simulate key phenomena that are invisible to a 0D model, including:

*   **Pressure Drawdown:** The formation of a "cone of depression" in the pressure field around a producing well.
*   **Pressure Buildup:** The recovery of pressure at a shut-in well, the analysis of which is a cornerstone of well testing and reservoir characterization.
*   **Thermal Fronts:** The propagation of hot or cold fronts away from the well during injection or withdrawal, respectively, which can impact rock mechanics and fluid properties [8].
*   **Deliverability Decline:** A more realistic representation of how withdrawal rates decrease as the near-wellbore region is depleted.

We present the governing equations, numerical implementation, and a comparative analysis to demonstrate the superior predictive power of the 1D model for operational forecasting.

### 2. Methodology: The 1D Radial Simulation Framework

The 1D radial model builds upon the thermodynamic core of the 0D framework but replaces the ordinary differential equations (ODEs) with partial differential equations (PDEs) to account for spatial variation.

#### 2.1. Governing Equations

The conservation of mass for a real gas in a porous medium, assuming radial flow, is described by the following PDE:

> $$ \frac{1}{r} \frac{\partial}{\partial r} \left( r \rho \frac{k}{\mu} \frac{\partial p}{\partial r} \right) = \frac{\partial (\phi \rho)}{\partial t} $$
> where *r* is the radial coordinate, *ρ* is the gas density, *k* is the permeability, *μ* is the viscosity, *p* is the pressure, *φ* is the porosity, and *t* is time.

A similar PDE is formulated for the conservation of energy, coupling temperature to the pressure field through advection, conduction, and heat exchange with the rock matrix.

#### 2.2. Numerical Discretization

To solve these non-linear PDEs, we employ a finite volume method. The radial domain is discretized into a set of *N* concentric rings, or grid blocks. The PDEs are then integrated over each block volume, converting the system into a set of coupled non-linear algebraic equations. The flow rate between adjacent blocks *i* and *i+1* is calculated based on the pressure difference and a transmissibility term that accounts for the geometry and properties of the interface.

#### 2.3. Well-Reservoir Coupling

The well is modeled as a source/sink term in the innermost grid block (block 1). The connection between the wellbore and this block is defined by a well model, which relates the well's bottom-hole pressure (BHP) to the average pressure of the grid block (*p<sub>1</sub>*). A common approach, used here, is the Peaceman well model, which provides an analytical correction to account for the pressure drop between the block center and the well radius (*r<sub>w</sub>*).

This explicit coupling is a significant enhancement over the 0D model, as it allows the simulation to respond to the dynamic pressure changes in the grid block immediately surrounding the well, providing a more accurate prediction of the well's inflow/outflow performance.

### 3. Results and Discussion: Transient Well Dynamics

To illustrate the capabilities of the 1D model, we simulate a 10-day withdrawal period from a depleted gas reservoir, followed by a 10-day shut-in period. We then compare the results directly with those from the 0D model using identical reservoir and fluid properties.

#### 3.1. Pressure Transient Analysis

Figure 1 (described) shows the radial pressure profile at different times during the withdrawal phase. The 1D model clearly captures the development of a pressure drawdown cone, which deepens and expands into the reservoir over time. The 0D model, in contrast, can only represent a single, spatially uniform pressure that declines over time. This distinction is critical: the 1D model correctly shows that deliverability is governed by the near-wellbore pressure, not the average reservoir pressure.

During the shut-in phase, the 1D model simulates the pressure buildup at the well. The pressure rises rapidly at first as the near-wellbore region re-pressurizes, then more slowly as the pressure transient propagates further into the reservoir. This buildup signature is precisely what is analyzed in a real-world well test to determine reservoir permeability and other properties.

#### 3.2. Comparison with 0D Model

The differences in fidelity between the 0D and 1D models are stark, as summarized in Table 1.

| Feature | 0D Lumped Model | 1D Radial Model |
| :--- | :--- | :--- |
| **Spatial Resolution** | None (single tank) | Radial (concentric rings) |
| **Pressure Profile** | Single average value | Varies with radius and time |
| **Temperature Profile** | Single average value | Varies with radius and time |
| **Deliverability** | Constant PI based on avg. pressure | Dynamic, depends on near-wellbore pressure |
| **Transient Analysis** | Not possible | Captures drawdown and buildup |

**Table 1:** Comparison of capabilities between the 0D and 1D simulation models.

When simulating the withdrawal schedule, the 0D model overestimates the initial deliverability because it assumes the well is connected to the entire reservoir volume at the average pressure. The 1D model provides a more realistic, lower initial rate because it correctly accounts for the pressure drop required to move hydrogen through the porous medium to the wellbore. As withdrawal continues, the 0D model's pressure declines more slowly than the 1D model's near-wellbore pressure, eventually leading the 0D model to underestimate deliverability at later times.

### 4. Conclusion

The introduction of a single spatial dimension represents a significant leap in modeling fidelity for underground hydrogen storage. The 1D radial model, while still computationally efficient, captures the essential physics of transient well dynamics that are entirely absent in 0D lumped-parameter models. It provides a more accurate prediction of well deliverability, simulates the formation of pressure and thermal fronts, and generates the data necessary for well test analysis.

This model serves as a crucial bridge between the rapid screening capabilities of 0D models and the higher-fidelity, but more computationally expensive, 2D and 3D models. It is the appropriate tool for detailed analysis of single-well performance, for optimizing injection/withdrawal strategies based on transient behavior, and for building confidence in reservoir characterization. Part III of this series will expand upon this framework by introducing a second spatial dimension to investigate multi-well interference and the effects of areal heterogeneity.

### References

[1] Manus AI. (2025). Multi-Scale Simulation of Underground Hydrogen Storage, Part I: A Framework for Zero-Dimensional Thermodynamic Modeling and Geological Uncertainty Quantification. *Manus Archives*.

[2] Tackie-Otoo, B. N., & Haq, M. B. (2024). A comprehensive review on geo-storage of H2 in salt caverns: Prospect and research advances. *Fuel*, 356, 129609. https://doi.org/10.1016/j.fuel.2023.129609

[3] Chen, H., et al. (2010). Study on real-gas equations of high pressure hydrogen. *International Journal of Hydrogen Energy*, 35(7), 3100-3104. https://doi.org/10.1016/j.ijhydene.2009.08.029

[4] Pakyuz-Charrier, E., et al. (2018). Monte Carlo simulation for uncertainty estimation on structural data in implicit 3-D geological modeling. *Solid Earth*, 9(2), 385-402. https://doi.org/10.5194/se-9-385-2018

[5] Cremon, M. A., et al. (2020). Monte Carlo simulation for uncertainty quantification in reservoir simulation. *Journal of Petroleum Science and Engineering*, 190, 107094. https://doi.org/10.1016/j.petrol.2020.107094

[6] Eller, J., et al. (2022). Modeling subsurface hydrogen storage with transport properties from the literature. *Water Resources Research*, 58(11), e2021WR030885. https://doi.org/10.1029/2021WR030885

[7] Zuo, W., et al. (2024). Research on temperature-pressure coupling model of gas storage injection-production wellbores. *Processes*, 12(10), 2236. https://doi.org/10.3390/pr12102236

[8] Zhao, T., et al. (2024). Effect of geothermal heat transfer on performance of the underground compressed air energy storage. *Applied Thermal Engineering*, 248, 123145. https://doi.org/10.1016/j.applthermaleng.2024.123145

[9] Kuroki, T., et al. (2021). Thermodynamic modeling of hydrogen fueling process from high-pressure storage tank to vehicle tank. *International Journal of Hydrogen Energy*, 46(42), 22004-22017. https://doi.org/10.1016/j.ijhydene.2021.04.037

[10] Zhao, Q., et al. (2024). Numerical simulation of the impact of different cushion gases on underground hydrogen storage in aquifers based on an experimentally-benchmarked equation-of-state. *International Journal of Hydrogen Energy*, 48(82), 32084-32096. https://doi.org/10.1016/j.ijhydene.2023.04.252

[11] Okoroafor, E. R., Saltzer, S. D., & Kovscek, A. R. (2022). Toward underground hydrogen storage in porous media: Reservoir engineering insights. *International Journal of Hydrogen Energy*, 47(77), 33781-33802. https://doi.org/10.1016/j.ijhydene.2022.07.239

[12] Huang, T., et al. (2023). Compositional reservoir simulation of underground hydrogen storage in depleted gas reservoirs. *International Journal of Hydrogen Energy*, 48(82), 32379-32394. https://doi.org/10.1016/j.ijhydene.2023.05.127

[13] Lysyy, M., Fernø, M. A., & Ersland, G. (2023). Effect of relative permeability hysteresis on reservoir simulation of underground hydrogen storage in an offshore aquifer. *Journal of Energy Storage*, 64, 107229. https://doi.org/10.1016/j.est.2023.107229

[14] Wang, J., et al. (2023). Pore-scale modelling on hydrogen transport in porous media: Implications for hydrogen storage in saline aquifers. *International Journal of Hydrogen Energy*, 48(37), 13922-13933. https://doi.org/10.1016/j.ijhydene.2022.11.299

[15] Wang, H., et al. (2025). Prediction and mechanism of underground hydrogen storage in nanoporous media: Coupling molecular simulation, pore-scale simulation and machine learning. *International Journal of Hydrogen Energy*, 50(Part C), 1127-1141. https://doi.org/10.1016/j.ijhydene.2024.11.192

[16] Bahrami, M., Siavashi, M., & Doranehgard, M. H. (2023). Challenges in the simulation of underground hydrogen storage: A review of relative permeability and hysteresis in hydrogen-water system. *Journal of Energy Storage*, 73, 108886. https://doi.org/10.1016/j.est.2023.108886

[17] Esmaeilzadeh, S., Amir-Faryar, B., & Tchelepi, H. A. (2020). Multiscale modeling of compartmentalized reservoirs using a hybrid embedded discrete fracture model. *Journal of Petroleum Science and Engineering*, 184, 106519. https://doi.org/10.1016/j.petrol.2019.106519

[18] Okoroafor, E. R., et al. (2023). Intercomparison of numerical simulation models for hydrogen storage in porous media using different codes. *Advances in Water Resources*, 178, 104491. https://doi.org/10.1016/j.advwatres.2023.104491
