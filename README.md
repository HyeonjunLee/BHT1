아래는 해당 코드를 설명하는 README 파일의 예시입니다. 필요에 따라 내용을 수정하여 사용하실 수 있습니다.

---

# Transient Heat Conduction in a Composite Cylinder with Variable Material Properties

This project simulates the transient heat conduction in a composite cylindrical barrel that consists of a thin chrome coating on the inner wall and a steel substrate. The simulation also incorporates the effect of convective heat transfer at the inner surface with high-temperature internal air, which is modeled using a lumped parameter approach.

## Overview

The simulation calculates the temperature evolution in both the composite barrel (chrome + steel) and the internal air over time by solving the transient heat conduction equation in cylindrical coordinates using an explicit finite difference method. The material properties (thermal conductivity and specific heat) for chrome and steel are modeled as temperature-dependent through simple linear approximations.

### Key Features

- **Composite Geometry:**  
  - **Chrome Coating:** Thin inner layer (e.g., 0.1 mm thick) with temperature-dependent properties.  
  - **Steel Substrate:** Main barrel material with its own temperature-dependent properties.
  
- **Boundary Conditions:**  
  - **Inner Surface (Chrome):** Convection with internal air (modeled with a convective heat transfer coefficient).  
  - **Outer Surface (Steel):** Dirichlet condition (fixed at ambient temperature).

- **Internal Air Modeling:**  
  - The air inside the barrel is treated as a lumped system. Its temperature evolves based on the convective energy exchange with the inner barrel surface.

- **Variable Material Properties:**  
  - Both thermal conductivity and specific heat for chrome and steel are functions of temperature. Linear approximations are used in this example.

- **Visualization:**  
  - Two animations are provided:  
    1. **Composite Domain Animation:** Displays the entire domain from the center of the air region (r = 0) to the outer steel surface, showing the constant internal air temperature for 0 ≤ r < r₁ and the transient temperature distribution in the barrel for r₁ ≤ r ≤ r₂.  
    2. **Shell-only Animation:** Shows only the barrel (chrome + steel) region.
  - An additional plot shows the evolution of the internal air temperature over time.

## File Structure

- **simulation_code.py** (or as a Jupyter Notebook):  
  Contains the complete code for the simulation and visualization, including:
  - Functions for temperature-dependent material properties:
    - `thermal_conductivity_chrome(T)`
    - `thermal_conductivity_steel(T)`
    - `specific_heat_chrome(T)`
    - `specific_heat_steel(T)`
  - The composite simulation function `transient_heat_conduction_composite_variable_properties()` which computes the transient temperature profiles.
  - Code for generating two animations (composite and shell-only) and a plot of internal air temperature evolution.

## Requirements

- **Python 3**
- **NumPy**
- **Matplotlib**

You can install the required packages using pip:

```bash
pip install numpy matplotlib
```

## How to Run

1. Save the code into a Python file (e.g., `simulation_code.py`) or open it in a Jupyter Notebook.
2. Adjust the simulation parameters as desired (ambient temperature, dimensions, time step, number of grid points, etc.).
3. Run the script. Two animation windows and one static plot will be generated:
   - The first animation shows the entire composite domain (internal air + barrel).
   - The second animation shows only the barrel (shell-only) region.
   - A separate plot displays the evolution of the internal air temperature over time.

## Simulation Parameters (Example)

- **Ambient Temperature (T_outside):** 25 °C  
- **Internal Air Initial Temperature (T_air_initial):** 1000 °C  
- **Inner Diameter:** 155 mm  
- **Steel Thickness:** 77.5 mm  
- **Chrome Coating Thickness:** 0.1 mm  
- **Convective Heat Transfer Coefficient (h_in):** 25 W/m²·K  
- **Total Simulation Time:** 300 s  
- **Time Step (dt):** 0.01 s  
- **Number of Radial Grid Points (N):** 150  
- **Number of Snapshots (num_snapshots):** 50  

Feel free to modify these parameters to match your specific simulation needs.

## Code Explanation

- **Temperature-dependent Functions:**  
  The code uses simple linear models for both chrome and steel to calculate thermal conductivity and specific heat as functions of temperature.

- **Composite Simulation Function:**  
  The function `transient_heat_conduction_composite_variable_properties()` sets up the composite geometry (from the inner surface of the chrome to the outer steel surface) and solves the heat conduction equation using an explicit finite difference method.  
  - It determines which grid points belong to chrome or steel.  
  - It applies a convective boundary condition at the inner surface and a Dirichlet condition at the outer surface.
  - It also updates the internal air temperature using a lumped energy balance.

- **Visualization:**  
  Two animations are created:
  1. **Composite Domain Animation:** Combines the internal air region (from r = 0 to r = r₁, shown as a constant value) and the barrel (r₁ to r₂) into one temperature profile.
  2. **Shell-only Animation:** Displays only the barrel region.
  
  Additionally, a separate plot shows the time evolution of the internal air temperature.

## Acknowledgments

This simulation is intended for educational purposes and demonstrates basic concepts of transient heat conduction in composite materials and the use of numerical methods to solve partial differential equations in cylindrical coordinates.

---

이 README 파일을 기반으로 프로젝트에 대한 개요와 사용법을 문서화할 수 있습니다.
