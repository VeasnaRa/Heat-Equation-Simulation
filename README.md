# Heat-Equation-Simulation-C++
A C++ project that simulates heat diffusion in 1D (bar) and 2D (plate) implementing the finite difference method. The simulation uses SDL2 for real-time visualization with a heatmap display.

## Overview

This project simulates heat diffusion through materials using the heat equation, a fundamental partial differential equation (PDE) in physics. The simulator supports both 1D (bar) and 2D (plate) scenarios with real-time visualization of temperature distribution across four different materials simultaneously.

**Key Features:**
- Unconditionally stable implicit scheme (Backward Euler)
- 1D Bar Simulation: Linear heat diffusion with 2 heat sources
- 2D Plate Simulation: Radial heat diffusion with 4 corner sources
- Real-time visualization at 60 FPS using SDL2
- Simultaneous 4-material comparison in 2×2 grid mode
- Complete Doxygen documentation with UML class diagram

---

## Mathematical Foundation

### The Heat Equation

The heat equation models temperature distribution over time:

$$\frac{\partial u}{\partial t} = \frac{\lambda}{\rho c} \Delta u + \frac{F}{\rho c}$$

Where:
- **u(t,x)**: Temperature [K] at time t and position x
- **α = λ/(ρc)**: Thermal diffusivity [m²/s]
- **λ**: Thermal conductivity [W/(m·K)]
- **ρ**: Density [kg/m³]
- **c**: Specific heat capacity [J/(kg·K)]
- **F**: Heat source [W/m³]

### Material Properties
| Material | $\lambda$ (W/(m·K)) | $\rho$ (kg/m³) | $c$ (J/(kg·K)) |
|----------|---------------------|----------------|----------------|
| Copper | 389 | 8940 | 380 |
| Iron | 80.2 | 7874 | 440 |
| Glass | 1.2 | 2530 | 840 |
| Polystyrene | 0.1 | 1040 | 1200 |

### Boundary Conditions

- **Neumann** (x=0, y=0): ∂u/∂n = 0 (insulated boundary)
- **Dirichlet** (x=L, y=L): u = u₀ (fixed temperature)

---

## Numerical Methods

### 1D Case: Implicit Finite Differences + Thomas Algorithm

**Mathematical preliminaries**

$$\frac{\partial u}{\partial t} = \frac{\lambda}{\rho c}\frac{\partial^2 u}{\partial x^2} + \frac{F}{\rho c}$$

#### Discretization

**Temporal discretization:**

$$\Delta t = \frac{t_{max}}{1000}, \quad t_n = n \cdot \Delta t$$

**Spatial discretization:**

$$\Delta x = \frac{L}{N-1}, \quad x_i = i \cdot \Delta x$$


#### Implicit Scheme (Backward Euler)

The backward Euler method discretizes the heat equation as:

$$-r \cdot u_{i-1}^{n+1} + (1 + 2r) \cdot u_i^{n+1} - r \cdot u_{i+1}^{n+1} = u_i^n + \Delta t \frac{F_i}{\rho c}$$

where **r = α·Δt/Δx²** is the diffusion number.

#### Thomas Algorithm (TDMA)

Solves tridiagonal systems in **O(n)** complexity:

For a tridiagonal system of form:    $$a_i u_{i-1} + b_i u_i + c_i u_{i+1} = d_i$$

**Forward sweep:**

$$c'_i = \frac{c_i}{b_i - a_i c'_{i-1}}, \quad d'_i = \frac{d_i - a_i d'_{i-1}}{b_i - a_i c'_{i-1}}$$

**Back substitution:**

$$u_{n-1} = d'_{n-1}, \quad u_i = d'_i - c'_i u_{i+1}$$

**Complexity:** This forms a **tridiagonal system** $A \mathbf{u}^{n+1} = \mathbf{b}$ solved by the **Thomas algorithm** in $O(n)$.

### 2D Case: 5-Point Stencil + Gauss-Seidel Iteration
**Mathematical preliminaries**

$$\frac{\partial u}{\partial t} = \alpha \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) + \frac{F}{\rho c}$$

#### Discretized Form

Using a 5-point stencil, each grid point depends only on its 4 neighbors:

$$-r \cdot u_{i-1,j}^{n+1} - r \cdot u_{i,j-1}^{n+1} + (1 + 4r) \cdot u_{i,j}^{n+1} - r \cdot u_{i+1,j}^{n+1} - r \cdot u_{i,j+1}^{n+1} = u_{i,j}^n + \Delta t\frac{F_{i,j}}{\rho c}$$

#### Gauss-Seidel Method

**Update formula:**

$$u_{i,j}^{new} = \frac{1}{1 + 4r} \left[ u_{i,j}^n + \Delta t \frac{F_{i,j}}{\rho c} + r(u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1}) \right]$$

**Convergence criterion:**

$$\max_{i,j} |u_{i,j}^{(k+1)} - u_{i,j}^{(k)}| < \varepsilon = 10^{-6}$$

Typically converges in less than 20 iterations for the heat equation.

---

## Installation

### Prerequisites

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install g++ libsdl2-dev doxygen
```

**macOS:**
```bash
brew install sdl2 doxygen
```

**Verify installation:**
```bash
pkg-config --modversion sdl2
```

### Build
```bash
# Clone repository
git clone https://github.com/VeasnaRa/Heat-Equation-Simulation.git
cd heat-equation-simulator

# Compile
g++ -g -Wall -Wextra -o heat_sim *.cpp $(pkg-config --cflags --libs sdl2)

# Run the Simulator
./heat_sim
```
--- 

## Usage
### Interactive Menu
```
========================================
   HEAT EQUATION SIMULATOR
   ENSIIE - Master 1
========================================

SELECT SIMULATION TYPE
----------------------
  1. 1D Bar  (All 4 Materials - 2x2 Grid)
  2. 2D Plate (All 4 Materials - 2x2 Grid)
  0. Quit
Choice: 1

PARAMETERS (Enter for default, 'b' to go back)
----------------------------------------------
Domain length L [1.0] m: 
Max time tmax [16.0] s: 
Initial temp u0 [13.0] C: 
Source amplitude f [80.0] C: 
```

### Keyboard Controls

| Key | Action |
|-----|--------|
| `SPACE` | Pause/Resume simulation |
| `R` | Reset to initial state |
| `↑` / `↓` | Increase/Decrease simulation speed |
| `ESC` | Quit simulation |

---

## Project Structure
```
heat-equation-simulator/
├── heat_equation_solver.hpp/cpp  # Numerical solvers (1D/2D)
├── material.hpp                  # Material properties
├── sdl_core.hpp/cpp              # SDL initialization
├── sdl_window.hpp/cpp            # Window management
├── sdl_heatmap.hpp/cpp           # Visualization engine
├── sdl_app.hpp/cpp               # Application controller
├── main.cpp                      # Entry point & menu
├── Doxyfile                      # Documentation config
├── UML_diagram.puml              # Class diagram source
├── rapport_PAP.pdf               # Report detail about the project
└── README.md                     # This file
```

---

## Documentation

### Generate Documentation
```bash
# Generate Doxygen documentation
doxygen Doxyfile

# View documentation in browser
xdg-open doc/html/index.html   # Linux
open doc/html/index.html       # macOS
```

### Algorithm Complexity

| Case | Method | Complexity | Speedup |
|------|--------|------------|---------|
| 1D | Thomas Algorithm | O(n) | 1000× vs O(n³) |
| 2D | Gauss-Seidel | O(k·n²) | 1,000,000× vs O(n⁶) |

## References

### Numerical Methods
- [Tridiagonal Matrix Algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
- [Gauss-Seidel Method (Wikipedia)](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
- [Finite Difference Method (Wikipedia)](https://en.wikipedia.org/wiki/Finite_difference_method)

### Libraries
- [SDL2 - Simple DirectMedia Layer](https://www.libsdl.org/)
- [Doxygen - Documentation Generator](https://www.doxygen.nl/)
