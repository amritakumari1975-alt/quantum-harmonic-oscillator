# Computational Study of the Quantum Harmonic Oscillator

## Overview
This project numerically solves the one-dimensional quantum harmonic oscillator using the finite-difference method. The numerical energy eigenvalues are compared with the analytical result

E_n = n + 1/2

in dimensionless units. The project also visualizes the first few normalized wave functions and their probability density distributions.

## Objectives
- Discretize the time-independent Schrödinger equation.
- Construct the Hamiltonian matrix using finite differences.
- Calculate the lowest energy eigenvalues and eigenstates.
- Compare numerical results with analytical energy levels.
- Visualize wave functions and probability densities.

## Technologies
- Python
- NumPy
- SciPy
- Matplotlib

## Installation
```bash
pip install -r requirements.txt
```

## Run
```bash
python quantum_harmonic_oscillator.py
```

The program creates:
- `energy_comparison.png` — numerical vs analytical energy levels
- `wavefunctions.png` — first five normalized wave functions
- `probability_densities.png` — corresponding probability densities

## Method
In dimensionless units (ℏ = m = ω = 1), the Hamiltonian is

H = -1/2 d²/dx² + 1/2 x²

The second derivative is approximated with a finite-difference scheme and the resulting Hamiltonian is diagonalized numerically.

## Expected Result
The lowest numerical energy levels should be close to:

n = 0 → 0.5
n = 1 → 1.5
n = 2 → 2.5
n = 3 → 3.5
n = 4 → 4.5

Small differences occur because the continuous problem is represented on a finite numerical grid.

## Sample Output

```text
Quantum Harmonic Oscillator
---------------------------------------------
n    Numerical         Analytical        Absolute Error 
0    0.49999444        0.50000000        5.565e-06      
1    1.49997218        1.50000000        2.782e-05      
2    2.49992766        2.50000000        7.234e-05      
3    3.49986087        3.50000000        1.391e-04      
4    4.49977183        4.50000000        2.282e-04      

Plots saved:
  energy_comparison.png
  wavefunctions.png
  probability_densities.png
```
