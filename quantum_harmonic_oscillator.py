import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


def solve_harmonic_oscillator(x_min=-8.0, x_max=8.0, points=1200, states=5):
    # Dimensionless units: hbar = m = omega = 1
    x = np.linspace(x_min, x_max, points)
    dx = x[1] - x[0]

    # Interior grid; Dirichlet boundary condition psi = 0 at both ends
    xi = x[1:-1]
    n = len(xi)

    # Finite-difference approximation:
    # d2 psi/dx2 ≈ (psi[i+1] - 2 psi[i] + psi[i-1]) / dx^2
    kinetic_diag = np.full(n, 1.0 / dx**2)
    kinetic_off = np.full(n - 1, -0.5 / dx**2)

    potential = 0.5 * xi**2

    hamiltonian = diags(
        [kinetic_off, kinetic_diag + potential, kinetic_off],
        offsets=[-1, 0, 1],
        format="csc"
    )

    # Lowest eigenvalues/eigenvectors
    energies, wavefunctions = eigsh(
        hamiltonian,
        k=states,
        which="SA"
    )

    order = np.argsort(energies)
    energies = energies[order]
    wavefunctions = wavefunctions[:, order]

    # Normalize wavefunctions using numerical integration
    for i in range(states):
        norm = np.sqrt(np.trapezoid(wavefunctions[:, i] ** 2, xi))
        wavefunctions[:, i] /= norm

        # Choose a consistent sign for visualization
        if wavefunctions[np.argmax(np.abs(wavefunctions[:, i])), i] < 0:
            wavefunctions[:, i] *= -1

    analytical = np.arange(states) + 0.5

    return xi, energies, analytical, wavefunctions


def create_plots(x, numerical, analytical, wavefunctions):
    states = len(numerical)

    # 1. Energy comparison
    plt.figure(figsize=(8, 5))
    n = np.arange(states)
    plt.plot(n, analytical, "o-", label="Analytical")
    plt.plot(n, numerical, "s--", label="Numerical")
    plt.xlabel("Quantum number n")
    plt.ylabel("Energy")
    plt.title("Quantum Harmonic Oscillator: Energy Comparison")
    plt.xticks(n)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("energy_comparison.png", dpi=300)
    plt.close()

    # 2. Wave functions
    plt.figure(figsize=(9, 6))
    for i in range(states):
        plt.plot(x, wavefunctions[:, i] + numerical[i], label=f"n={i}")

    plt.plot(x, 0.5 * x**2, "k--", alpha=0.5, label="Potential V(x)")
    plt.xlabel("Position x")
    plt.ylabel("Energy / Wavefunction")
    plt.title("First Five Harmonic Oscillator Wave Functions")
    plt.xlim(-5, 5)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("wavefunctions.png", dpi=300)
    plt.close()

    # 3. Probability densities
    plt.figure(figsize=(9, 6))
    for i in range(states):
        plt.plot(x, wavefunctions[:, i] ** 2, label=f"n={i}")

    plt.xlabel("Position x")
    plt.ylabel(r"$|\psi_n(x)|^2$")
    plt.title("Probability Density Distributions")
    plt.xlim(-5, 5)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("probability_densities.png", dpi=300)
    plt.close()


def main():
    x, numerical, analytical, wavefunctions = solve_harmonic_oscillator()

    print("\nQuantum Harmonic Oscillator")
    print("-" * 45)
    print(f"{'n':<5}{'Numerical':<18}{'Analytical':<18}{'Absolute Error':<15}")

    for n, (num, ana) in enumerate(zip(numerical, analytical)):
        error = abs(num - ana)
        print(f"{n:<5}{num:<18.8f}{ana:<18.8f}{error:<15.3e}")

    create_plots(x, numerical, analytical, wavefunctions)

    print("\nPlots saved:")
    print("  energy_comparison.png")
    print("  wavefunctions.png")
    print("  probability_densities.png")


if __name__ == "__main__":
    main()
