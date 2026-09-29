import glob
import os
import re
import numpy as np
from scipy.optimize import curve_fit

# Opens a QE output file and extracts the lattice parameter a and total energy E, and converts the lattice parameter from Bohr to Angstroms, then return both values as Python floats.

def parse_qe_output(filepath):

    """Extracts lattice parameter 'a' (Å) and total energy 'E' (Ry) from QE output."""
    a_val, total_energy = None, None

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:

        for line in f:
            if "celldm(1)" in line:
                match = re.search(r"celldm\(1\)\s*=\s*([0-9\.]+)", line)
                if match:
                    a_val = float(match.group(1)) * 0.52917721092
            if line.strip().startswith("!") and "total energy" in line:
                match = re.search(r"=\s*([-\d\.]+)", line)
                if match:
                    total_energy = float(match.group(1))

    return a_val, total_energy

# The third-order Birch-Murnaghan EOS function, returns internal energy E as a function of a bunch of different parameters. Ultimately the goal is the compare the calculated energies so that we can find the equilibrium lattice parameter, bulk modulus, and its pressure derivative.

def birch_murnaghan(V, E0, B0, BP, V0):

    """3rd-order Birch-Murnaghan EOS."""
    eta = (V0 / V) ** (2.0 / 3.0)
    term1 = 9.0 * V0 * B0 / 16.0
    term2 = (eta - 1.0) ** 3 * BP
    term3 = (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    return E0 + term1 * (term2 + term3)

def main():

    output_files = sorted(glob.glob("calculations/bulk_eos_final_corrected_k14/scf_a_*.out"))

    volumes, energies = [], []
    print(f"Found {len(output_files)} QE output files. Parsing...\n")

    for fp in output_files:
        a_val, energy = parse_qe_output(fp)
        if a_val is not None and energy is not None:
            vol = (a_val**3) / 2.0   # BCC single-atom primitive cell volume formula
            volumes.append(vol)
            energies.append(energy)
            print(f"Parsed: {os.path.basename(fp)} | a = {a_val:.4f} Å | E = {energy:.6f} Ry")

    if len(energies) < 4:
        print("ERROR: Less than 4 valid data points parsed.")
        return

    vols, ens = np.array(volumes), np.array(energies)

    # Here I'm passing all the (Volume, Energy) pairs to the Birch-Murnaghan EOS function
    popt, _ = curve_fit(
        birch_murnaghan,
        vols,
        ens,
        p0=[np.min(ens), 0.01, 4.0, vols[np.argmin(ens)]],
    )

    # Now trying to extract properties (of things like a, V0, B0) from fitted parameters
    E0_fit, B0_fit, BP_fit, V0_fit = popt

    # This part basically reconstructs bulk lattic parameter
    a_eq = (2.0 * V0_fit) ** (1.0 / 3.0)
    # B0_GPa = B0_fit * 14710.5 This is where the mistake was, I was using the wrong conversion factor for B0 (Ry/Å³ → GPa)
    B0_GPa = B0_fit * 2179.87

    print("\n" + "=" * 50)
    print("BIRCH-MURNAGHAN EOS FIT RESULTS")
    print("=" * 50)
    print(f"Equilibrium Lattice Constant (a_eq): {a_eq:.4f} Å")
    print(f"Equilibrium Cell Volume (V0): {V0_fit:.4f} Å³")
    print(f"Bulk Modulus (B0): {B0_GPa:.2f} GPa")
    print(f"Equilibrium Energy (E0): {E0_fit:.8f} Ry")
    print(f"Pressure Derivative (B0'): {BP_fit:.3f}")
    print("=" * 50)

    with open("eos_results.txt", "w", encoding="utf-8") as f:
        f.write("=" * 50 + "\n")
        f.write("BIRCH-MURNAGHAN EOS FIT RESULTS\n")
        f.write("=" * 50 + "\n")
        f.write(f"Equilibrium Lattice Constant (a_eq): {a_eq:.4f} Angstroms\n")
        f.write(f"Equilibrium Cell Volume (V0): {V0_fit:.4f} Angstroms^3\n")
        f.write(f"Bulk Modulus (B0): {B0_GPa:.2f} GPa\n")
        f.write(f"Equilibrium Energy (E0): {E0_fit:.8f} Ry\n")
        f.write(f"Pressure Derivative (B0'): {BP_fit:.3f}\n")
        f.write("=" * 50 + "\n")

# So the above gives us a_eq, V0, E0, B0, and B0' from the Birch-Murnaghan EOS fit.

if __name__ == "__main__":
    main()

# Updates:
# After running I noticed that calculations were alright, but Bulk Modulus gave a value of 1175.13 GPa, too high to be reasonable.
# Now I'm going to run a K-Points convergence test because insufficient K-Point mesh could distort the EOS curvature even if the energy curve appears smooth