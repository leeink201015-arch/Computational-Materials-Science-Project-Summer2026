# Now scaling the lattice parameter of the BCC Fe structure by a factor of a certain range
# Then generate QE input files for each volume of the unit cell

import os
import numpy as np
from pymatgen.core import Structure, Lattice

# This is the base (initial) BCC Iron Structure (setting a = 2.866 Å)
a_initial = 2.866

# Parameters for EOS Scaling and QE Calculation
scale_factors = np.linspace(0.95, 1.05, 11)  # Created 11 evenly spaced numbers between 0.95 and 1.05, picked 0.95 and 1.05 to cover a range of ±5% around the initial lattice parameter
ecutwfc = 65.0
ecutrho = 782.0
k_grid = (8, 8, 8)

# Create output directory...
output_dir = os.path.join("calculations", "bulk_eos_final_corrected")  # Ik the calculcation won't be 100% "final corrected" but this name helps me identify the file names for the calculations because I ran QE before that already and got pre off results that are completely wrong so yea
os.makedirs(output_dir, exist_ok=True)

def write_qe_bulk_scf(a, filename, ecutwfc, ecutrho, k_grid):

    # a = structure.lattice.a

    with open(filename, "w") as file:
        file.write("&CONTROL\n")
        file.write("calculation = 'scf'\n")
        file.write("restart_mode = 'from_scratch'\n")
        file.write("prefix = 'fe_bulk'\n")
        file.write("pseudo_dir = '/usr/share/espresso/pseudo/'\n")
        file.write("outdir = './tmp/'\n")
        file.write("/\n\n")

        file.write("&SYSTEM\n")
        file.write("ibrav = 3\n")
        file.write(f"celldm(1) = {a / 0.52917721092:.6f}  ! converted Å to Bohr\n")
        file.write("nat = 1\n")
        file.write("ntyp = 1\n")
        file.write(f"ecutwfc = {ecutwfc}\n")
        file.write(f"ecutrho = {ecutrho}\n")
        file.write("occupations = 'smearing'\n")
        file.write("smearing = 'mv'\n")
        file.write("degauss = 0.02\n")
        file.write("nspin = 2\n")
        file.write("starting_magnetization(1) = 0.5\n")
        file.write("/\n\n")

        file.write("&ELECTRONS\n")
        file.write("conv_thr = 1.0d-8\n")
        file.write("mixing_beta = 0.7\n")
        file.write("/\n\n")

        file.write("ATOMIC_SPECIES\n")
        file.write("Fe 55.845 Fe.pbe-spn-rrkjus_psl.0.2.1.UPF\n")
        file.write("\n")
        file.write("ATOMIC_POSITIONS (crystal)\n")
        file.write("Fe 0.0 0.0 0.0\n")
        file.write("\n\n")

        file.write("K_POINTS (automatic)\n")
        file.write(f"{k_grid[0]} {k_grid[1]} {k_grid[2]} 0 0 0\n")

# Loop over scaling factors and write inputs
if __name__ == "__main__":

    print(f"Generating bulk EOS input files in '{output_dir}'...")
    
    for scale in scale_factors:
        a_scaled = a_initial * scale

        file_name = f"scf_a_{a_scaled:.4f}.in"
        file_path = os.path.join(output_dir, file_name)

        write_qe_bulk_scf(
            a=a_scaled,
            filename=file_path,
            ecutwfc=ecutwfc,
            ecutrho=ecutrho,
            k_grid=k_grid
        )
        print(f" Created: {file_name} (a = {a_scaled:.4f} Å)")