# Apparently there were some problems with my bulk modulus calculation and I suspected that the K-Point mesh was insufficient, so Ima run a K-Points convergence test real quick
# This file follows the Fit EOS.py file

import os

a = 2.8331

ecutwfc = 65.0
ecutrho = 782.0

k_points_tests = [
    # (8, 8, 8),
    # (10, 10, 10),
    # (12, 12, 12),
    (14, 14, 14),
    (16, 16, 16),
]

output_dir = os.path.join("calculations", "bulk_kpoints_convergence")
os.makedirs(output_dir, exist_ok=True)

def write_qe_kpoint_scf(a, filename, ecutwfc, ecutrho, k_grid):

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

if __name__ == "__main__":

    for k_grid in k_points_tests:
    
        k = k_grid[0]
        filename = f"scf_k_{k}x{k}x{k}.in"
        filepath = os.path.join(output_dir, filename)

        write_qe_kpoint_scf(
            a=a,
            filename=filepath,
            ecutwfc=ecutwfc,
            ecutrho=ecutrho,
            k_grid=k_grid
        )

        print(f"{filename} created")