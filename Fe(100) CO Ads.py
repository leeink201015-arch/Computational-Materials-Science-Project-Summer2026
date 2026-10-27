from pathlib import Path

# The model:
# CO on the top surface of relaxed 2x2x3 Fe(100) slab
# C-down, O-up, CO molecule initially perpendicular to the surface

# Lattice parameters from the completed slab calculation
celldm1 = 10.70756
celldm3 = 4.1479

bohr_to_angstrom = 0.529177210903
cell_height = celldm1 * celldm3 * bohr_to_angstrom

# Final relaxed coordinates (from October 13)
fe_positions = [
    (0.0000000000, 0.0000000000, 0.4401518937),
    (0.0000000000, 0.5000000000, 0.4401600566),
    (0.5000000000, 0.0000000000, 0.4401600566),
    (0.5000000000, 0.5000000000, 0.4401514763),

    (0.2500000000, 0.2500000000, 0.5000000000),
    (0.2500000000, 0.7500000000, 0.5000000000),
    (0.7500000000, 0.2500000000, 0.5000000000),
    (0.7500000000, 0.7500000000, 0.5000000000),

    (0.0000000000, 0.0000000000, 0.5598481063),
    (0.0000000000, 0.5000000000, 0.5598399434),
    (0.5000000000, 0.0000000000, 0.5598399434),
    (0.5000000000, 0.5000000000, 0.5598485237),
]

# Selecting the first top-layer Fe atom
x_fe, y_fe, z_fe = fe_positions[8]

# Initial adsorption geometry (angstroms)
fe_c_distance = 2.0
co_bond_length = 1.1411

# Converting vertical distances to fractional coordinates
z_c = z_fe + fe_c_distance / cell_height
z_o = z_c + co_bond_length / cell_height

assert len(fe_positions) == 12
assert z_o < 1.0, "CO extends beyond the cell boundary"

# Building atomic positions and relaxation flags here
atomic_lines = []

for index, (x, y, z) in enumerate(fe_positions):
    flags = "0 0 0" if 4 <= index <= 7 else "1 1 1"
    atomic_lines.append(
        f"Fe {x:.10f} {y:.10f} {z:.10f} {flags}"
    )

atomic_lines.append(
    f"C {x_fe:.10f} {y_fe:.10f} {z_c:.10f} 1 1 1"
)
atomic_lines.append(
    f"O {x_fe:.10f} {y_fe:.10f} {z_o:.10f} 1 1 1"
)

atomic_positions = "\n".join(atomic_lines)

qe_input = f"""&CONTROL
  calculation = 'relax'
  restart_mode = 'from_scratch'
  prefix = 'fe100_co_ontop'
  pseudo_dir = '/usr/share/espresso/pseudo/'
  outdir = './out_co_ontop/'
  forc_conv_thr = 1.0D-3
/
&SYSTEM
  ibrav = 6
  celldm(1) = {celldm1}
  celldm(3) = {celldm3}
  nat = 14
  ntyp = 3

  ecutwfc = 65.0
  ecutrho = 782.0

  occupations = 'smearing'
  smearing = 'marzari-vanderbilt'
  degauss = 0.02

  nspin = 2
  starting_magnetization(1) = 0.5
/
&ELECTRONS
  conv_thr = 1.0D-8
  mixing_beta = 0.3
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
Fe 55.845 Fe.pbe-spn-rrkjus_psl.0.2.1.UPF
C  12.011 C.pbe-n-kjpaw_psl.0.1.UPF
O  15.999 O.pbe-kjpaw.UPF

ATOMIC_POSITIONS (crystal)
{atomic_positions}

K_POINTS automatic
2 2 1 0 0 0
"""

filename = Path("fe100_2x2x3_co_ontop.in")
filename.write_text(qe_input)

print(f"Created: {filename}")
print(f"Number of atoms: {len(atomic_lines)}")
print(f"Cell height: {cell_height:.4f} Å")
print(f"Initial Fe-C distance: {fe_c_distance:.4f} Å")
print(f"Initial C-O distance: {co_bond_length:.4f} Å")
print(f"Carbon fractional z: {z_c:.10f}")
print(f"Oxygen fractional z: {z_o:.10f}")