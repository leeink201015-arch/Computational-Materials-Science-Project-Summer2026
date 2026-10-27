import os

os.makedirs("out_2x2x3", exist_ok=True)

original_positions = [
    (0.0, 0.0, 0.4398482736, 1, 1, 1),
    (0.5, 0.5, 0.5000000000, 0, 0, 0),
    (0.0, 0.0, 0.5601517264, 1, 1, 1),
]

positions_2x2x3 = []

for x, y, z, fx, fy, fz in original_positions:
  for i in range(2):
    for j in range(2):
      positions_2x2x3.append(f"Fe {(x+i)/2:.10f} {(y+j)/2:.10f} " f"{z:.10f} {fx} {fy} {fz}")

atomic_positions = "\n".join(positions_2x2x3)
assert len(positions_2x2x3) == 12

slab_input = f"""&CONTROL
  calculation    = 'relax'
  restart_mode   = 'from_scratch'
  prefix         = 'fe100_slab_2x2x3'
  pseudo_dir     = '/usr/share/espresso/pseudo/'
  outdir         = './out_2x2x3/'
  forc_conv_thr = 1.0D-3
/
&SYSTEM
  ibrav       = 6
  celldm(1)   = 10.70756
  celldm(3)   = 4.1479
  nat         = 12
  ntyp        = 1
  ecutwfc     = 65.0
  ecutrho     = 782.0
  occupations = 'smearing'
  smearing    = 'marzari-vanderbilt'
  degauss     = 0.02
  nspin       = 2
  starting_magnetization(1) = 0.5
/
&ELECTRONS
  conv_thr    = 1.0D-8
  mixing_beta = 0.3
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  Fe  55.845  Fe.pbe-spn-rrkjus_psl.0.2.1.UPF

ATOMIC_POSITIONS (crystal)
{atomic_positions}

K_POINTS automatic
  2 2 1 0 0 0
"""

# List of things changed from 2x2x7 to 2x2x3:
# 1) From 2x2x7 = 28 to 2x2x3 = 12
# 2) Remove top2 and bottom2 layers, only keep middle 3 layers
# 3) Change calculation from SCF to relaxed. The original 2x2x7 model used coordinates from the initially relaxed OG 1x1x7 model, now since we changed the number of layers, the remaining 3 layers change physical environment
# 4) So now I needed to add the IONS block as well
# 5) K-POINTS from 5 5 1 0 0 0 to 2 2 1 0 0 0 to reduce computation time
 
filename = "fe100_slab_2x2x3.in"

with open(filename, "w") as f:
  f.write(slab_input)

print(f"Created {filename}")