# Oct 8: So now I looked at my Fe(100) slab and it was a 1x1 slab. Now I need a 2x2 slab, so copy pasted the code from that file and going to do some modifications to get ready for calculating E_ads of CO on Fe(100) surface.
# Note that I originally had a 1x1 because I was just characterizing the clean Fe(100) surface, NOT INVESTIGATING CO ADSORTION
# Now I'm investigating the ads and introducing CO changed the requirement of the model

# QE uses periodic boundary conditions and because I don't want the CO molecules to be "densely packed" in the simulation cell, I need to make the slab larger

import os

os.makedirs("out_2x2", exist_ok=True)

original_positions = [
    (0.0, 0.0, 0.3189370858, 1, 1, 1),
    (0.5, 0.5, 0.3778786719, 1, 1, 1),
    (0.0, 0.0, 0.4398482736, 1, 1, 1),
    (0.5, 0.5, 0.5000000000, 0, 0, 0),
    (0.0, 0.0, 0.5601517264, 1, 1, 1),
    (0.5, 0.5, 0.6221213281, 1, 1, 1),
    (0.0, 0.0, 0.6810629142, 1, 1, 1),
]

positions_2x2 = []

for x, y, z, fx, fy, fz in original_positions:
  for i in range(2):
    for j in range(2):
      positions_2x2.append(f"Fe {(x+i)/2:.10f} {(y+j)/2:.10f} " f"{z:.10f} {fx} {fy} {fz}")

atomic_positions = "\n".join(positions_2x2)
assert len(positions_2x2) == 28

slab_input = f"""&CONTROL
  calculation  = 'scf'
  restart_mode = 'from_scratch'
  prefix       = 'fe100_slab_2x2'
  pseudo_dir   = '/usr/share/espresso/pseudo/'
  outdir       = './out_2x2/'
/
&SYSTEM
  ibrav       = 6
  celldm(1)   = 10.70756
  celldm(3)   = 4.1479
  nat         = 28
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
ATOMIC_SPECIES
  Fe  55.845  Fe.pbe-spn-rrkjus_psl.0.2.1.UPF

ATOMIC_POSITIONS (crystal)
{atomic_positions}

K_POINTS automatic
  5 5 1 0 0 0  # Results obtained from convergence test
"""

# List of things changed so far:
# 1) prefix changed to fe100_slab_2x2
# 2) outdir changed to ./out_2x2/
# 3) nat changed to 28 because there are now 28 Fe atoms in the 2x2 slab
# 4) K_POINTS changed to 5 5 1 because the slab is now larger, so fewer k-points are needed to sample the Brillouin zone

filename = "fe100_slab_2x2.in"

with open(filename, "w") as f:
  f.write(slab_input)

print(f"Created {filename}")