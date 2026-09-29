import os

os.makedirs("out", exist_ok=True)

# k_point_tests = [
#     (8, 8, 1),
#     (10, 10, 1),
#     (12, 12, 1),
#     (14, 14, 1),
# ]

# Writing the converged fe100_slab.in file, note that:
# 1) 7-layer symmetric Fe(100) clean-surface slab
# 2) a_eq = 2.8326 Å from the final bulk EOS calculation
# 3) ~15 Å total vacuum between periodic slab images, let this be the initial choice as it is quite reasonable to pick 15 Å (because for metallic surface calculations, roughly 10–20 Å is a practical starting range)

# 7 Fe layers contain 6 interlayer spacings of a/2, so slab thickness = 3a = 8.4978 Å, adding 15 Å vacuum gives c = 23.4978 Å.
# I set celldm(1) = 5.35378 because celldm(1) = a_eq / Bohr_radius = 2.8326 Å / 0.529177 Å ≈ 5.35378
# I set celldm(3) = 8.2958 because c/a = 23.4978 / 2.8326 ≈ 8.2958

slab_input = """&CONTROL
  calculation  = 'relax'
  restart_mode = 'from_scratch'
  prefix       = 'fe100_slab'
  pseudo_dir   = '/usr/share/espresso/pseudo/'
  outdir       = './out/'
  forc_conv_thr= 1.0D-3
/
&SYSTEM
  ibrav       = 6
  celldm(1)   = 5.35378
  celldm(3)   = 8.2958
  nat         = 7
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
  Fe  0.0000  0.0000  0.319175  1 1 1
  Fe  0.5000  0.5000  0.379450  1 1 1
  Fe  0.0000  0.0000  0.439725  1 1 1
  Fe  0.5000  0.5000  0.500000  0 0 0
  Fe  0.0000  0.0000  0.560275  1 1 1
  Fe  0.5000  0.5000  0.620550  1 1 1
  Fe  0.0000  0.0000  0.680825  1 1 1

K_POINTS automatic
  10 10 1 0 0 0  # Results obtained from convergence test
"""

filename = "fe100_relax.in"

with open(filename, "w") as f:
  f.write(slab_input)

print(f"Created {filename}")