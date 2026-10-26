import os
import numpy as np

VSEPR_VECTORS = {
    "linear_diatomic": np.array([
        [1.0, 0.0, 0.0]
    ])
}

reacting_molecules = {
    "CO": {
        "centre": "C",
        "ligands": ["O"],
        "geometry": "linear_diatomic",
        "bond_length": [1.128]
    }
}

involved_elements_informations = {
    "C": {"mass": 12.011, "pseudopotential": "C.pbe-n-kjpaw_psl.0.1.UPF"},
    "O": {"mass": 15.999, "pseudopotential": "O.pbe-kjpaw.UPF"}
}

def molecule_generated(centre, ligands, geometry, bond_length):

    molecular_geometry = VSEPR_VECTORS[geometry]

    symbols = [centre]
    coordinates = [[0.0, 0.0, 0.0]]

    for element, vector, length in zip(
        ligands, molecular_geometry, bond_length
    ):
        symbols.append(element)
        coordinates.append((vector * length).tolist())

    return symbols, np.array(coordinates)

CO_data = reacting_molecules["CO"]

symbols, coordinates = molecule_generated(
    centre=CO_data["centre"],
    ligands=CO_data["ligands"],
    geometry=CO_data["geometry"],
    bond_length=CO_data["bond_length"]
)

print(symbols)
print(coordinates)

box_size = 15.0  # So that the CO molecules are spread out from each other and the periodic boundary conditions that QE uses don't cause interactions between CO molecules.

molecule_center = np.mean(coordinates, axis=0)

# Second thing is that the CO molecule has to be centered in the box

box_center = np.array([
    box_size / 2,
    box_size / 2,
    box_size / 2
])

coordinates_in_box = coordinates - molecule_center + box_center

# Here just defining the simulation cell

print(coordinates_in_box)

cell = np.eye(3) * box_size

print(cell)

def write_qe_input(symbols, coordinates, cell, output_file):

    with open(output_file, "w") as file:

        file.write("&CONTROL\n")
        file.write("calculation = 'relax'\n")
        file.write("prefix = 'CO_gas'\n")
        file.write("pseudo_dir = '/usr/share/espresso/pseudo/'\n")
        file.write("outdir = './tmp/'\n")
        file.write("/\n\n")

        file.write("&SYSTEM\n")
        file.write(f"ibrav = 0\n")
        file.write("nat = 2\n")
        file.write("ntyp = 2\n")
        file.write("ecutwfc = 65.0\n")
        file.write("ecutrho = 782.0\n")
        file.write("occupations = 'fixed'\n")
        file.write("nspin = 1\n")
        file.write("/\n\n")

        file.write("&ELECTRONS\n")
        file.write("conv_thr = 1.0d-8\n")
        file.write("mixing_beta = 0.7\n")
        file.write("/\n\n")

        file.write("&IONS\n")
        file.write("ion_dynamics = 'bfgs'\n")
        file.write("/\n\n")

        file.write("ATOMIC_SPECIES\n")
        for element in dict.fromkeys(symbols):
            data = involved_elements_informations[element]
            file.write(f"{element} {data['mass']} {data['pseudopotential']}\n")
        file.write("\n")

        file.write("CELL_PARAMETERS angstrom\n")  # First time I forgot to add "angstrom" so the units were wrong
        for vector in cell:
            file.write(f"{vector[0]:12.6f} {vector[1]:12.6f} {vector[2]:12.6f}\n")
        file.write("\n")

        file.write("ATOMIC_POSITIONS angstrom\n")  # Same thing as line 109, forgot to include the units now lemme just regenerate the input file before running QE again
        for element, position in zip(symbols, coordinates):
            file.write(f"{element} {position[0]:12.6f} {position[1]:12.6f} {position[2]:12.6f}\n")
        file.write("\n")

        file.write("K_POINTS gamma\n")  # Write an instruction to QE to only use the gamma point for the calculation

write_qe_input(symbols, coordinates_in_box, cell, "CO_gas.in")
print("QE input file 'CO_gas.in' generated successfully")