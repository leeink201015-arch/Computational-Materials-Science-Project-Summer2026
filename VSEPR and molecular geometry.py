import numpy as np

reacting_molecules = {
    "CO": {"centre": "C", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.128]},
    "NO": {"centre": "N", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.150]},

    "H2": {"centre": "H", "ligands": ["H"], "geometry": "linear_diatomic", "bond_length": [0.741]},
    "O2": {"centre": "O", "ligands": ["O"], "geometry": "linear_diatomic", "bond_length": [1.208]},
    "N2": {"centre": "N", "ligands": ["N"], "geometry": "linear_diatomic", "bond_length": [1.098]},
    
    "CO2": {"centre": "C", "ligands": ["O", "O"], "geometry": "linear_triatomic", "bond_length": [1.162, 1.162]},
    "N2O": {"centre": "N", "ligands": ["N", "O"], "geometry": "linear_triatomic", "bond_length": [1.120, 1.190]},

    "NO2": {"centre": "N", "ligands": ["O", "O"], "geometry": "bent", "bond_length": 1.197, "bond_angle": [134.1, 134.1]},
    "SO2": {"centre": "S", "ligands": ["O", "O"], "geometry": "bent", "bond_length": 1.430, "bond_angle": [119.5, 119.5]},
    "H2O": {"centre": "O", "ligands": ["H", "H"], "geometry": "bent", "bond_length": 0.958, "bond_angle": [104.5, 104.5]},
    
    "CH4": {"centre": "C", "ligands": ["H", "H", "H", "H"], "geometry": "tetrahedral", "bond_length": [1.087, 1.087, 1.087, 1.087]},
    "NH3": {"centre": "N", "ligands": ["H", "H", "H"], "geometry": "trigonal_pyramidal", "bond_length": [1.012, 1.012, 1.012], "lp_centre_ligand_angle": 112.0}
}

VSEPR_VECTORS = {
    "linear_diatomic": np.array([
        [1.0, 0.0, 0.0]
        #[-1.0, 0.0, 0.0]
    ]),

    "linear_triatomic": np.array([
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0]
    ]),

    "trigonal_planar": np.array([
        [ 1.0, 0.0, 0.0],
        [-0.5, np.sqrt(3)/2, 0.0],
        [-0.5, -np.sqrt(3)/2, 0.0]
    ]),
    
    "tetrahedral": np.array([
        [0.0, 0.0, 1.0],
        [2*np.sqrt(2)/3, 0.0, -1/3],
        [-np.sqrt(2)/3,  np.sqrt(2/3), -1/3],
        [-np.sqrt(2)/3, -np.sqrt(2/3), -1/3]
    ]),

    "square_planar": np.array([
        [np.sqrt(2)/2, np.sqrt(2)/2, 0.0],
        [-np.sqrt(2)/2, np.sqrt(2)/2, 0.0],
        [np.sqrt(2)/2, -np.sqrt(2)/2, 0.0],
        [-np.sqrt(2)/2, -np.sqrt(2)/2, 0.0],
    ]),

    "t_shaped": np.array([
        [0.0, 1.0, 0.0],
        [np.cos(np.radians(86.2)), np.sin(np.radians(86.2)), 0.0],
        [-np.cos(np.radians(86.2)), np.sin(np.radians(86.2)), 0.0]
    ]),

    "trigonal_bipyramidal": np.array([
        [0.0, 0.0, 1.05],
        [0.0, 0.0, -1.05],
        [1.0, 0.0, 0.0],
        [np.cos(np.radians(120)), np.sin(np.radians(120)), 0.0],
        [np.cos(np.radians(240)), np.sin(np.radians(240)), 0.0]
    ]),

    "octahedral": np.array([
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0]
    ]),
}

def bent_vectors (bond_angle):
    half_bond_angle = np.radians(bond_angle/2)

    return np.array([
        [np.sin(half_bond_angle), np.cos(half_bond_angle), 0.0],
        [-np.sin(half_bond_angle), np.cos(half_bond_angle), 0.0],
    ])

def trigonal_pyramidal_vectors (lp_centre_ligand_angle):
    supplementary_angle = np.radians(180)-np.radians(lp_centre_ligand_angle)

    return np.array([
        [1.0 * np.sin(supplementary_angle), 0.0, np.cos(supplementary_angle)],
        [-0.5 * np.sin(supplementary_angle), np.sqrt(3)/2 * np.sin(supplementary_angle), np.cos(supplementary_angle)],
        [-0.5 * np.sin(supplementary_angle), -np.sqrt(3)/2 * np.sin(supplementary_angle), np.cos(supplementary_angle)]
    ])

def molecule_generated (centre, ligands, geometry, bond_length, bond_angle=None, lp_centre_ligand_angle=None):

    if geometry == "bent":
        molecular_geometry = bent_vectors(bond_angle)

    elif geometry == "trigonal_pyramidal":
        molecular_geometry = trigonal_pyramidal_vectors(lp_centre_ligand_angle)

    else:
        molecular_geometry = VSEPR_VECTORS[geometry]

    symbols = [centre]
    coordinates = [[0.0, 0.0, 0.0]]

    for element, vector in zip(ligands, molecular_geometry):
        symbols.append(element)
        coordinates.append((vector * bond_length).tolist())
        
    return symbols, np.array(coordinates)