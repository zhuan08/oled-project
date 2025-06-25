import os
import ase.io
from octahedral_embed import octahedral_embed
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Draw
import pandas as pd
import numpy as np
IPythonConsole.ipython_3d = True

# Create an output folder for the geometries, if it doesn't already exist
geom_dir_name = 'geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass

# -------- Read csv files for smiles strings --------
data = pd.read_csv('oled_smiles_energy.csv', sep=',')
data = data[["mol_id", "smiles", "absorption_energy", "emission_energy"]]

# -------- list of smiles strings --------
smiles = []
for i in data["smiles"]:
    smiles.append(i)
mol_id = []
for i in data["mol_id"]:
    mol_id.append(i)
error_mol = []
for i in data["mol_id"]:
    error_mol.append(i)
# -------- smiles to rdkit molecule --------
diff_energy_list = []
final_mol_id = []
error_msg = []

for mol_id, smile in zip(mol_id, smiles):
    # Draw.MolToImage(molecule).show()
    # Geometry optimization
    # Path name, based on the molecule identifier
    geom_path = os.path.join(geom_dir_name, f'{mol_id}.xyz')
    # First, check if there's already a geometry saved, and if so, just load it
    if os.path.exists(geom_path):
        atom = ase.io.read(filename=geom_path)
        print(f'Found geometry for file: {mol_id}')
        final_mol_id.append(mol_id)
    else:
        print(f'Did not find geometry for file: {mol_id}')
        # Geometry optimization
        molecule = Chem.MolFromSmiles(smile)
        molecule = Chem.AddHs(molecule)
            # -------- rdkit molecule positions --------
        try:
            octahedral_embed(molecule, isomer="tridentate")
            print(f'Used octahedral_embed on: {mol_id}')
        except Exception as e: 
            print(f'Did not use octahedral_embed on: {mol_id}')
            error_msg.append(e)
            continue
        try:
            conf_mol = molecule.GetConformer()
            print(f'Conformed molecule on: {mol_id}')
        except Exception as e:
            print(f'Did not conform molecule on: {mol_id}')
            error_msg.append(e)
            continue
        pos_mol = conf_mol.GetPositions()
        # -------- rdkit.atom obj to get symbols --------
        atom_mol = molecule.GetAtoms()
        atom_sym = []
        for char in atom_mol:
            atom_sym.append(char.GetSymbol())
        print(atom_sym)
        atom = Atoms(atom_sym, positions=pos_mol)  # atom object to calculate triplet - singlet gap (in triplet geometry)
        atom.calc = TBLite(multiplicity=3)
        opt = BFGS(atom, logfile=None, trajectory=None)
        try:
            opt.run(fmax=0.05)
            print(f'Optimized run on: {mol_id}')
        except Exception as e:
            print(f'Did not optimize run on: {mol_id}')
            error_msg.append(e)
            continue
        # Write the geometry to a file
        ase.io.write(filename=geom_path, images=atom)
    diff_energy = 0
    for mult_e in [1, 3]:
        atom.calc = TBLite(multiplicity=mult_e)
        if mult_e == 1:
            energy = atom.get_potential_energy()
            diff_energy -= energy
            print(f'{mol_id}', ' molecule singlet energy (triplet geometry)', ': %5.2f eV' % energy)
        if mult_e == 3:
            energy = atom.get_potential_energy()
            diff_energy += energy
            print(f'{mol_id}', ' molecule triplet energy (triplet geometry)', ': %5.2f eV' % energy)
    diff_energy_list.append(diff_energy)
    error_msg.append('No Error Message')

i = 0
for i in range(len(final_mol_id)):
    print(final_mol_id[i],',',diff_energy_list[i],',',error_msg[i])
