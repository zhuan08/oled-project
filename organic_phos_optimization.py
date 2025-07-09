import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

# Create an output folder for the geometries, if it doesn't already exist
geom_dir_name = 'organic_phos_geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass
# -------- Read csv files for smiles strings --------
data = pd.read_csv('smiles_energy.csv', sep=',')
data = data[["mol_id", "smiles", "absorption_energy", "emission_energy"]]
# -------- list of smiles strings --------
smiles = []
for i in data["smiles"]:
    smiles.append(i)
# -------- list of mol_ids --------
mol_ids = []
for i in data["mol_id"]:
    mol_ids.append(i)
# -------- smiles strings to rdkit molecules --------
diff_energy_list = []
mol_error_pair = []

for mol_id, smile in zip(mol_ids, smiles):
    # Path name, based on the molecule identifier
    geom_path = os.path.join(geom_dir_name, f'{mol_id}.xyz')
    # First, check if there's already a geometry saved, and if so, just load it
    if os.path.exists(geom_path):
        atom = ase.io.read(filename=geom_path)
        print(f'Found geometry for file: {mol_id}')
    else:
        print(f'Did not find geometry for file: {mol_id}')
        molecule = Chem.MolFromSmiles(smile, removeHs=False, sanitize=False)
        if molecule is None:
            raise ValueError(f'MolFromMol2File returned None on {mol_id}.mol2')
        molecule = Chem.AddHs(molecule)
        # Embed molecule
        try:
            rdDistGeom.EmbedMolecule(molecule)
            print(f'Embeded molecule on {mol_id}')
        except Exception as e:
            print(f'Did no embeded molecule on {mol_id}')
            mol_error_pair.append((mol_id, e))
            continue
        # conform molecule
        try:
            conf_mol = molecule.GetConformer()
            print(f'Conformed molecule on: {mol_id}')
        except Exception as e:
            print(f'Did not conform molecule on: {mol_id}')
            mol_error_pair.append((mol_id, e))
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
            mol_error_pair.append((mol_id, e))
            continue
            # Write the geometry to a file
        ase.io.write(filename=geom_path, images=atom)
        mol_error_pair.append((mol_id, "No Error"))
i = 0
for i in range(len(mol_error_pair)):
    print(mol_error_pair[i])