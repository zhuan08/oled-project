import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from ase.calculators.gamess_us import GAMESSUS
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

tblite_calc_singlet = TBLite(multiplicity=1)
tblite_calc_triplet = TBLite(multiplicity=3)
gamess_calc_singlet = GAMESSUS(contrl={'mult': 1}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
gamess_calc_triplet = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
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
    else:
        molecule = Chem.MolFromSmiles(smile, removeHs=False, sanitize=False)
        if molecule is None: 
            raise ValueError(f'MolFromSmiles returned None on {mol_id}')
        molecule = Chem.AddHs(molecule)
        # Embed molecule
        try:
            rdDistGeom.EmbedMolecule(molecule)
        except Exception as e:
            mol_error_pair.append((mol_id, e))
            continue
        # conform molecule
        try:
            conf_mol = molecule.GetConformer()
        except Exception as e:
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
        atom.calc = tblite_calc_triplet
        opt = BFGS(atom, logfile=None, trajectory=None)
        try:
            opt.run(fmax=0.05)
        except Exception as e:
            mol_error_pair.append((mol_id, e))
            continue
            # Write the geometry to a file
        ase.io.write(filename=geom_path, images=atom)
        mol_error_pair.append((mol_id, "No Error"))
i = 0
for i in range(len(mol_error_pair)):
    print(mol_error_pair[i])