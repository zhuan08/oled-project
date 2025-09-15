import os
import glob
import ase.io
from octahedral_embed import octahedral_embed
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from ase.calculators.gamess_us import GAMESSUS
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Draw
IPythonConsole.ipython_3d = True

tblite_calc_singlet = TBLite(multiplicity=1)
tblite_calc_triplet = TBLite(multiplicity=3)
gamess_calc_singlet = GAMESSUS(contrl={'mult': 1}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
gamess_calc_triplet = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')

geom_dir_name = 'oled_geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass

mol_error_pair = []
mol_ids = []
path_id = []
diff_energy_list = []

# Creating list of mol_ids
for mol in glob.glob('structures_mol2/*.mol2'):
    if '.mol2' in mol:
        updated_file = mol.replace('structures_mol2', '')
        updated_file = updated_file.replace('.mol2', '')
        updated_file = updated_file.replace('/', '')
        mol_ids.append(updated_file)

# Creating list of paths for mol2 molecules
for path in glob.glob('structures_mol2/*.mol2'):
    if '.mol2' in path:
        path_id.append(path)

for mol_id, path_id in zip(mol_ids, path_id):
    geom_path = os.path.join(geom_dir_name, f'{mol_id}.xyz')
    if os.path.exists(geom_path):
        atom = ase.io.read(filename=geom_path)
    else:
        molecule = Chem.MolFromMol2File(path_id, removeHs=False, sanitize=False)
        if molecule is None:
            raise ValueError(f'MolFromMol2File returned None on {mol_id}.mol2')
        try:
            rdDistGeom.EmbedMolecule(molecule)
        except Exception as e:
            mol_error_pair.append((mol_id, e))
            continue
        try:
            conf_mol = molecule.GetConformer()
        except Exception as e:
            mol_error_pair.append((mol_id, e))
            continue
        pos_mol = conf_mol.GetPositions()
        atom_mol = molecule.GetAtoms()
        atom_sym = []
        for char in atom_mol:
            atom_sym.append(char.GetSymbol())
        atom = Atoms(atom_sym, positions=pos_mol)
        atom.calc = tblite_calc_triplet
        opt = BFGS(atom, logfile=None, trajectory=None)
        try:
            opt.run(fmax=0.05)
        except Exception as e:
            mol_error_pair.append(mol_id, e)
            continue
        ase.io.write(filename=geom_path, images=atom)
        mol_error_pair.append((mol_id, "No Error"))
i = 0
for i in range(len(mol_error_pair)):
    print(mol_error_pair[i])