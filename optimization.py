import os
import glob
import ase.io
from octahedral_embed import octahedral_embed
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Draw
IPythonConsole.ipython_3d = True

geom_dir_name = 'geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass

mol_id = []
path_id = []
final_mol_id = []
error_msg = []
diff_energy_list = []

# Creating list of mol_ids
for mol in glob.glob('structures_mol2/*.mol2'):
    if '.mol2' in mol:
        updated_file = mol.replace('structures_mol2', '')
        updated_file = updated_file.replace('.mol2', '')
        updated_file = updated_file.replace('/', '')
        mol_id.append(updated_file)
print(mol_id)

# Creating list of paths for mol2 molecules
for path in glob.glob('structures_mol2/*.mol2'):
    if '.mol2' in path:
        path_id.append(path)
print(path_id)

for mol_id, path_id in zip(mol_id, path_id):
    geom_path = os.path.join(geom_dir_name, f'{mol_id}.xyz')
    if os.path.exists(geom_path):
        atom = ase.io.read(filename=geom_path)
        final_mol_id.append(mol_id)
    else:
        molecule = Chem.MolFromMol2File(path_id, removeHs=False, sanitize=False)
        if molecule is None:
            raise ValueError(f'MolFromMol2File returned None on {mol_id}.mol2')
        molecule = Chem.AddHs(molecule)
        try:
            octahedral_embed(molecule, isomer="tridentate")
        except Exception as e: 
            error_msg.append(e)
            continue
        try:
            conf_mol = molecule.GetConformer()
        except Exception as e:
            error_msg.append(e)
            continue
        pos_mol = conf_mol.GetPositions()
        atom_mol = molecule.GetAtoms()
        atom_sym = []
        for char in atom_mol:
            atom_sym.append(char.GetSymbol())
        print(atom_sym)
        atom = Atoms(atom_sym, positions=pos_mol)
        atom.calc = TBLite(multiplicity=3)
        opt = BFGS(atom, logfile=None, trajectory=None)
        try:
            opt.run(fmax=0.05)
        except Exception as e:
            error_msg.append(e)
            continue
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