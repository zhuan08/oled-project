from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

# -------- Read csv files for smiles strings --------
data = pd.read_csv('smiles_energy.csv', sep=',')
data = data[["mol_id", "smiles", "absorption_energy", "emission_energy"]]

# -------- list of smiles strings --------
smiles = []
for i in data["smiles"]:
    smiles.append(i)
# -------- smiles to rdkit molecule --------
csv_list = []
for smile in smiles:
    molecule = Chem.MolFromSmiles(smile)
    molecule = Chem.AddHs(molecule)
    # -------- rdkit molecule positions --------
    rdDistGeom.EmbedMolecule(molecule)
    conf_mol = molecule.GetConformer()
    pos_mol = conf_mol.GetPositions()
    # -------- rdkit.atom obj to get symbols --------
    atom_mol = molecule.GetAtoms()
    atom_sym = []
    print(smile)
    for char in atom_mol:
        atom_sym.append(char.GetSymbol())
    print(atom_sym)
    atom = Atoms(atom_sym, positions=pos_mol)  # atom object to calculate triplet - singlet gap (in triplet geometry)
    atom.calc = TBLite(multiplicity=3)
    opt = BFGS(atom, logfile=f'opt.log_{smile}')
    opt.run(fmax=0.05)
    # atom = ase.io.read(filename='geometry.xyz')
    diff_energy = 0
    for mult_e in [1, 3]:
        atom.calc = TBLite(multiplicity=mult_e)
        if mult_e == 1:
            energy = atom.get_potential_energy()
            diff_energy -= energy
            print(f'{smile}', ' molecule singlet energy (triplet geometry)', ': %5.2f eV' % energy)
        if mult_e == 3:
            energy = atom.get_potential_energy()
            diff_energy += energy
            print(f'{smile}', ' molecule triplet energy (triplet geometry)', ': %5.2f eV' % energy)
    csv_list.append(diff_energy)
    print("Energy difference: ", diff_energy, "\n")
print(csv_list)

for element in csv_list:
    print(float(element))

