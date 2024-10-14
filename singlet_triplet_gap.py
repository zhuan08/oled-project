# import ase
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.gamess_us import GAMESSUS
# from ase.visualize import view
# from ase.io import read
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

# -------- Constants --------
# h = 6.626*(10**(-34))
# c = 2.998*(10**8)
# energy = 3.3*(1.6*10**(-19))
# wavelength = "{:2e}".format((c*h)/energy)
# print('wavelength conversion: ', f'{wavelength}m')

# -------- Read csv files for smiles strings --------
data = pd.read_csv('smiles_energy.csv', sep=',')
data = data[["mol_id", "smiles", "absorption_energy", "emission_energy"]]

# -------- list of smiles strings --------
smiles = []
for i in data["smiles"]:
    smiles.append(i)
# -------- smiles to rdkit molecule --------
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
    atom.calc = GAMESSUS(contrl={'mult': 3, 'scftyp': 'ROHF'},
                         label='geometry_op',
                         command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
    opt = BFGS(atom, logfile=f'opt.log_{smile}')
    opt.run(fmax=0.05)
    # atom = ase.io.read(filename='geometry.xyz')
    diff_energy = 0
    for mult_e in [1, 3]:
        atom.calc = GAMESSUS(contrl={'mult': mult_e, 'scftyp': 'UHF', 'runtyp': 'energy',
                                     'dfttyp': 'B3LYP', 'maxit': 200},
                             label=f'molecule{mult_e}',
                             system={'mwords': 50},
                             basis={'gbasis': 'STO', 'ngauss': 6, 'npfunc': 3, 'npfunc': 3, 'diffsp': True},
                             command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
        if mult_e == 1:
            energy = atom.get_potential_energy()
            diff_energy -= energy
            print(f'{smile}', ' molecule singlet energy (triplet geometry)', ': %5.2f eV' % energy)
        if mult_e == 3:
            energy = atom.get_potential_energy()
            diff_energy += energy
            print(f'{smile}', ' molecule triplet energy (triplet geometry)', ': %5.2f eV' % energy)
    print("Energy difference: ", diff_energy, "\n")
# view(atom)
