import ase
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

# -------- Methanol --------
inpath = 'misc_files/C1H4O1.xyz'
calc = TBLite()
atoms = ase.io.read(inpath)
atoms.calc = calc
energy = atoms.get_potential_energy()
print(f'Energy: {energy} eV')
print(atoms.get_positions())
