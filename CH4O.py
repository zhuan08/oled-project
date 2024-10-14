from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.gamess_us import GAMESSUS
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

# -------- Methanol --------
# inpath = 'C1H4O1.xyz'
# calc = GAMESSUS(contrl={'mult': 3}, label='methanol',
#                 command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
# atoms = ase.io.read(inpath)
# atoms.calc = calc
# energy = atoms.get_potential_energy()
# print(f'Energy: {energy} eV')
# print(atoms.get_positions())
