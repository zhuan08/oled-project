import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Draw
from ase.visualize import view
import pandas as pd
IPythonConsole.ipython_3d = True

files = ['CUYRAS', 'CUYZZ03', 'CUYZZ07', 'CUYZZ08', 'CUYZZ09', 'IJUZZ04', 'LOKZZ01', 'QAYZZ04', 'UVEFAE']
for name in files:
    atom = ase.io.read(filename=f'geometries/{name}.xyz')
    view(atom)
