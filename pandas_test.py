from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Draw
import pandas as pd
IPythonConsole.ipython_3d = True

data = pd.read_csv('oled_smiles_energy.csv', sep=',')
data = data[["mol_id", "smiles", "absorption_energy", "emission_energy"]]

# print(data["mol_id"])
for i in data["mol_id"]:
    print(i)

for i in data["emission_energy"]:
    print(i)
# print(data[["emission_energy"]])