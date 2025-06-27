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

geom_id = []
diff_energy_list = []

for geom_path in glob.glob('geometries/*.xyz'):
    geom_id.append(geom_path)
print(geom_id)

for path in geom_id:
    atom = ase.io.read(filename=path)
    diff_energy = 0
    for mult_e in [1, 3]:
        atom.calc = TBLite(multiplicity=mult_e)
        if mult_e == 1:
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if mult_e == 3:
            energy = atom.get_potential_energy()
            diff_energy += energy
    diff_energy_list.append(diff_energy)

for i in range(len(geom_id)):
    print(geom_id[i], diff_energy_list[i])