import glob
import ase.io
import pandas as pd
from pandas import merge_ordered
from tblite.ase import TBLite
from ase.calculators.gamess_us import GAMESSUS
from ase.calculators.psi4 import Psi4
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_3d = True

path_id = []
geom_id = []
diff_energy_list = []

data = pd.read_csv('smiles_energy.csv', sep=',')
data = data[["mol_id", "smiles", "absorption_energy", "emission_energy"]]
actual_energy = pd.DataFrame({"mol_id": data["mol_id"],
                              "actual_energy": data["emission_energy"]})

tblite_calc_singlet = TBLite(multiplicity=1)
tblite_calc_triplet = TBLite(multiplicity=3)
gamess_calc_singlet = GAMESSUS(contrl={'mult': 1}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
gamess_calc_triplet = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')

CHARGE = 0
METHOD = 'b3lyp'
BASIS = 'def2-svp'
psi4_calc_singlet = Psi4(method=METHOD, basis=BASIS, charge=CHARGE, multiplicity=1)
psi4_calc_triplet = Psi4(method=METHOD, basis=BASIS, charge=CHARGE,
                  multiplicity=3, reference='uks')

for geom_path in glob.glob('organic_phos_geometries/*.xyz'):
    path_id.append(geom_path)

for path in path_id:
    atom = ase.io.read(filename=path)
    diff_energy = 0
    for multiplicity in ['singlet', 'triplet']:
        if multiplicity == 'singlet':
            atom.calc = psi4_calc_singlet
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if multiplicity == 'triplet':
            atom.calc = psi4_calc_triplet
            energy = atom.get_potential_energy()
            diff_energy += energy
    diff_energy_list.append(diff_energy)

for geom_path in glob.glob('organic_phos_geometries/*.xyz'):
    if '.xyz' in geom_path:
        updated_file = geom_path.replace('.xyz', '')
        updated_file = updated_file.replace('organic_phos_geometries', '')
        updated_file = updated_file.replace('/', '')
        geom_id.append(updated_file)
        
organic_phos_data = pd.DataFrame({'mol_id': geom_id,
                     'predicted_energy': diff_energy_list})
merged = (organic_phos_data.merge(actual_energy))
merged.to_csv('organic_phos_data', index=False)
