import glob
import ase.io
from ase.calculators.gamess_us import GAMESSUS
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_3d = True

path_id = []
geom_id = []
diff_energy_list = []

for geom_path in glob.glob('organic_phos_geometries/*.xyz'):
    path_id.append(geom_path)
print(path_id)

for path in path_id:
    atom = ase.io.read(filename=path)
    diff_energy = 0
    for mult_e in [1, 3]:
        atom.calc = GAMESSUS(contrl={'mult': mult_e}, label='molecule',
                       command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
        if mult_e == 1:
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if mult_e == 3:
            energy = atom.get_potential_energy()
            diff_energy += energy
    diff_energy_list.append(diff_energy)

for geom_path in glob.glob('organic_phos_geometries/*.xyz'):
    if '.xyz' in geom_path:
        updated_file = geom_path.replace('.xyz', '')
        updated_file = updated_file.replace('organic_phos_geometries', '')
        updated_file = updated_file.replace('/', '')
        geom_id.append(updated_file)

for i in range(len(geom_id)):
    print(geom_id[i], diff_energy_list[i])