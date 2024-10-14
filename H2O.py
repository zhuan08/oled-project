from ase import Atoms
from ase.visualize import view
from ase.calculators.gamess_us import GAMESSUS

# -------- Hydrogen Gas --------
H2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.7414)])
H2.calc = GAMESSUS(contrl={'mult': 1}, label='h2', command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
E_H2 = H2.get_potential_energy()

# -------- Oxygen Gas --------
O2 = Atoms('O2', [(0, 0, 0), (0, 0, 1.2075)])
O2.calc = GAMESSUS(contrl={'mult': 3}, label='o2', command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
E_O2 = O2.get_potential_energy()

# -------- Water Atom --------
H2O = Atoms('H2O', [(0, 0.7572, -0.4692), (0, -0.7575, -0.4692), (0, 0, 0.1173)])
H2O.calc = GAMESSUS(contrl={'mult': 1}, label='h2o', command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
E_H2O = H2O.get_potential_energy()

E_combustion = (E_H2O - E_H2 - 0.5*E_O2)

print(f'Hydrogen gas energy: {E_H2:.2f} eV')
print(f'Oxygen gas energy: {E_O2:.2f} eV')
print(f'Water molecule energy: {E_H2O:.2f} eV')
print(f'Combustion energy: {E_combustion:.2f} eV')
# view(H2O_molecule)
