from ase import Atoms
from ase.visualize import view
from ase.calculators.gamess_us import GAMESSUS

# -------- Nitrogen Atom 1 --------
N = Atoms('N')
N.calc = GAMESSUS(contrl={'mult': 2}, label='N', command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
E_N = N.get_potential_energy()

# -------- Nitrogen Atom 2 --------
N2 = Atoms('N2', [(0., 0., 0.), (0., 0., 1.1)])
N2.calc = GAMESSUS(contrl={'mult': 1}, label='molecule', command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
E_N2 = N2.get_potential_energy()

E_atomization = E_N2 - 2 * E_N

print('Nitrogen atom energy: %5.2f eV' % E_N)
print('Nitrogen molecule energy: %5.2f eV' % E_N2)
print('Atomization energy: %5.2f eV' % -E_atomization)
view(N2)
