from ase import Atoms
from ase.visualize import view
from tblite.ase import TBLite
from tblite.ase import GAMESSUS

tblite_calc_sing = TBLite(multiplicity=1)
tblite_calc_trip = TBLite(multiplicity=3)
gamess_calc_sing = GAMESSUS(contrl={'mult': 1}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
gamess_calc_trip = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')

calculator = tblite_calc_sing

# -------- Nitrogen Atom 1 --------
N = Atoms('N')
N.calc = TBLite()
E_N = N.get_potential_energy()

# -------- Nitrogen Atom 2 --------
N2 = Atoms('N2', [(0., 0., 0.), (0., 0., 1.1)])
N2.calc = calculator
E_N2 = N2.get_potential_energy()

E_atomization = E_N2 - 2 * E_N

print('Nitrogen atom energy: %5.2f eV' % E_N)
print('Nitrogen molecule energy: %5.2f eV' % E_N2)
print('Atomization energy: %5.2f eV' % -E_atomization)
view(N2)
