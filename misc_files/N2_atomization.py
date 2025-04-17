from ase import Atoms
from ase.visualize import view
from tblite.ase import TBLite

# -------- Nitrogen Atom 1 --------
N = Atoms('N')
N.calc = TBLite()
E_N = N.get_potential_energy()

# -------- Nitrogen Atom 2 --------
N2 = Atoms('N2', [(0., 0., 0.), (0., 0., 1.1)])
N2.calc = TBLite()
E_N2 = N2.get_potential_energy()

E_atomization = E_N2 - 2 * E_N

print('Nitrogen atom energy: %5.2f eV' % E_N)
print('Nitrogen molecule energy: %5.2f eV' % E_N2)
print('Atomization energy: %5.2f eV' % -E_atomization)
view(N2)
