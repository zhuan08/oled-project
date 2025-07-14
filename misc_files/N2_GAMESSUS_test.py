import ase
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.gamess_us import GAMESSUS
from ase.io import read, write
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_3d = True

# -------- Comparing Energies based on Positions --------
N2 = Atoms('N2', [(0., 0., 0.), (0., 0., 1.1)])
N2.calc = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
opt = BFGS(N2, logfile='opt.log')
opt.run(fmax=0.05)
ase.io.write(filename='geometry.xyz', images=N2, format='xyz')

for mult_e in [1, 3]:
    N2.calc = GAMESSUS(contrl={'mult': mult_e}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
    if mult_e == 1:
        energy = N2.get_potential_energy()
        print('N2, molecule singlet energy (triplet geometry)', ': %5.2f eV' % energy)
    if mult_e == 3:
        energy = N2.get_potential_energy()
        print('N2, molecule triplet energy (triplet geometry)', ': %5.2f eV' % energy)
