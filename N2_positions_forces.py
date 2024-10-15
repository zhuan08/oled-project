from ase import Atoms
from ase.visualize import view
from ase.calculators.gamess_us import GAMESSUS
from ase.optimize import BFGS
from rdkit import Chem

# -------- Nitrogen Atom --------
# N2 = Atoms('N2', [(0., 0., 0.), (0., 0., 0.8)])
# N2.calc = GAMESSUS(contrl={'mult': 1}, label='molecule',
#                    command='rungms PREFIX.inp 2023.R1.intel > PREFIX.log 2> PREFIX.err')
# P_N2 = N2.get_positions()
# F_N2 = N2.get_forces()
# print("Positions: \n", P_N2)
# print("Forces: \n", F_N2)
#
# opt = BFGS(N2, trajectory='opt.traj', logfile='opt.log')
# opt.run(fmax=0.05)
# P_N2 = N2.get_positions()
# F_N2 = N2.get_forces()
# print("Positions: \n", P_N2)
# print("Forces: \n", F_N2)