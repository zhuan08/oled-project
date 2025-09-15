from ase import Atoms
from ase.visualize import view
from tblite.ase import TBLite
from ase.optimize import BFGS
from rdkit import Chem

# -------- Nitrogen Atom --------
N2 = Atoms('N2', [(0., 0., 0.), (0., 0., 0.8)])
N2.calc = TBLite()
P_N2 = N2.get_positions()
F_N2 = N2.get_forces()
print("Positions: \n", P_N2)
print("Forces: \n", F_N2)

opt = BFGS(N2, trajectory='opt.traj', logfile='opt.log')
opt.run(fmax=0.05)
P_N2 = N2.get_positions()
F_N2 = N2.get_forces()
print("Positions: \n", P_N2)
print("Forces: \n", F_N2)
