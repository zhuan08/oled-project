from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
IPythonConsole.ipython_3d = True
import numpy as np

# -------- Working with Conformers --------
# nitrogen = Chem.MolFromSmiles('N#N')
# print(nitrogen.GetNumConformers())
# rdDepictor.Compute2DCoords(nitrogen) # <- generated 2D coordinates, added as conformer
# print(nitrogen.GetNumConformers())
# print(nitrogen.GetConformer().Is3D()) # <- flagged for not being 3D

# -------- Generating 3D Structure --------
# nitrogen = Chem.AddHs(nitrogen)
# rdDistGeom.EmbedMolecule(nitrogen)
# print(nitrogen.GetNumConformers(),nitrogen.GetConformer().Is3D())

# -------- Multiple Conformers One Molecule --------
# rdDistGeom.EmbedMultipleConfs(nitrogen,10, randomSeed=0xf00d) # <- generate 10 conformers
# print(nitrogen.GetNumConformers(),nitrogen.GetConformer().Is3D())

# -------- Marking 3D conformer as 2D --------
# tmol = Chem.Mol(nitrogen)
# tmol.GetConformer().Set3D(False)
# Draw.MolToImage(tmol).show()

# -------- Conformer IDs --------
# print([x.GetId() for x in nitrogen.GetConformers()])
# cp = Chem.Mol(nitrogen)
# cp.RemoveConformer(0)
# cp.RemoveConformer(3)
# cp.RemoveConformer(7)

# -------- Default Conformer is 1 --------
# print([x.GetId() for x in cp.GetConformers()])
# print((nitrogen.GetConformer().GetId(), cp.GetConformer().GetId()))
# cp.GetConformer(12) # <- error if conformer DNE

# -------- Getting Atom Positions --------
# conf = nitrogen.GetConformer()
# pos = conf.GetAtomPosition(0)
# print(pos)
# print(pos.x) # <- x,y,z, axes
# print(pos[0]) # <- vectors
# print(list(pos)) # <- converted to list

# -------- Get Positions --------
# conf = nitrogen.GetConformer()
# print(conf.GetPositions())

# -------- list of atoms --------
# [N#N, O=O, Cl-Cl]
# -------- list of geometry --------
# [positions of N#N, O=O, Cl-Cl]
# -------- list of multiplicity for energy --------
# [mult=1, mult=3]

nitrogen = Chem.MolFromSmiles('N#N')
nitrogen = Chem.AddHs(nitrogen)
rdDistGeom.EmbedMolecule(nitrogen)
conf_n = nitrogen.GetConformer()
n_pos = conf_n.GetPositions()

oxygen = Chem.MolFromSmiles('O=O')
oxygen = Chem.AddHs(oxygen)
rdDistGeom.EmbedMolecule(oxygen)
conf_o = nitrogen.GetConformer()
o_pos = conf_n.GetPositions()

chlorine = Chem.MolFromSmiles('Cl-Cl')
chlorine = Chem.AddHs(chlorine)
rdDistGeom.EmbedMolecule(chlorine)
conf_cl = nitrogen.GetConformer()
cl_pos = conf_n.GetPositions()

elements = ['N2', 'O2', 'Cl2']
positions = [n_pos, o_pos, cl_pos]
h = 6.626*(10**(-34))
c = 2.998*(10**8)
print("Singlet Triplet Energy Calculations")
print("-----------------------------------")
for element in elements:
    print(f'[{element}]')
    for position in positions:
        molecule = Atoms(symbols=f'{element}', positions=position)
        for mult_g in [1, 3]:
            molecule.calc = TBLite(multiplicity=mult_g)
            opt = BFGS(molecule, logfile='opt.log')
            opt.run(fmax=0.05)
            molecule_pos = molecule.get_positions()
            diff_energy = 0
            print("Positions/Geometry (mult =", f'{mult_g}) \n', molecule_pos)
            for mult_e in [1, 3]:
                molecule.calc = TBLite(multiplicity=mult_e)
                if mult_g == 1 and mult_e == 1:
                    energy = molecule.get_potential_energy()
                    print(f'{element}', ' molecule singlet energy (mult =', f'{mult_e})', ': %5.2f eV' % energy)
                elif mult_g == 1 and mult_e == 3:
                    energy = molecule.get_potential_energy()
                    print(f'{element}', ' molecule triplet energy (mult =', f'{mult_e})', ': %5.2f eV' % energy)
                elif mult_g == 3 and mult_e == 1:
                    energy = molecule.get_potential_energy()
                    diff_energy -= energy
                    print(f'{element}', ' molecule singlet energy (mult =', f'{mult_e})', ': %5.2f eV' % energy)
                elif mult_g == 3 and mult_e == 3:
                    energy = molecule.get_potential_energy()
                    diff_energy += energy
                    lam = "{:2e}".format((h * c) / (diff_energy * (1.6 * 10 ** (-19))))
                    print(f'{element}', ' molecule triplet energy (mult =', f'{mult_e})', ': %5.2f eV' % energy)
                    print(f'{element}', ' energy difference: %5.2f eV' % diff_energy)
                    print(f'{element}', ' wavelength conversion: ', f'{lam}m \n')
print("------------------------------------")
