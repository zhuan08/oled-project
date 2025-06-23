import os
import ase.io
from octahedral_embed import octahedral_embed
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Draw
import pandas as pd
import numpy as np
IPythonConsole.ipython_3d = True