# Joshua Ashby 6/22/23 this program performs NEB on a system
# where an atom is removed and a nearby atom is changed
# to Cu and moved to the other atom's nearby position
from ase.build import fcc111, add_adsorbate
from ase.spacegroup import crystal
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, ExpCellFilter
from ase.build import make_supercell
from ase.optimize import QuasiNewton, GPMin
import matplotlib.pyplot as plt
from ase.optimize import BFGS
from ase.io import write
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
from ase.constraints import Filter, ExpCellFilter
from ase.optimize import BFGS
from ase.build import fcc110
from ase.calculators.vasp import Vasp
from ase.calculators.lj import LennardJones
from ase import Atoms, Atom
from ase.build import fcc111, add_adsorbate
from ase.spacegroup import crystal
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms, ExpCellFilter
from ase.build import make_supercell
from ase.optimize import QuasiNewton, GPMin
import matplotlib.pyplot as plt
from ase.optimize import BFGS
from ase.io import write
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
from ase.constraints import Filter, ExpCellFilter
from ase.optimize import BFGS
from ase.build import fcc110
from ase.calculators.vasp import Vasp
from ase.calculators.lj import LennardJones
import numpy as np

from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize.fire import FIRE as QuasiNewton
from ase.io import read

calc1 = Vasp(
    encut=400.000000,
    sigma=0.100000,
    ismear=0,
    ediff=1.00e-06,
    algo='Normal',
    pp='PBE',
    nelmin=5,
    lasph='.TRUE.',
    ivdw=13,
    ispin=2,
    nelm=200,
    ncore=8,
    kpts=[1, 1, 1],
    ldau='.TRUE.',
    ldauj=[0, 0, 0],
    ldaul=[2, 0, 0],
    ldautype=2,
    ldauu=[3.1, 0, 0],
    lmaxmix=4
    )

slab = read('LI_POS_INIT.vasp')
slab.calc = calc1
slab[168].position += [0.15,0.18,0.73]
slab.rattle(0.01)
print('energy:', slab.get_potential_energy())
print(slab.get_positions())