from ase.build import fcc111
from ase.visualize import view

surf = fcc111(symbol="Ru", size=[2,2,2], a=2.7, vacuum=10.0)
view(surf)
