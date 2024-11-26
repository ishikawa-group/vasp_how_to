from ase import Atom
from ase.io import write
from ase.build import hcp0001, add_adsorbate, molecule
from ase.visualize import view

surf = hcp0001(symbol="Ru", size=[3, 3, 4], a=2.7, c=4.2, vacuum=10.0, periodic=True)
clean_surface = surf.copy()

add_adsorbate(slab=surf, adsorbate=Atom("H"), height=1.2, position="fcc")

# output H2 POSCAR
h2 = molecule("H2")
h2.cell = [10, 10, 10]
h2.center()
write("POSCAR1", h2)

# output surface POSCAR
write("POSCAR2", clean_surface)

# output H-adsorbed surface POSCAR
write("POSCAR3", surf)