from ase import Atoms
from ase.build import hcp0001, add_adsorbate
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.visualize import view

# Set up the Ru(0001) surface
slab = hcp0001('Ru', size=(2, 2, 3), vacuum=10.0)

# Set up the N atom on the surface
adsorbate = Atoms('N')
add_adsorbate(slab, adsorbate, height=1.8, position='ontop')

# Set up the calculator for the N/Ru(0001) system
calc = Vasp(xc='PBE', encut=400, kpts=(4, 4, 1), ediff=1e-5, ediffg=-0.02, ibrion=2, nsw=100, isif=2, lwave=False, lcharg=False)
slab.set_calculator(calc)

# Optimize the structure
opt = BFGS(slab)
opt.run(fmax=0.02)

# Get the total energy of the N/Ru(0001) system
E_N_Ru0001 = slab.get_potential_energy()

# Calculate the energy of the clean Ru(0001) surface
clean_slab = fcc111('Ru', size=(2, 2, 3), vacuum=10.0)
clean_slab.set_calculator(calc)
E_Ru0001 = clean_slab.get_potential_energy()

# Calculate the energy of the isolated N atom
N_atom = Atoms('N')
N_atom.set_calculator(calc)
E_N = N_atom.get_potential_energy()

# Calculate the adsorption energy
E_ads = E_N_Ru0001 - E_Ru0001 - E_N

print(f'Adsorption energy: {E_ads} eV')