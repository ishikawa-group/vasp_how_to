# Adsorption and surface reaction
* In this section, we will see how to calculate the adsorption energy.
* When we considering the catalytic reaction using surface, adsorption is the first event.
* Chemical species adsorbed on a surface is called **adsorbates**.
* Adsorbates will undergo the **surface reaction**, which is a reaction that employes the surface. Usually, a surface reaction is easier to occur than the corresponding gas-phase reaction.
* There are two types of the surface reaction; one is called **Langmuir-Hinshelwood** type and the other is **Eley-Rideal** type. The former is the reaction that occurs *between the adsorbed species*. In the latter, the reaction occurs *between the adsorbed species and the gas-phase species*.
* Therefore, we need to consider following four types of reactions along the adsorption event (the asterisk $*$ means the surface, and species with asterisk is the adsorbate).

1. Adsorption
  ```math
  \ce{A + $*$ -> A$*$}
  ```
2. Surface reaction (Langmuir-Hinshelwood)
  ```math
  \ce{A$*$ + B$*$ -> C$*$}
  ```
3. Surface reaction (Eley-Rideal) 
  ```math
  \ce{A$*$ + D -> E}
  ```
4. Desorption
  ```math
  \ce{F$*$ -> F + $*$}
  ```

## Calculating adsorption energy
* As an example, we will see how to calculate the hydrogen adsorption energy on the Ru surface.
* This reaction is important, for example, in the ammonia (NH3) formation reaction.
* We use the ASE to build the bare, adsorbate, and adsorbed surface.

### ASE scripts
* The followings are the Python script using ASE.

#### 1. Simplest approach
```python{cmd}
from ase import Atom
from ase.io import write
from ase.build import hcp0001, add_adsorbate, molecule
from ase.visualize import view

surf = hcp0001(symbol="Ru", size=[3,3,4], a=2.7, c=4.2, vacuum=10.0, periodic=True)
clean_surface = surf.copy()

add_adsorbate(slab=surf, adsorbate=Atom("H"), height=1.2, position="fcc")

# output H2 POSCAR
h2 = molecule("H2")
h2.cell = [10,10,10]
h2.center()
write("POSCAR1", h2)

# output surface POSCAR
write("POSCAR2", clean_surface)

# output H-adsorbed surface POSCAR
write("POSCAR3", surf)
```
* Advantage: easy
* Disadvantage: non-generalizable (especally for high Miller index surfaces)

#### 2. General approach
```python
from ase import Atom
from ase.io import write
from ase.build import bulk, surface, add_adsorbate, molecule
from ase.visualize import view

bulk = bulk(name="Ru", crystalstructure="hcp", a=2.7, b=2.7, c=4.2)
surf = surface(lattice=bulk, indices=[0,0,1], layers=2, vacuum=10.0)

# make surface a little bigger
surf = surf*[3,3,1]

add_adsorbate(slab=surf, adsorbate=Atom("H"), height=1.2, offset=(0.22, 0.11))
...
```
* Advantage: a bit complex
* Disadvantage: general

#### Executing VASP
* After generating the POSCAR files (here POSCAR1, POSCAR2, and POSCAR3), the VASP calculation can be done.
* Do not forget to make directories for each calculations, setup the KPOINTS and POTCAR files (by copying VASP potentials files; see previous sections).
* Geometry optimization calculation should be done.
* After the calculation, get the energies by takeing the value of `energy(sigma->0)`.
* Then the adsorption energy ($E_{\rm ads}$) can be calculated as
```math
E_{\rm ads} = E_{\rm surf-ads} - \left(E_{\rm surf} + E_{\rm adsorbate}\right)
```
* Therefore in this case, the adsortion energy is
```math
E_{\rm ads} = E_\ce{Ru-H} - \left(E_\ce{Ru} + \frac{1}{2}E_\ce{H2}\right)
```
* The reference adsorption energy is -0.63 eV (fcc hollow). Compare with your calculated value.
  https://www.researchgate.net/figure/Adsorption-sites-for-H-on-the-Ru0001-surface-with-corresponding-energies-of-adsorption_fig1_340125395


### Advanced
* You can load several files using `ase.io.read`, to make a variety of surface-adsorbate pairs.
* In the following, we consider the benzene adsorption on the Si surface.
* Assuming that we have the coordinate of benzene ("benzene.xyz") and that of bulk Si ("si.cif"). There are available from many sources (such as PubChem, Materials Project).

```python
from ase import Atom
from ase.io import read, write
from ase.build import bulk, surface, add_adsorbate, molecule
from ase.visualize import view

bulk = read("./files/si.cif")
surf = surface(lattice=bulk, indices=[0,0,1], layers=2, vacuum=10.0)
surf = surf*[2,2,1]

adsorbate = read("benzene.xyz")
add_adsorbate(slab=surf, adsorbate=adsorbate, height=2.0, offset=(0.5, 0.5))

view(surf)
```

## Exercise
1. Calculate the adsorption energy on other sites. The "hollow" keyword in `add_adsorbate` for hcp0001 puts adsorbate on the fcc hollow site. Use other keywords such as "ontop" or "bridge", and compare with the reference.
2. Adjust the `offset` value in the `add_adsorbate`, and put the H atom to the hcp hollow site. Then calculate the $E_{\rm ads}$.