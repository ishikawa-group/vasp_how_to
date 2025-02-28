# Geometry optimization
* In this section, we learn how to do geometry optimization with VASP.
* Geometry optimization is also called *ion relaxation*, as it is the procedure to move the ion (or nuclei) position to have *lower potential energy*.
* Thus the *geometry* means the configuration space spanned by many nuclei.
* Lower potential energy means molecules, bulk, or surfaces are more stable than the high-potential-energy state. Thus we can consider the stable molecules or bulk (observed experimentally) have the optimized structures.
* The optimized geometry thus can be compared with exprimental molecular or solid strucutures, which is identified by the X-ray or nuclear magnetic resonance (NMR) methods.

## Optimization of bulk
### POSCAR
* POSCAR is the same with the energy calculation.
* Let's take the Pt example ([file](./Pt_POSCAR)).
```
Pt 
 1.0000000
     3.9231600    0.0000000    0.0000000
     0.0000000    3.9231600    0.0000000
     0.0000000    0.0000000    3.9231600
 Pt 
   4
Cartesian
  0.0000000  0.0000000  0.0000000
  1.9615800  1.9615800  0.0000000
  0.0000000  1.9615800  1.9615800
  1.9615800  0.0000000  1.9615800
```

### INCAR
* The INCAR file for the geometry optimization is as follows.
```
SYSTEM = Pt bulk
ISMEAR =  1
NSW    = 10
IBRION =  2
```
* `IBRION = 1 or 2` specifies to do geometry optimization. The difference of 1 and 2 is the optimization algorithm.
  + `IBRION = 1`: RMM-DIIS algorithm
  + `IBRION = 2`: Congugate-Gradient algorithm
  + Usually, it is recommended to use 2.
  + For conjugate-gradient, see: https://qiita.com/Dason08/items/27559e192a6a977dd5e5
* `NSW` is the maximum number of optimization step. Usually, set this value to 10-200 steps.

### KPOINTS and POTCAR
* These files should be set as usual. Since this this the bulk calculation, it is better to set k-points = 3x3x3 or larger for meaningful calculation.
* POTCAR for Pt atom should be set.

### Execure
* Execute VASP calculation as usual.

### Checking output
* The simplest method to see the optimization process is to use ASE; `ase gui vasprun.xml` or `ase gui OUTCAR`.
* If the GUI of your remote environment (TSUBAME etc.) is slow, download your `vasprun.xml` or `OUTCAR` file to your local computer (via `scp` or any) then use `ase gui vasprun.xml` or `ase gui OUTCAR`.
* By using above, **you can see the atoms are moving**. This is what the geometry optimizaitons are doing: the atoms are relaxed, in order to lower the energy.
* Also check the energy in OUTCAR (the latter value of the following line):
`energy without entropy = -1070.27460902  energy(sigma->0) = -1070.45373832`
    + `grep` is useful for this purpose: `grep "energy  without entropy" OUTCAR`
* You can see the energy is (mostly) decreasing. If it reached the minimum, the optimization is finished.

## Unit cell optimization
* To change the unit cell size and shape, insert the following tag in `INCAR`.
  + `ISIF = 3`: changing cell volume, shape, and atomic positions
  + `ISIF = 4`: changing cell shape and atomic positions
  + Default `ISIF` is 2 (changing atomic positions)
  + Details can be seen in https://www.vasp.at/wiki/index.php/ISIF

## Exercise
* Download the FCC metals (nickel, platinum, silver) and do the cell optimization. Compare the calculated lattice constants.

---

## Optimization of surface
* The geometry optimization of surface is mostly same with the bulk optimization, but several differences existss.

### Fixing the surface
* In surface optimization, part of the surface (or called slab) is fixed or frozen from its initial position.
* This is because we often mimic the surface material by fixing (or freezing) the lower part of the system. This is reasonable, because, in reallity, only surface atoms moves significantly and lower part have similar structure with bulk material.
* To do such calculation, we fix the lower part of the surface model and only relax the upper layers. This can be done with POSCAR setting.
* Let's take the Pt surface as example ([file](./Pt_surface_POSCAR)). This is made by just changing the unit cell size by extending the z-direction (check it by `ase gui` or VESTA).
```
Pt 
 1.0000000
     3.9231600    0.0000000    0.0000000
     0.0000000    3.9231600    0.0000000
     0.0000000    0.0000000   20.0000000
 Pt 
   4
Cartesian
  0.0000000  0.0000000  0.0000000
  1.9615800  1.9615800  0.0000000
  0.0000000  1.9615800  1.9615800
  1.9615800  0.0000000  1.9615800
```
* Then, let's change it as follows:
```
Pt 
 1.0000000
     3.9231600    0.0000000    0.0000000
     0.0000000    3.9231600    0.0000000
     0.0000000    0.0000000   20.0000000
 Pt 
   4
Selective dynamics
Cartesian
  0.0000000  0.0000000  0.0000000   F   F   F
  1.9615800  1.9615800  0.0000000   F   F   F
  0.0000000  1.9615800  1.9615800   T   T   T
  1.9615800  0.0000000  1.9615800   T   T   T
```
* The `Selective dynamics` indicates we have both fixed and relax atoms.
* `F F F` or `T T T` of each atom is the answer to the question: *Do you want this atom to move?*. `F` means *False* thus does not move, and `T` means *True* so it allows to move atom.
* `ase gui` of this POSCAR distinguishes fixed and free atoms (with check mark for atoms) while vesta does not make distinction.
* After changing POSCAR, execute VASP and see the relaxed structure.

#### Using ASE (Advanced topic)
* With ASE, fixing the surface can be done by setting `constraint` property of `Atoms`.
* To use this function, `tag` of Atoms should be used; tags are automatically set when surface is made with ASE.
  + adsorbate: tag = 0, surface layer: tag = 1, etc.
  + Lower (i.e. smaller in z-coordinate) surface layers have large tag values.
* By using this tag and `ase.constraints.FixAtoms`, the constraint property can be set.
```python
from ase.constraints import FixAtoms
from ase.build import fcc111
from ase.io import write
from ase.visualize import view

surf = fcc111(symbol='Pd', size=[3, 3, 4], a=3.99, vacuum=10.0)
constraints = FixAtoms(indices=[atom.index for atom in surf if atom.tag >= 3])
surf.constraints = constraints
write("POSCAR", surf)
```

## Exercise
* Confirm the optimization process by `ase gui`.
