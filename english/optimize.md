# Geometry optimization
* In this section, we learn how to do geometry optimization with VASP.
* Geometry optimization is also called *ion relaxation*, as it is the procedure to move the ion (or nuclei) position to have *lower potential energy*.
* Lower potential energy means molecules, bulk, or surfaces are more stable than the high-potential-energy state. Thus we can consider the stable molecules or bulk (observed experimentally) have the optimized structures.
* The optimized geometry thus can be compared with exprimental molecular or solid strucutures, which is identified by the X-ray or NMR methods.

## POSCAR
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

## INCAR
* The INCAR file for the geometry optimization is as follows.
```
SYSTEM = Pt bulk
ISMEAR =  1
NSW    = 10
IBRION =  2
```
* `IBRION = 1 or 2` specifies to do geometry optimization. The difference of 1 and 2 is the optimization algorithm. Usually, it is recommended to use 2.
* `NSW` is the maximum number of optimization step. Usually, set this value to 10-200 steps.

## KPOINTS and POTCAR
* These files should be set as usual. Since this this the bulk calculation, it is better to set k-points = 3x3x3 or larger for meaningful calculation.
* POTCAR for Pt atom should be set.

## Execure
* Execute VASP calculation as usual.

## Checing output
* The simplest method to see the optimization process is to use ASE; `ase gui vasprun.xml` or `ase gui OUTCAR`.
* If the GUI of your remote environment (TSUBAME etc.) is slow, download your vasprun.xml or OUTCAR file to your local computer (via `scp` or any) then use `ase gui`.

# Geometry optimization of surface
* The geometry optimization of surface is different from those of molecule and bulk materials.
* This is because we often mimic the surface material by fixing (or freezing) the lower part of the system. This is reasonable, because, in reallity, only surface atoms moves significantly and lower part have similar structure with bulk material.
* To do such calculation, we fix the lower part of the surface model and only relax the upper layers. This can be done with POSCAR setting.
* Let's take the Pt surface as example ([file](./Pt_surface_POSCAR)). This is made by just changing the unit cell size by extending the z-direction (check it by `ase gui` or vesta).
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

## Adsorption
* 