# Energy calculation (bulk and surface)
* In this section, we will see how to do the energy calculation of the bulk material and the surface of the bulk material.
* Please see *the energy calculation for molecules* first.
* The bulk material is crystal surface, which is periodic in three dimensions (i.e. x-, y-, and z-directions).
* On the other hand, the surface loses one periodic direction. Usually this axis is set as z-direction.
* To cut the periodicity of z-direction, a vacant region called **vacuum layer** in inserted in the surface calculaiton.

# Downloading CIF files
* To do the bulk/surface calcuation, Crystal information file (CIF) file is often needed.
* The CIF file can be downloaded from **Materials project**, **ICSD**, or others.
* The CIF file contains the information of atomic position and unit cell information.

# Making POSCAR file
## VESTA
* VESTA is the free software for visualization and editing the molecular/bulk/surface structure.
* You can visualize CIF, POSCAR, xyz (and other) files with VESTA.
* The surface file can be made with VESTA, but it is rather complicated. So we will use Atomic simulation environment (ASE) instead.

## ASE
* To make the POSCAR file for surface (e.g. fcc 111 surface), write and execute the following python script.
    ```python{cmd}
    from ase.build import fcc111
    from ase.io import write

    surf = fcc111("Pt", a=3.92, size=[3, 3, 3], vacuum=10.0)
    write("POSCAR", surf)
    ```
* You can check the structure with VESTA, or using `ase gui` if you've installed the ASE.

## k-ponits
* For the bulk calculation,
```
bulk
0
Monkhorst Pack
3 3 3
0 0 0
```
* For the surface calculation, just one k-point is OK for the non-periodic direction. So when the vaccum layer is in the z-direction, following is fine.
```
surface
0
Monkhorst Pack
3 3 1
0 0 0
```
* Using more k-points gives higher accuracy, but larger computational cost.

# Executing VASP
* Execution of VASP is the same with the energy calculation of molecule.
    + edit `run.sh` script (or any name)
    + `qsub run.sh`

## Exercise
* In the Pt bulk case using ASE, `a=3.92` in `fcc111` means the lattice constant (size of the unit cell) is 3.92 Angstrom. Change this value to other value (like 3.5 or 4.5) and compare the calculated energy with that of 3.92 Angstrom.