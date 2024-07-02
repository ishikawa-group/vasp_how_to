# Energy calculation (bulk and surface)
* In this section, we will see how to do the energy calculation of the bulk material and the surface of the bulk material.
* Please see [the energy calculation for molecules](./energy_molecule.md) if not yet.
* The bulk material is crystal surface, which is periodic in three dimensions (i.e. x-, y-, and z-directions).
* If you want to do calculations for surfaces, the situation is a little different because the surface loses one periodic direction. Usually this axis is set as z-direction.
* To cut the periodicity of z-direction, a vacant region called **vacuum layer** in inserted in the surface calculaiton.

# Downloading CIF files
* To do the bulk/surface calcuation, Crystal information file (CIF) file is often needed.
* The CIF file can be downloaded from **Materials project**, **ICSD**, or others.
* Here, a simple CIF file (`Pt.cif`) is prepared so please use it.
* The CIF file contains the information of atomic position and unit cell information.

# Making POSCAR file
## VESTA
* VESTA is the free software for visualization and editing the molecular/bulk/surface structure: https://jp-minerals.org/vesta/jp/
* You can visualize CIF, POSCAR, xyz (and other) files with VESTA.

### Bulk calculation
* To make the POSCAR file for bulk calculation, open the CIF file and then
    1. `File` -> `Export Data` and then choose `VASP file` in `File type`.
    2. Any name is OK.
    3. Either "fractional coordinate" or "Cartesian coordinate" is OK.
    4. Rename the file to `POSCAR`.

### Surface calculation
* It is possible to make the surface file with VESTA, but it is rather complicated. So we will use Atomic simulation environment (ASE) instead.

#### ASE
* To make the POSCAR file for surface (e.g. fcc 111 surface), write and execute the following Python command (or make a `.py` script and execute).
    1. `python`
    2. `>>> from ase.io import read, write`
    3. `bulk = read("Pt.cif")`
    4. `write("POSCAR", bulk)`
    5. `quit()`
* You can check the structure with VESTA, or using `ase gui POSCAR` if you've installed the ASE.

## KPOINTS
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

## POTCAR
* Do not forget to make POTCAR file by `cat {POTCAR_directory}/Pt/POTCAR > POTCAR`.

# Executing VASP
* Execution of VASP is the same with the energy calculation of molecule.
    + edit `run.sh` script (or any name)
    + `qsub run.sh`

## Exercise
* In the Pt bulk case using ASE, `a=3.92` in `fcc111` means the lattice constant (size of the unit cell) is 3.92 Angstrom. Change this value to other value (like 3.5 or 4.5) and compare the calculated energy with that of 3.92 Angstrom.