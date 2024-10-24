# Energy calculation (bulk and surface)
* In this section, we will see how to do the energy calculation of the bulk material and the surface of the bulk material.
* Please see [the energy calculation for molecules](./energy_molecule.md) if not yet.
* The bulk material is crystal surface, which is periodic in three dimensions (i.e. x-, y-, and z-directions).

## Downloading CIF files
* To do the bulk/surface calcuation, Crystal information file (CIF) file is often needed.
* The CIF file can be downloaded from **Materials project**, **ICSD**, or others.
* Here, a simple CIF file (`Pt.cif`) is prepared so please use it.
* The CIF file contains the information of atomic position and unit cell information.

## Bulk
* Similar to the calculations for molecules, we need to make `POSCAR`, `POTCAR`, `INCAR`, and `KPOINTS` files.

### Making POSCAR
* Here we will see how to make a POSCAR file from a cif file using VESTA.
* VESTA is the free software for visualization and editing the molecular/bulk/surface structure: https://jp-minerals.org/vesta/jp/
* You can visualize CIF, POSCAR, xyz (and other) files with VESTA.
* To make the POSCAR file for bulk calculation, open the CIF file and then
  1. `File` -> `Export Data` and then choose `VASP file` in `File type`.
  2. Any name is OK.
  3. Either "fractional coordinate" or "Cartesian coordinate" is OK.
  4. Rename the file to `POSCAR`.

### POTCAR
* Do not forget to make POTCAR file: `cat {POTCAR_directory}/Pt/POTCAR > POTCAR`.

### INCAR
* Actually, INCAR for bulk/surface calculations can be very similar to those of molecules.
* However, we are going to use `ISMEAR = 1` instead of `ISMEAR =0` in molecular calculation. This is because Pt is metals, so use of `ISMEAR = 1` is recommended for metallic systems.
  ```
  SYSTEM = Pt bulk
  ISMEAR =  1
  NSW    =  0
  IBRION = -1
  ```

### KPOINTS
* For the bulk calculation,
```
bulk
0
Monkhorst Pack
3 3 3
0 0 0
```

### Executing VASP
* Execution of VASP is the same with the energy calculation of molecule.
  + edit `run.sh` script (or any name)
  + `qsub run.sh`

### Exercise
* Download the Au (gold) cif file from web (Materials Project etc.) and perform the bulk structure calculation.

---

## Surface
* If you want to do calculations for surfaces, the situation is a little different because the surface loses one periodic direction. Usually this axis is set as z-direction.
* To cut the periodicity of z-direction, a vacant region called **vacuum layer** in inserted in the surface calculaiton.

### Making POSCAR
* It is possible to make the surface file with VESTA; see https://qiita.com/h-nabata/items/290a575e07a2e56c7c94
* It is rather complicated. So we will use Atomic simulation environment (ASE) instead.

#### ASE
* ASE is useful Python library for theoretical/computational atomic simulations.
* To setup the Python environment, see; https://github.com/ishikawa-group/python_introduction/blob/main/setup.md.

* To make the POSCAR file for surface (e.g. fcc 111 surface), write and execute the following Python script.
  ```python{cmd}
  from ase.io import write
  from ase.build import fcc111

  surf = fcc111(a=3.92, symbol="Pt", size=[2, 2, 4], vacuum=10.0)
  write("POSCAR", surf)
  ```
* You can check the structure with VESTA, or using `ase gui POSCAR` if you've installed the ASE.
* The above script, using the ASE function `fcc111`, makes the Pt surface with the supercell size of 2x2x4, and introduces the vacuum layer of 10 Angstrom in z-direction. As noted, this vacuum layer is necessary to cut the interaction between upper and lower periodic slabs.
* More details are given in [adsorption.md](./adsorption.md).

#### Fixing some atoms
* In surface calculations, the lower part of the slab is usually fixed to mimic the bulk structure.
* To fix some atoms, make the POSCAR as following:
```
Pt
 1.000
    11.074    0.000    0.000
     0.000   11.495    0.000
     0.000    0.000   29.090
 Pt
  32
Selective dynamics
Cartesian
  0.000  2.874  0.100   F F F
  2.768  0.000  0.100   F F F
  0.000  8.621  0.100   T T T
...
```
* Put **Selective dynamics** in after the atomic number section.
* Put `F F F` at the end of atomic species, to freeze x-, y-, and z-coordinate of that atom.
* `F` means *.FALSE. to move that atom*, and `T` means *.TRUE. to move*. So open some editor, and put `F F F` at the atoms you want to freeze.

### POTCAR
* Same with bulk or molecular calculations.

### KPOINTS
* For the surface calculation, just one k-point is OK for the non-periodic direction. So when the vaccum layer is in the z-direction, following is fine.
```
surface
0
Monkhorst Pack
3 3 1
0 0 0
```
* Using more k-points gives higher accuracy, but larger computational cost.

### Executing VASP
* Same with bulk or molecular calculations.

## Exercise
1. Perform the Pt surface calculation by yourself.
2. Perform the Pt bulk calculation, and compare the Pt bulk and surface calculations by taking the **energy per Pt atom**. Which is lower (lower is more stable), bulk or surface?