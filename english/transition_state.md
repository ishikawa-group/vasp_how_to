# Transition state
* In chemical reaction, we can define the reactant state as the starting state of some reaction and the product state as the state that the reaction is completed.
* The transition state lays in between the reactant state and product state. Its the exact position is defined as the energy-maximum-point on the minium energy path.
* Here we will learn how to calculate this transition state (TS) with VASP.
* A method called **nudged elastic band (NEB)** is often used for this purpose.
* To run the NEB, some Perl scripts provided by the Henkelman group from University of Texas is often used:
  + https://theory.cm.utexas.edu/vtsttools/

## Procedure
1. Find the structures of reactant and product states by carrying out the geometry optimization.
2. Rename the reactant and product final POSCAR files (or CONTCAR) to *POSCAR1* and *POSCAR2*, respectively. Then put these files to the current directory.
3. Make interpolated structure as the initial structures of NEB. Here, we make 4 interpolated structures.
  + `nebmake.pl POSCAR1 POSCAR2 4`
4. Find the directories 00, 01, 02, 03, 04, and 05 is made. Put the POSCAR1 and POSCAR2 files to 00 and 05, respectively.
5. Make INCAR file
```
...
IBRION  = 2
POTIM   = 0.1

IMAGES  =  4  # Change according to the image number
SPRING  = -5

ISYM    = 0   # Better to turn off the symmetry
...
```
6. Execute the VASP.

## Advanced topics
* Usually, the climbing image NEB (CINEB) is more stable than the standard NEB. To perform the CINEB, you need to download some files from the Henkelman group website above then need to compile VASP again.
* We don't mention this procedure. Interested readers should visit their website.
* Also, one of the major algorithm for the TS search is the *dimer method*. This is also included in the Henkelman's group code so need to compile VASP.