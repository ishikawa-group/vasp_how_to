# Density of state (DOS)
* DOS represents how densely the orbitals are packed in each energy range.
* A larger value of DOS means that the states are more densely packed in that energy range.

## Input file
* No specifical tag is necessary to output the total DOS, but `LORBIT` need to be set when PDOS is needed (mentioned later).
* Note that the correct DOS cannot be obtained unless the number of k-points is sufficiently large.
* The smearing parameters in the `INCAR` (ISMEAR, SIGMA) also affect the accuracy of the DOS.

## After VASP calculation
* A file called `DOSCAR` will be output.
* Be aware that the output differs between spin-polarized and non-polarized calculations.
* For non-polarized calculations, the meaning of each column in `DOSCAR` is as follows:
    `energy    DOS    integrated-DOS`

## Making plottable DOS file
* `DOSCAR` has many sections, so it is not easy to make plot from it.
* It is useful to use `vaspkit` for this purpose.

### Installing vaspkit
1. Go to vaspkit website: https://vaspkit.com/index.html
2. Go to the "latest release page" (in SourceForge).
3. Find the latest tar.gz file for Linux, and get the URL.
4. Go to your VASP-installed-computer (supercomputer).
5. Download the tar.gz file by using `wget` (`wget URL`).
6. Extract with `tar zxvf vaspkit-xxx.tar.gz`.
7. `cd vaspkit.x.x.x`
8. `source setup.sh`
9. `source ~/.bashrc`

### Using vaspkit
* vaspkit is quite easy to use. Just type `vaspkit` in terminal, and follow the instruction.
* Total DOS
  + To generate the DOS-related functions, type `11`.
  + To get the total DOS, type `111`.
  + Total DOS file `TDOS.dat` is generated. The energy position is shifted so as the Fermi energy becomes 0.

## Plot
* Copy `*.dat` file generated from vaspkit. Now it is easy to make plot.
* Any software (Excel, Gnuplot, etc) is OK.
* When using gnuplot, the command becomes like 
```bash
gnuplot
> plot "TDOS.dat" using 1:2 with lines
> plot "TDOS.dat" using 1:2 with lines, "TDOS.dat" using 1:3 with lines # when spin-polarized case
```

## PDOS
* The projected DOS (PDOS) refers to the DOS decomposed by angular momentum (s, p, d, etc.).
* To get the PDOS, use `LORBIT` tag in `INCAR`.
    + `LORBIT = 10`: makes s-, p-, d-decomposed DOSs.
    + `LORBIT = 11`: makes s-, px-, py-, pz-, dxy- ... DOSs.
* Plotting procedure is same with the total DOS.

## LDOS
* The localized DOS (LDOS) refers to the DOS calculated separately for each atom.
* LDOS can be obtained from `vaspkit`.