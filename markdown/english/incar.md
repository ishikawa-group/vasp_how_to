# Some important INCAR tags
* This page summarizes some important INCAR tags in the VASP calculation.

## ENCUT
* Energy cutoff threshold for the plane-wave basis set.
* Plane waves satisfying following condition is included:
$$
|{\bf G} + {\bf k}| \lt G_{cut}\text{ \ where \ }E_{cut} = \frac{\hbar^2}{2m}G_{cut}^2
$$
* Lager ENCUT leads higher accuracy, but lager computational cost.
* In POTCAR, you can find ENCUT-like values `ENMAX`, `ENMIN`: These corresponds to the recommended maximum and minimum ENCUT values. Therefore, it is usually sufficient to set *ENCUT to the maximum value of ENMAXs for all the elements in the system*.

## EDIFF
* Convergence threshold (in eV unit) for the SCF calculation.
* Default value is $10^{-4}$.

## IBRION
* Tag related to the ionic motions.
    + -1: No ionic motion (single point calculation).
    + 0: Do molecular dynamics (MD).
    + 1: Do geometry optimization with RMM-DIIS algorithm.
    + 2: Do geometry optimization with conjugate gradient algorithm.
    + 3: Do geometry optimization with dampled MD algorithm.

## POTIM
* Controls the step size in geometry optimization, and also the timestep in the MD calculation (in fs).
