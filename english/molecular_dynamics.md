# Molecular dynamics
* In this section, we learn how to do molecular dynamics (MD) simulation with VASP.
* MD is the simulation in which *all the particles (in computational unit cell) moves according to the equation of motion*. The equation of motion is the Newton's $F = ma$ in classical MD simulation.
* In the classical MD simulation, the nuclei are the only particle move because there is no electron's degree-of-freedom in classical mechanics.
* In quantum MD simulation, mostly the nuclei are the only particle move because motion of electrons and nuclei are quite different in time scale (Born-Oppenheimer approximation).
* Therefore, in most of the quantum MD calculation, atoms move on the **potential energy surface (PES)**, and PES is calculated by solving the Schrodinger equation at each nuclei geometry.
* This procedure is very similar to the geometry optimization. In optimization, one is trying to find the minimum on the PES, while in MD one can move on the PES.

## Ensembles in MD calculation
* In MD simulation, we need to choose **conserved variables**. Possible variables are
    + $N$: number of particles (atoms, molecules)
    + $V$: volume of the simulation cell
    + $T$: temperature
    + $P$: pressure
    + $\mu$: chemical potential
* The conserved variables are the variables kept constant during the MD simulation.
* Usually, three variables are chosen for the MD calculation and that calculation is called like the *NVT ensemble*, by taking conserved variables. In this case, N, V, and T are kept constant during the MD.
* Two types of ensembles are possible in VASP: *NVE* or *NVT*. These are also called *microcanonical ensemble* or *canonical ensemble*, respectively.
* NVT simulation is more usual so we will cover it here. Keeping $N$ means we have no generation or elimination of the atoms or electsons (as usual simulations does). Keeping $V$ means cell size is fixed during the simulation. So additional factor is the temperature of the system.
* To control the temperature, we use *thermostat*. This is the heat bath connected to the system, and it supplies or extract the temperature to adjust the system's temperature to the target temperature value.
* For $NVT$ simulation, we have two options in VASP.
    1. Berendsen thermostat
    2. Nose-Hoover thermostat
* One can switch these two by changing the `SMASS` tag in INCAR.
    + `SMASS = -1` for Berendsen
    + `SMASS > 0`  for Nose-Hoover
* Usually Berendsen thermostat is easier to use but Nose-Hoover is better from theoretical viewpoint. So we learn to use the Berendsen thermostat in this page.
* In $NVT$ simulation, you need to specify the temperature by `TEBEG` and `TEEND` tags. These are the target temperature at the beggining of the simulation and end of the simulation. The units are in Kelvin (K).

## INCAR setting
* INCAR tag related to the NVT MD simulation is as follows. Other parts can be same with the previous INCAR setting.
* The full INCAR file can be found in `INCAR_MD`.
```
IBRION =  0
POTIM  =  1.0
SMASS  = -1
NSW    =  100
NBLOCK =  10
TEBEG  =  300
TEEND  =  300
```
* The MD-specific keywords are as follows:
    + `POTIM`: In MD calculation, this means the timestep in femtosecond (fs, 1.0e-15 s). Usually, 0.5-2.0 should be used.
    + `NBLOCK`: Frequency to control the temperature when using Berendsen thermostat. `NBLOCK` = 10 means temperature is scaled every 10 steps during the MD.
    + `TEBEG`: Target temperature at the beggining of the MD.
    + `TEEND`: Target temperature at the end of the MD. If this value is different from `TEBEG`, gradual increase/decrease of temperature during MD is taken.

## POSCAR
* Any POSCAR file is fine. A file with 32 water molecules are prepared; see `POSCAR_MD`.
* Copy it to `POSCAR` and do the calculation.

## KPOINTS
* Nothing special with the KPOINTS file. For above water case, 1x1x1 k-point is fine.

## POTCAR
* Prepare as usual.

## OUTPUT
* In standard output (`.out`) or OSZICAR, the simulation temperature is written (`T=300`). See how these temperatures are changing.
* The trajectory can be seen in `vasprun.xml` or `OUTCAR` file. `vasurun.xml` is lighter in size so analyze it by `ase gui vasprun.xml` (after getting .xml file to your local PC).

## Exersice
* Perform the MD calculation above.