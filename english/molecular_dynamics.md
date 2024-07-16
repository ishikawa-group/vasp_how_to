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
    + $H$: enthalpy
* The conserved variables are the variables kept constant during the MD simulation.
* Usually, three variables are chosen for the MD calculation, and the MD simulation is labeled by these conserved variables.
* Four types of ensembles are possible in VASP: *NVE*, *NVT*, *NpT*, and *NpH*. Among them, formar two ensembles are widely used, and also called *microcanonical ensemble* and *canonical ensemble*, respectively.
* NVT simulation is more usual so we will cover it here. Keeping $N$ means we have no generation or elimination of the atoms or electsons (as usual simulations does). Keeping $V$ means cell size is fixed during the simulation. So additional factor is the temperature of the system.
* To control the temperature, we use *thermostat*. This is the heat bath connected to the system, and it supplies or extract the temperature to adjust the system's temperature to the target temperature value.
* For $NVT$ simulation, we have two options in VASP.
    1. Berendsen thermostat
    2. Andersen thermostat
    3. Nose-Hoover thermostat
    4. Langevion thermostat
    5. Multiple Andersen thermostat
* Among them, Berendesen and Nose-Hoover is popular so we will cover them here. How to set these thermostats with INCAR will be introduced later.
* Usually Berendsen thermostat is easier to use but Nose-Hoover is better from theoretical viewpoint. So we learn to use the Berendsen thermostat in this page.
* The full details of the MD calculation in VASP is shown in https://www.vasp.at/wiki/index.php/Molecular_dynamics_calculations

## INCAR setting
* The MD-specific keywords are as follows:
1. Common to all MD calculation
    + Set `IBRION` to 0.
    + `NSW`: The number of steps (same with geometry optimization).
    + `POTIM`: In MD calculation, this means the timestep in femtosecond (fs, $10^{-15}$ s). Usually, 0.5-2.0 should be used.
2. Common to all NVT MD calculation
    + `TEBEG`: Target temperature (in Kelvin) at the beggining of the MD.
    + `TEEND`: Target temperature at the end of the MD. If this value is different from `TEBEG`, gradual increase/decrease of temperature during MD is taken.
3. When using Berendsen thermostat
    + Set `SMASS` to -1.
    + `NBLOCK`: Frequency of the velocity scaling. Setting this value to, for example 10, means temperature is scaled every 10 steps.
4. When using Nose-Hoover thermostat
    + Set `SMASS` to some positive value (0.5~2.0 is standard).
    + `MDALGO = 2`.

* So the INCAR part for the Berendsen NVT MD becomes like
```
IBRION =  0
POTIM  =  1.0
SMASS  = -1
NSW    =  100
NBLOCK =  10
TEBEG  =  300
TEEND  =  300
```
* In Nose-Hoover, temperature fluctuates around the target value. Adjusting the `SMASS` value may resolve this to some extent.
* `SMASS` in Nose-Hoover controls the frequency of the coupling to the heat bath, so it is a parameter we need to fix. Usually, larger `SMASS` leads slow increase/decrease of temperature.

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