# Report for Materials Simulation
* Choose one of the following two tasks, depending on your background.
* After finishing the task, **make one-page presentation slide (PowerPoint, KeyNote, etc.)** summarizing your computational result.
* Please send it by e-mail in any format (.ppt, .pdf, .jpg etc.)

# 1. Surface diffusion
## Background
* Diffusion of atoms or molecules on the surface is an important process for surface chemistry, catalytic chemistry, etc.
* Knowing the energetic propoerty of surface diffusion is helpful to understand this diffusion process.

## Task
* Use NEB method to find the activation energy of diffusion of oxygen atom (O) on Pt 111 surface.

## Hint
* First, perform the geometry optimizations to locate the initial and final positions.
* Adsorption sites like fcc or hcp hollow sites are favorable site for initial and final positions.
* Any computational level (cutoff, k-points, number of imagess) is OK.
* Draw and paste the whole potential energy curve obtained with the NEB method.

# 2. Density of state and adsorption
## Background
* Energetics of the adsorption of atoms or molecules on a surface depends on the electronic character of surface.
* Upon adsorption, electron transfer between surface and adsorbate would occur. This considerably affects the adsorption strength.
* Extent of charge transfer is difficult to measure with experiments, so it is sometimes a good idea to do calculation.

## Task
* Adsorb the O atom on any **two** metal surfaces of fcc111, and compare the adsorption energy of these two metals.
* Then, perform the calculation for **surface** (not surface + adsorbate), and get the `DOSCAR` file.
* Compare the DOSs of two metals, by making plots.

## Hint
* Choosing two from Ag, Pd, Rh, Ru would be fine.
* When calculating the adsorption energy, use $ E_{ads} = E_{\rm O-Surf} - \left(E_{\rm surf} + 1/2 E_{\rm O_2} \right) $. You need to calculate O2 molecule. Since O2 molecule has triplet spin state, use spin-polarized calculation (by setting `ISPIN = 2` in `INCAR`).
* Non-spin-polarized calculation is OK for surfaces and O-adsorbed surfaces.
* Plotting d-DOS instead of total DOS makes more clear difference of two metals.
