# Oxygen Vacancy Analysis with VASP
* With VASP (or any DFT code), any atomic defects can be calculated. Here, we will focus on the oxygen vacancy of transition metal oxides, as this is widely studied case.

## Introduction
* Oxygen vacancies are point defects that occur when oxygen atoms are missing from their regular lattice positions in transition metal oxide.
* These defects play a crucial role in determining the physical, chemical, and electronic properties of materials.

## Theoretical background
### Vacancy formation energy
* The formation energy of an oxygen vacancy is a key parameter that quantifies the energy required to create a vacancy. It is calculated as:
$$
E_{\rm vac} = E_{\rm defect} - E_{\rm pristine} + n_O\cdot\mu_O
$$
where:
  - $E_{\rm defect}$: Total energy of the system with the oxygen vacancy
  - $E_{\rm pristine}$: Total energy of the prinstine (i.e., without defect) system
  - $n_O$: Number of oxygen atoms removed (typically 1)
  - $\mu_O$: Chemical potential of oxygen
* According to this definition, **a small $E_{\rm vac}$ value means oxygen vacancy can be formed easily**.
* The chemical potential of oxygen ($\mu_O$) depends on temperature and partical pressure of $O_2$ gas. For simple analysis, it can be approximated by the DFT-calculated energy of $O_2$ molecule by
$$
\mu_O = \frac{1}{2}E(O₂)
$$

## Kr&ouml;ger–Vink notation
* For defects, **Kr&ouml;ger-Vink notations** are often used to depict the atomic defects with charges. Following tables provide the most common notations.
* Notation: $\rm{A_B^C}$
  - $A$: element ($V$ for vacancy)
  - $B$: defect site ($O$ for oxygen lattice position)
  - $C$: charge ($•$ for positive, ' for negative, and $\times$ for neutral)

|           Meaning           |    Notation     |
| :-------------------------: | :-------------: |
|       Oxygen vacancy        |   $\rm{V_O}$    |
|        Metal vacancy        |   $\rm{V_M}$    |
| Charged oxygen vacancy (+1) | $\rm{V_O^{•}}$  |
| Charged oxygen vacancy (+2) | $\rm{V_O^{••}}$ |
|          electrons          |       e'        |
|            holes            |       h•        |

## Computational setup
### Procedure
* To calculate the defects, following procedures are typically taken:
1. Create a sufficiently large supercell of the pristine crystal
2. Remove one oxygen atom to create the vacancy
3. Perform structural relaxation to allow neighboring atoms to adjust their positions

### INCAR
* Following is an example INCAR for oxygen vacancy calculations:

```
ENCUT = 400        # Energy cutoff
ISMEAR = 0         # Gaussian smearing
SIGMA = 0.05       # Smearing width
LREAL = Auto       # Real space projection

# Structural relaxation
IBRION = 2         # Conjugate gradient algorithm
NSW = 100          # Maximum number of ionic steps
EDIFFG = -0.05     # Force convergence criterion (eV/Å)

LCHARG = .TRUE.    # Write charge density
```

## Analysis workflow
### 1. Pristine crystal calculation
* First, calculate the total energy and electronic structure of the pristine crystal:
1. Prepare POSCAR with the pristine crystal structure
2. Set up INCAR, KPOINTS, and POTCAR files
3. Run VASP calculation
4. Extract the total energy ($E_{\rm pristine}$)

### 2. Defect calculation
* Next, create and calculate the system with an oxygen vacancy:
1. Create a supercell from the pristine crystal
2. Remove one oxygen atom to create the vacancy
3. Adjust the INCAR file for defect calculation
4. Run VASP calculation
5. Extract the total energy ($E_{\rm defect}$)

### 3. Formation energy calculation
* Calculate the vacancy formation energy using the formula provided earlier.

### 4. Electronic structure analysis
* Analyze the impact of the vacancy on the electronic structure:
1. Density of States (DOS) analysis
  - Compare DOS of pristine and defective systems
  - Identify defect states in the band gap
2. Charge density analysis
  - Visualize the charge density difference between pristine and defective systems
  - Identify charge redistribution around the vacancy
  - Charge density analysis method is written in other file ("charge.md")
  - Depending on the vacancy type, the vacancy position has some charge density.
    + For example, when O atom is removed from MgO, two electron should be the vacancy position because O atom in MgO is actually $\ce{O^{2-}}$.
    + When removing $\ce{O^{2-}}$ from MgO, no electron is on the vacancy position. In this case, the positive charge (from $\ce{Mg^{2+}}$) tend to be on the vacancy position (thus $\rm{V_O^{••}}$).
  - To identify such charge state, charge density plot is useful.

### 5. Structural analysis
* Analyze the structural relaxation around the vacancy:
1. Measure the displacement of atoms neighboring the vacancy
2. Identify any symmetry breaking or structural distortions

## Practical considerations
### Hybrid functionals
* For more accurate electronic structure and formation energies, hybrid functionals (e.g., HSE06) can be used:
```
LHFCALC = .TRUE.
HFSCREEN = 0.2
AEXX = 0.25
```
* Using standard GGA functionals (such as PBE) is known to lead delocalized charge density. In reality, the charge density of oxygen vacancy material should localzed to the vacancy position. The hybrid functionals (or DFT + U method) is known to fix this problem.
* Note that hybrid functional calculations are computationally more demanding.

## Typical oxygen vacancy formation energies
* Here are some examples of oxygen vacancy formation energies ($E_{vac}$) in typical metal oxides calculated using DFT:
* You can see the oxygen vacancy formation energy depends on the material, and also on the calculation method.

| Material                            | E_vac (eV) | Method   | Reference |
| ----------------------------------- | ---------- | -------- | --------- |
| $\ce{TiO2}$ (rutile, (110) surface) | 2.16       | PBE      | [1]       |
| $\ce{TiO2}$ (rutile, (110) surface) | 3.66       | PBE + U  | [1]       |
| $\ce{CeO2}$ ((111) surface)         | 2.60       | PW91 + U | [2,3]     |
| $\ce{SrTiO3}$                       | 6.0        | HSE06    | [4]       |
| $\ce{LaAlO3}$                       | 8.3        | HSE06    | [4]       |

1. Morgan, B. J., & Watson, G. W., A DFT+U description of oxygen vacancies at the TiO₂ rutile (110) surface. Surface Science, 604(9-10), 715-721 (2010)
2. Ganduglia-Pirovano, M. V., et al., Oxygen vacancies in transition metal and rare earth oxides: Current state of understanding and remaining challenges. Surface Science Reports, 62(6), 219-270 (2007)
3. Nolan, M, Watson, G. W. et al., The electronic structure of oxygen vacancy defects at the low index surfaces of ceria, Surf. Sci. 595, 223 (2005)
3. Mitra, C., et al. Electronic structure of oxygen vacancies in SrTiO₃ and LaAlO₃. Physical Review B, 86(15), 155105 (2012)
