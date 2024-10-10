## Why quantum mechanics?
* The behavior of the particles is governed by the equation of motion, and its classical mechanical version is known as Newton's law.
* The proper description of atoms, molecules, and electrons is given by the laws of quantum mechanics. For this reason, we need to consider the **Schrödinger equation**, which is a quantum-mechanical equation of motion.
* If the solutions of the Schrödinger equations are generated without reference to experimental data, the methods are usually called "ab initio" (Latin: "from the beginning") or "first principle".

## The Schrödinger equation
* In solid-state physics, the fundamental interaction we are interested in is the electrostatic interaction.
* Here we introduce three assumptions to the SE
    1. time-independent
    2. non-relativistic
    3. Born-Oppenheimer approximation
* Under these approximations, the system of nuclei and electrons is described with a Hamiltonian below
```math
\hat{H} = T_{\rm nuc} + T_{\rm el} + V_{\rm nuc-nuc} + V_{\rm nuc-el} + V_{\rm el-el}
```

# Hamiltonian
* The terms in Hamiltonian are written in the atomic unit as
```math
\begin{align*}
T_{\rm nuc} &= \sum_{I=1}^L\frac{\nabla_I^2}{2M_I} &\text{(kinetic energy of nuclei)} \\
T_{\rm el}  &= \sum_{i=1}^N\frac{\nabla_i^2}{2} &\text{(kinetic energy of electrons)} \\
V_{\rm nuc-nuc} &= \frac{1}{2}\sum_{I \ne J}\frac{Z_I Z_J}{|{\bf R}_I - {\bf R}_J|} &\text{(nuclei-nuclei repulsion)} \\
V_{\rm nuc-el} &= -\sum_{i,I}\frac{Z_I}{|{\bf r}_i - {\bf R}_I|} &\text{(nuclei-electron attraction)} \\
V_{\rm el-el} &= \frac{1}{2}\sum_{i \ne j}\frac{1}{|{\bf r}_i - {\bf r}_j|} &\text{(electron-electron repulsion)}
\end{align*}
```

## Dirac's braket notation
* It is convenient to use Dirac's "bra-ket notation" for wave functions and multi-dimensional integrals in electronic structure theory to simplify the notation. The equivalences are defined as
```math
\begin{align*}
\ket{\Psi} \equiv \Psi, &\hspace{5pt} \bra{\Psi} \equiv \Psi^{*} \\
\int{d{\bf r} \Psi^{*}\Psi} &= \braket{\Psi|\Psi} \\
\int{d{\bf r} \Psi^{*}\hat{H}\Psi } &= \braket{\Psi|\hat{H}|\Psi}
\end{align*}
```
* The ket $\ket{\Psi}$ denotes a wave function while the bra $\bra{\Psi}$ denotes a complex conjugate wave function $\Psi^{*}$. The combined bracket denotes that the whole expression should be integrated over all coordinates.

# Effective Hamiltonian
* Except for the simplest cases, there are no simple ways to solve the SE in a closed analytical form so we have to solve it numerically.
* The SE is a second-order partial differential equation (PDE), so in principle, it can be directly solved. However, this needs integration over a large number of dimensions ($3\times N_{\rm elec}$), which is impossible.
* This difficulty can be solved by two approaches
    1. Approximate the electron-electron interaction by the effective one-electron problem. This reduces the $3N$ dim. integration to a sum of 3 dim. integrations.
    2. Expanding the wave function in some suitable basis set, because this will convert the PDE into a set of algebraic equations.
* The methodology based on the above approaches is to use the **effective hamiltonian**.
```math
\left[\hat{h} + v_{\rm eff}({\bf r})\right]\psi_i({\bf r}) = \varepsilon_i \psi_i({\bf r})
```

## Hartree product
* We are mainly interested in the electronic ground state energy $E_0$.
* There is an important quantum mechanical principle - the *Rayleigh-Ritz variational principle* - that provides a route to find approximate solutions for $E_0$.
* It states that the expectation value of $\hat{H}$ of any $\Psi$ is always higher than or equal to the exact $E_0$, i.e.
```math
E_0 \le \frac{\braket{\Psi|\hat{H}|\Psi}}{\braket{\Psi|\Psi}}
```
* So, the expectation value calculated by the wave function at your hand $\Psi$  will always be an upper bound for the true ground state energy. By improving $\ket{\Psi}$, you will have a lower expectation value that is closer to the true ground state energy.

##
* Since $V_{\rm nuc-el}$ is an effective external potential for an electron, we write it as a function of $r$ and also define the one-electron Hamiltonian $\hat{h}$ as
```math
\hat{v}_{\rm ext}({\bf r}) = -\sum_I \frac{Z_I}{|{\bf r} - {\bf R}_I|} \\
\hspace{2mm} \hat{h}({\bf r}) = -\frac{\nabla^2}{2} + \hat{v}_{\rm ext}({\bf r})
```
* Then one can form the one-electron SE as
```math
\hat{h}({\bf r})\psi_i({\bf r}) = \epsilon_i\psi_i({\bf r})
```
* In this case, the $N$-electron wave function can be expressed by the product of $\psi_i$ as
```math
\Psi_{\rm HP}({\bf r}_1, {\bf r}_2, \cdots, {\bf r}_N) = \psi_1({\bf r}_1)\psi_2({\bf r}_2)\cdots\psi_N({\bf r}_N)
```
* This wave function is called the **Hartree product**, and it is the first crude guess for the true $N$-electron wave function.
* Note that $\psi_i$ is orthonormal thus $\braket{\psi_i({\bf r})|\psi_j({\bf r})} = \delta_{ij}$.

## Spatial and spin orbitals
* The one-electron wave function $\psi_i({\bf r})$ is called the *orbital*, and the Hartree product means that $N$-electron wave function is expressed by the product of orbitals.
* Up to now, we assumed that the orbital depends only on ${\bf r}$, but an electron has the spin degree of freedom. We write this spin variable by $\omega$, and combine it with the spatial coordinate ${\bf r}$ as ${\bf x} = ({\bf r}, \omega)$.

##
* Let the one-electron wave function in ${\bf x}$ as $\chi_i({\bf x})$.
* Assuming that ${\bf r}$ and $\omega$ are independent, we have $\chi_i({\bf x}) = \psi_i({\bf r})\sigma_i(\omega)$, where $\psi$ and $\sigma$ denote the spatial and spin parts.
* $\chi$, $\phi$, $\sigma$ are a spin-orbital, spatial orbital, and spin function.
* Since an electron has no chance to take both $\alpha$ and $\beta$ spin simultaneously, following integration over the spin variable holds.
```math
\begin{align*}
\int d\omega \alpha^{*}(\omega)\alpha(\omega) = 1 \\
\int d\omega \beta^{*}(\omega)\beta(\omega) = 1 \\
\int d\omega \alpha^{*}(\omega)\beta(\omega) = 0 \\
\int d\omega \beta^{*}(\omega)\alpha(\omega) = 0
\end{align*}
```

## Hartree equation
* Using the spin orbital $\chi$ above, we determine the expectation value of the Hamiltonian
```math
\hat{H} = \hat{h}({\bf r}) + \frac{1}{2}\sum_{i\ne j}^N\frac{1}{|{\bf r}_i - {\bf r}_j|}
```
with respect to the Hartree product. This becomes
```math
\begin{align*}
\braket{\Psi_{\rm HP}|\hat{H}|\Psi_{\rm HP}} = & \sum_{i=1}^N\int d{\bf x} \chi_i^{*}({\bf x}) \hat{h({\bf r})} \chi_i({\bf x}) \\
&+ \frac{1}{2}\sum_{i=1}^N\int d{\bf x} d{\bf x}' \chi_i^{*}({\bf x})\chi_j^{*}({\bf x}') \frac{1}{|{\bf r} - {\bf r}'|} \chi_i({\bf x})\chi_j({\bf x}') \\
&+ V_{\rm nuc-nuc}
\end{align*}
```

##
* Now we minimize this w.r.t. $\chi_i({\bf x})$ under the constraint that $\chi_i^{*}({\bf x})$ is normalized.
* This is a typical variational problem with the constraint taken into account via *Lagrange multipliers*, which gives
```math
\frac{\delta}{\delta\chi_i^{*}}\left[\braket{\Psi_{\rm HP}|\hat{H}|\Psi_{\rm HP}}-\sum_{i=1}^N\left\{\epsilon_i\left(1-\braket{\chi_i|\chi_i}\right)\right\}\right] = 0
```

##
* The $\epsilon_i$ act as Lagrange multipliers ensuring the normalization of $\chi_i({\bf x})$. This leads to the so-called **Hartree equation** as
```math
\left[\hat{h} + \sum_{j=1}^N\int d{\bf x}'\chi_j^{*}({\bf x}')\chi_j({\bf x}')\frac{1}{|{\bf r}-{\bf r}'|}\right]\chi_i({\bf x}) = \epsilon_i\chi_i({\bf x})
```
* This shows that an effective one-electron SE is solved for an electron embedded in the electrostatic field of all electrons (including itself).

## Hartree potential
* Using the electron density
```math
\rho({\bf r}) = \sum_{i=1}^N \psi_i^{*}({\bf r})\psi_i({\bf r})
```
the *Hartree potential* $\hat{v}_H$ can be defined as
```math
\hat{v}_H({\bf r}) = \int d{\bf r}' \frac{\rho({\bf r}')}{|{\bf r} - {\bf r}'|}
```
which corresponds to the electrostatic potential of all electrons. With $\hat{v}_H$, the Hartree equation can be written as
```math
\left[\hat{h} + \hat{v}_H \right]\chi_i({\bf x}) = \epsilon_i\chi_i({\bf x})
```

## Self-consistent field
* The Hartree equation has the form of one-electron SE. However, the solutions $\chi_i({\bf x})$ enter the effective one-particle Hamiltonian via $\hat{v}_H$.
* This dilemma can be resolved by using an iterative algorithm: One starts with some initial guess for the wave functions that form the effective one-particle Hamiltonian. The Hartree equations are then solved and a new set of solutions are determined.
* This cycle is repeated so often until the iterations no longer modify the solutions, i.e. self-consistency is reached. Such a method is known as the **self-consistent field (SCF)** method.

## Kohn-Sham equation
* In physics and quantum chemistry, specifically density functional theory, the Kohn–Sham equation is the non-interacting Schrödinger equation (more clearly, Schrödinger-like equation) of a fictitious system (the "Kohn–Sham system") of non-interacting particles (typically electrons) that generate the same density as any given system of interacting particles.
* In the Kohn–Sham theory the introduction of the noninteracting kinetic energy functional Ts into the energy expression leads, upon functional differentiation, to a collection of one-particle equations whose solutions are the Kohn–Sham orbitals.
* The Kohn–Sham equation is defined by a local effective (fictitious) external potential in which the non-interacting particles move, typically denoted as vs(r) or veff(r), called the Kohn–Sham potential. If the particles in the Kohn–Sham system are non-interacting fermions (non-fermion Density Functional Theory has been researched), the Kohn–Sham wavefunction is a single Slater determinant constructed from a set of orbitals that are the lowest-energy solutions to
```math
\left[\hat{h} + v_{\rm KS}({\bf r})\right]\psi_i({\bf r}) = \varepsilon_i \psi_i({\bf r})
```
* The Kohn–Sham equations are found by varying the total energy expression with respect to a set of orbitals, subject to constraints on those orbitals, to yield the Kohn–Sham potential as
```math
v_{\rm eff}({\bf r}) = v_{\rm ext}({\bf r}) + \int\frac{\rho({\bf r'})}{|{\bf r}-{\bf r'}|}d{\bf r'} + \frac{\delta E_{\rm XC}[\rho]}{\delta \rho({\bf r})}
```
* $E_{\rm XC}[\rho]$ in the last term is called the **exchange–correlation functional**.

# Basis set expansion
## Localized basis set (Gaussian)
## Plane wave basis set

# Periodic systems
## Reciprocal lattice
* The reciprocal lattice is a mathematical object used in solid-state physics to describe the periodicity of a crystal structure in momentum space rather than real space.
* The reciprocal lattice is used instead of the real lattice due to its advantages in analyzing wave-like phenomena in solids.
* The reciprocal lattice represents the periodicity of a crystal in momentum space, which is particularly relevant when dealing with wave-like behaviors, such as electrons in solids.
* We consider systems made up of a unit cell periodically repeated throughout all space by rigid translation along the lattice vectors ${\bf a_1}, {\bf a_2}$ and ${\bf a_3}$.
* The corresponding reciprocal lattice is spanned by the vectors
```math
\begin{align*}
{\bf b_1} = \frac{2\pi}{\Omega}{\bf a_2}\times{\bf a_3} \\
{\bf b_2} = \frac{2\pi}{\Omega}{\bf a_3}\times{\bf a_1} \\
{\bf b_3} = \frac{2\pi}{\Omega}{\bf a_1}\times{\bf a_2} \\
\end{align*}
```
where $\Omega={\bf a_1}\cdot {\bf a_2}\times{\bf a_3}$ is the volume of the unit cell.

## Bravais lattice
* A Bravais lattice is a discrete set of points in three-dimensional space that represents the periodic arrangement of atoms in a solid material.
* It is defined by a set of points generated by translating a set of basis vectors in three-dimensions.
* There are 14 types of Bravais lattices categolized into seven crystal systems:
* Each Bravais lattice type is characterized by its symmetry and the angles and lengths of the lattice vectors.

## Bloch theorem
* 
## Periodic Schrodinger equation
## Periodic Kohn-Sham equation
