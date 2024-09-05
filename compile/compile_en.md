# VASP compile
* Here, we will see how to compile the VASP package.

## Compling VASP on Ubuntu with GNU compilers
* Below is how to compile the VASP on Ubuntu Linux with GNU compilers.

### 1. Send VASP files to server
* You need to send `vasp.X.X.X.tar.gz` to the remote environment by scp after establishing ssh connection.
```bash
scp vasp.6.4.3.tgz your_name@X.X.X:/home/your_name -i ~/.ssh/your_private_key
```

### 2. Install compilers
* You need to install GNU compilers (gcc, g++, gfortran) and MPI compilers.
```bash
sudo apt install build-essential
sudo apt install gfortran
sudo apt install openmpi-bin libopenmpi-dev
```

### 3. Install libraries
* Install the mathematical libraries (FFTW3, openBLAS, SCALAPACK).
* When using apt, libraries are installed to `/usr/lib/x86_64-linux-gnu/`.
* Include files (\*.h) are placed to `/usr/include`.
* Check for the library file name (libxxx.so), then use it when modifying the makefile.
```bash
sudo apt install libfftw3-3 libfftw3-dev libfftw3-doc
sudo apt install libopenblas-dev
sudo apt install libscalapack-openmpi2.1 libscalapack-openmpi-dev
```

### 4. Prepare makefile
#### copy pre-defined makefile
```bash
cd vasp.x.x.x
cp ./arch/makefile.include.gnu ./makefile.include
```

#### Edit makefile.include
* You need to modify `makefile.include` to change FFTW, openBLAS, and SCALAPACK parts.
* An example is put in this directory: `makefile.include.gnu`.

### 5. make
```bash
make veryclean
make all
```
