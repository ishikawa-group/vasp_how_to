# VASP compile
* Here, we will see how to compile the VASP package.


## Compiling VASP with Intel compilers
* See Japanese document (and translate it if necessary).


## Compiling VASP on Ubuntu with GNU compilers
* Below is how to compile the VASP on Ubuntu Linux with GNU compilers.
* You can prepare this combination (Ubuntu + GNU) on cloud resources such as Google Computing Engine (GCE).
* After compiling VASP on GCE, you can log-in (via ssh) to that virtual machine (VM) by knowning its external IP address, then you can calculate something.

### 1. Send VASP files to server
* You need to send `vasp.X.X.X.tar.gz` to the remote environment by scp after establishing ssh connection.
* Also send paw-potential files if you want.
```bash
scp vasp.6.4.3.tgz your_name@X.X.X:/home/your_name -i ~/.ssh/your_private_key
```

### 2. Install compilers
* You need to install GNU compilers (gcc, g++, gfortran, openmpi) and other tools.
```bash
sudo apt -y update
sudo apt -y install rsync
sudo apt -y install build-essential
sudo apt -y install gfortran
sudo apt -y install openmpi-bin libopenmpi-dev
```

### 3. Install libraries
* Install the mathematical libraries (FFTW3, openBLAS, SCALAPACK).
* When using apt, libraries are installed to `/usr/lib/x86_64-linux-gnu/`.
* Include files (\*.h) are placed to `/usr/include`.
* Check for the library file name (libxxx.so), then use it when modifying the makefile.
* The exact version number of scalapack (2.2 below) may change; use `apt search scalapack`.
```bash
sudo apt -y install libfftw3-bin libfftw3-dev libfftw3-mpi-dev
sudo apt -y install libopenblas-dev
sudo apt -y install libscalapack-openmpi2.2 libscalapack-openmpi-dev
```

### 4. Prepare makefile
#### copy pre-defined makefile
```bash
tar zxvf vasp.x.x.x
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
