# Default precompiler options
CPP_OPTIONS = -DHOST=\"LinuxGNU\" \
              -DMPI -DMPI_BLOCK=8000 -Duse_collective \
              -DscaLAPACK \
              -DCACHE_SIZE=4000 \
              -Davoidalloc \
              -Dvasp6 \
              -Duse_bse_te \
              -Dtbdyn \
              -Dfock_dblbuf

CPP         = gcc -E -C -w $*$(FUFFIX) >$*$(SUFFIX) $(CPP_OPTIONS)

FC          = mpif90
FCL         = mpif90

FREE        = -ffree-form -ffree-line-length-none

FFLAGS      = -w -ffpe-summary=none

OFLAG       = -O2
OFLAG_IN    = $(OFLAG)
DEBUG       = -O0

OBJECTS     = fftmpiw.o fftmpi_map.o fftw3d.o fft3dlib.o
OBJECTS_O1 += fftw3d.o fftmpi.o fftmpiw.o
OBJECTS_O2 += fft3dlib.o

# For what used to be vasp.5.lib
CPP_LIB     = $(CPP)
FC_LIB      = $(FC)
CC_LIB      = gcc
CFLAGS_LIB  = -O
FFLAGS_LIB  = -O1
FREE_LIB    = $(FREE)

OBJECTS_LIB = linpack_double.o

# For the parser library
CXX_PARS    = g++
LLIBS       = -lstdc++

##
## Customize as of this point! Of course you may change the preceding
## part of this file as well if you like, but it should rarely be
## necessary ...
##

# When compiling on the target machine itself, change this to the
# relevant target when cross-compiling for another architecture
VASP_TARGET_CPU ?= -march=native
FFLAGS     += $(VASP_TARGET_CPU)

# For gcc-10 and higher (comment out for older versions)
# FFLAGS     += -fallow-argument-mismatch
FFLAGS     += -Wno-argument-mismatch

# BLAS and LAPACK (mandatory)
# OPENBLAS_ROOT ?= /path/to/your/openblas/installation
OPENBLAS_ROOT ?= /usr/lib/x86_64-linux-gnu
BLASPACK    = -L$(OPENBLAS_ROOT) -lopenblas

# scaLAPACK (mandatory)
# SCALAPACK_ROOT ?= /path/to/your/scalapack/installation
# SCALAPACK   = -L$(SCALAPACK_ROOT) -lscalapack
SCALAPACK_ROOT ?= /usr/lib/x86_64-linux-gnu
SCALAPACK   = -L$(SCALAPACK_ROOT) -lscalapack-openmpi

LLIBS      += $(SCALAPACK) $(BLASPACK)

# FFTW (mandatory)
# FFTW_ROOT  ?= /path/to/your/fftw/installation
# INCS       += -I$(FFTW_ROOT)/include
FFTW_ROOT  ?= /usr/lib/x86_64-linux-gnu
LLIBS      += -L$(FFTW_ROOT) -lfftw3
INCS       += -I/usr/include

# HDF5-support (optional but strongly recommended)
#CPP_OPTIONS+= -DVASP_HDF5
#HDF5_ROOT  ?= /path/to/your/hdf5/installation
#LLIBS      += -L$(HDF5_ROOT)/lib -lhdf5_fortran
#INCS       += -I$(HDF5_ROOT)/include

# For the VASP-2-Wannier90 interface (optional)
#CPP_OPTIONS    += -DVASP2WANNIER90
#WANNIER90_ROOT ?= /path/to/your/wannier90/installation
#LLIBS          += -L$(WANNIER90_ROOT)/lib -lwannier
