
# ===== Compiler selection =====
FC  = gfortran
#F77 = gfortran

## ===== Flags per compiler =====
ifeq ($(FC),gfortran)
    FFLAGS = -O3 -march=native -ffast-math -funroll-loops -fopenmp
    LAPACK = -llapack -lblas
endif

ifeq ($(FC),ifort)
    FFLAGS = -O3 -xHost -qopenmp -ipo -fp-model fast=2
#   LAPACK = -L/s/fred/lapack-3.11 -llapack -L/s/fred/lapack-3.11 -lrefblas
    LAPACK = -mkl
endif

ifeq ($(FC),ifx)
#   FFLAGS = -O3 -xHost -qopenmp -ipo -fp-model fast=2
    FFLAGS = -O3 -xHost -qopenmp -fp-model fast=2
    LAPACK = -qmkl
endif

#SRC = src
SRC_ = ../hpc_fedvr_main_ref/src
SRC = src


# --- Core modules (NO main here) ---
CORE_OBJS = \
       lobatto.o \
       constants.o \
       structure_parameters.o \
       dynamic_parameters.o \
       exploit_parameters.o \
       space_time_ops.o \
       math_util.o \
       util.o \
       fedvr.o \
       fedvr_topology.o \
       fedvr_conf_struct.o \
       fedvr_derivative_ops.o \
       global_assembly.o \
       propagation.o \
       observables.o \
       io_module.o 
#     conv_tests_v2.o \

# --- Executables ---
STRUCTURE_EXE = struct
CDYN_EXE     = cdyn
QDYN_EXE     = qdyn
EXPLOIT_EXE   = exploit

all: $(STRUCTURE_EXE)

# --- Structure build ---
$(STRUCTURE_EXE): $(CORE_OBJS) main_structure.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Dynamic -----------
$(C_DYN_EXE): $(CORE_OBJS) main_dynamic.o
	$(FC) $^ -o $@ $(LAPACK)

$(DYNAMIC_EXE): $(CORE_OBJS) main_dynamic_v2.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Exploitation -----------
$(EXPLOIT_EXE): $(CORE_OBJS) main_exploit.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Compilation rule ---
%.o: $(SRC_)/%.f
	$(FC) $(FFLAGS) -c $<

%.o: $(SRC_)/%.f90
	$(FC) $(FFLAGS) -c $<

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

# --- Dependency include ---
-include *.d

.PHONY: clean
clean:
	rm -f *.o *.d *.mod fort.* \
		$(STRUCTURE_EXE) $(CDYN_EXE) $(EXPLOIT_EXE) $(QDYN_EXE)
