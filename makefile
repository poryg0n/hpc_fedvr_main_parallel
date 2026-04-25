FC = gfortran
#FFLAGS = -O3 -Wall -cpp -MMD -MP -fcheck=all -g
FFLAGS = -O3 -march=native -ffast-math -funroll-loops -fopenmp
LAPACK = -llapack -lblas
#LAPACK = -L/s/fred/lapack-3.11 -llapack -L/s/fred/lapack-3.11 -lrefblas
SRC = src

# --- Core modules (NO main here) ---
CORE_OBJS = \
       lobatto.o \
       constants.o \
       structure_parameters.o \
       dynamic_parameters.o \
       exploit_parameters.o \
       math_util.o \
       util.o \
       fedvr.o \
       fedvr_topology.o \
       fedvr_conf_struct.o \
       fedvr_derivative_ops.o \
       space_time_ops_v2.o \
       global_assembly.o \
       propagation_v2.o \
       observables_v2.o \
       conv_tests_v2.o \
       io_module_v2.o 

# --- Executables ---
STRUCTURE_EXE = structure
DYNAMIC_EXE  = dynamic
EXPLOIT_EXE  = exploit

all: $(STRUCTURE_EXE)

# --- Structure build ---
$(STRUCTURE_EXE): $(CORE_OBJS) main_structure_v2.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Future ---
$(DYNAMIC_EXE): $(CORE_OBJS) main_dynamic_v2.o
	$(FC) $^ -o $@ $(LAPACK)

$(EXPLOIT_EXE): $(CORE_OBJS) main_exploit_v2.o
	$(FC) $^ -o $@ $(LAPACK)

# --- Compilation rule ---
%.o: $(SRC)/%.f
	$(FC) $(FFLAGS) -c $<

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

# --- Dependency include ---
-include *.d

.PHONY: clean
clean:
	rm -f *.o *.d *.mod fort.* $(STRUCTURE_EXE) $(DYNAMIC_EXE) $(EXPLOIT_EXE)

