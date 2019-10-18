export FFLAGS=-traceback -debug -O2 -static_intel -xhost -qopt-report=5 -qopt-report-phase:vec -qopenmp
#export FFLAGS=-traceback -debug -O2 -static_intel -qopenmp
export LDFLAGS=-L/usr/local/spglib/lib -lsymspg -qopenmp
export MPIFC=mpiifort
export MKLROOT=/opt/intel/mkl
MKL=$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group	\
$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a				\
 $(MKLROOT)/lib/intel64/libmkl_sequential.a				\
 $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
export LAPACK=$(MKL)
export LIBS=$(LAPACK)
