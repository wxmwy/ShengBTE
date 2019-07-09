include arch.make

LIBS+=-L/usr/local/cuda-10.0/lib64 -lcudart

objects:=cuda_wrapper.o cuda_processes.o conductivity.o config.o data.o input.o iterations.o misc.o gruneisen.o phonon_routines.o integrals.o processes.o scaling.o symmetry.o wedgetc.o dos_routines.o ShengBTE.o

all: ../ShengBTE
../ShengBTE: $(objects)
	$(MPIFC) $(FFLAGS) -g -o $@ $^ $(LDFLAGS) $(LIBS)
cuda_wrapper.o: cuda_wrapper.c
	icc -g -O3 -I/usr/local/cuda-10.0/include -c -o $@ $^ 
cuda_processes.o: cuda_processes.cu
	nvcc -g -O2 -gencode=arch=compute_60,code=sm_60 -c -o $@ $^
conductivity.o conductivity.mod: config.mod data.mod
config.o config.mod:  misc.mod data.mod symmetry.mod
input.o input.mod:  config.mod
iterations.o iterations.mod: data.mod config.mod wedgetc.mod
phonon_routines.o phonon_routines.mod: data.mod config.mod input.mod
dos_routines.o dos_routines.mod: config.mod
gruneisen.o gruneisen.mod: data.mod config.mod input.mod
processes.o processes.mod:  data.mod config.mod
scaling.o scaling.mod: config.mod data.mod
wedgetc.o wedgetc.mod: symmetry.mod config.mod misc.mod
ShengBTE.o: misc.mod config.mod input.mod processes.mod iterations.mod conductivity.mod scaling.mod gruneisen.mod dos_routines.mod phonon_routines.mod wedgetc.mod integrals.mod
%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../ShengBTE *.mod *.o
