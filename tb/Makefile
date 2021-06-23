include arch.make

objects:=util.o config.o input.o param.o hamiltonian.o circular.o surface.o device.o velocity_mat.o write_data.o infortran.o

all: infortran
infortran: $(objects)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
input.o input.mod: config.mod util.mod
param.o param.mod: input.mod
hamiltonian.o hamiltonian.mod: config.mod input.mod param.mod util.mod 
circular.o circular.mod: config.mod util.mod
surface.o surface.mod: config.mod input.mod param.mod util.mod circular.mod
device.o device.mod: util.mod config.mod surface.mod write_data.mod input.mod param.mod
write_data.o write_data.mod: config.mod surface.mod input.mod
velocity_mat.o velocity_mat.mod: input.mod util.mod config.mod surface.mod param.mod hamiltonian.mod write_data.mod circular.mod
infortran.o: config.mod input.mod param.mod hamiltonian.mod surface.mod device.mod write_data.mod circular.mod util.mod
%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f ../infortran *.mod *.o
