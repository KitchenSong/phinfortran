include arch.make

objects:=util.o qe.o config.o input.o fc.o dyn.o circular.o surface.o write_data.o device.o lesser.o velocity_mat.o phinfortran.o

all: phinfortran
phinfortran: $(objects)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
input.o input.mod: util.mod config.mod
fc.o fc.mod: input.mod qe.mod
dyn.o dyn.mod: config.mod input.mod util.mod qe.mod
circular.o circular.mod: config.mod util.mod
surface.o surface.mod: config.mod input.mod circular.mod util.mod
write_data.o write_data.mod: config.mod surface.mod input.mod
device.o device.mod: util.mod config.mod surface.mod write_data.mod input.mod
lesser.o lesser.mod: util.mod config.mod device.mod
velocity_mat.o velocity_mat.mod: input.mod util.mod config.mod fc.mod surface.mod dyn.mod write_data.mod circular.mod
phinfortran.o: qe.mod util.mod config.mod input.mod fc.mod dyn.mod circular.mod surface.mod write_data.mod device.mod lesser.mod velocity_mat.mod
%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f phinfortran *.mod *.o
