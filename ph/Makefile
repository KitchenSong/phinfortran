include arch.make

objects:=util.o config.o view_struct.o dyn.o slph.o 

all: slph
slph: $(objects)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)
config.o config.mod: util.mod
view_struct.o view_struct.mod: util.mod config.mod
dyn.o dyn.mod: util.mod config.mod view_struct.mod
slph.o: util.mod config.mod view_struct.mod
%.o %.mod: %.f90
	$(MPIFC) $(FFLAGS) -c -o $*.o $<
	touch $*.mod
clean:
	rm -f slph *.mod *.o
