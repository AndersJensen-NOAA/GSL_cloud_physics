physics_path = ../..
cldmp_path = ..
VPATH = ${physics_path}:${cldmp_path}

FC      = gfortran
FCFLAGS = -fdefault-real-8 -fdefault-double-8
SRC     = machine.F90,${cldmp_path}/module_mp_thompson.F90,${physics_path}/module_mp_radar.F90
OBJS    = machine.o,module_mp_thompson.o,module_mp_radar.o

clean:
	rm -f *.o *.mod driver

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^

driver.o: module_mp_thompson.o machine.o
module_mp_thompson.o: machine.o module_mp_radar.o
module_mp_radar.o: machine.o
