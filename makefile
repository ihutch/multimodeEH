#########################################################################
# This makefile should build on linux provided it has access to the X11 library
# which typically requires the development package, a fortran compiler,
# and git. It won't work on MSWindows.
# Decide the FORTRAN compiler and create the accis graphics routines:
include ACCIS.mk
#########################################################################
ifneq ("$(FORTRAN)","mpif90")
BLAH:=\
$(shell echo "WARNING Compiler $(FORTRAN) probably not MPI-capable.">&2;\
echo "May need to substitute mpiloops->mpiloopsdummy to compile." >&2;)
endif
#########################################################################
LIBRARIES := $(LIBRARIES) -lmodbess
LIBDEPS := $(LIBDEPS) libmodbess.a
COMPILE-SWITCHES:=$(COMPILE-SWITCHES) -Wno-unused-dummy-argument
#########################################################################
MODULES=shiftgen.o iterfind.o mpiloops.o
#########################################################################
# Patterns for compilation etc.
%.o : %.f makefile ;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f

%.o : %.f90 makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f90

%.o : %.F makefile;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.F

% : %.f $(ACCISX) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f $(LIBPATH) $(LIBRARIES)

% : %.f90  makefile $(ACCISX) $(MODULES) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f90 $(MODULES) $(LIBPATH) $(LIBRARIES)

# Defeat the Modula-2 make booby-trap.
% : %.mod


% : %.F $(ACCISX) $(LIBDEPS);
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBPATH) $(LIBRARIES)
#########################################################################
# A specific library of modified Bessel functions.
libmodbess.a : bessmodIs.f
	$(FORTRAN) -c $(COMPILE-SWITCHES) bessmodIs.f
	ar -crs libmodbess.a bessmodIs.o

altlibmodbess : toms715.f90
	$(FORTRAN) -c $(COMPILE-SWITCHES) toms715.f90
	ar -crs libmodbess.a toms715.o

mpiloopsdummy :
	ln -f -s mpiloopsdummy.f90 mpiloops.f90
	touch mpiloops.f90
	make modules
	make osolvomk

mpiloopstrue :
	ln -f -s mpiloopstrue.f90 mpiloops.f90
	touch mpiloops.f90
	make modules
	make osolvomk

modules : $(MODULES)

clean :
	rm -f *.o *.mod  bessmodsums omsolve fhgfuncmain libmodbess.a plot000*.ps testshiftgen slowstabgen osolvomk



