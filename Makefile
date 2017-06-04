#########################################


FSRC = constants.f95 \
	options.f95 \
	calls.f95 \
	main.f95 \
	grid.f95 \
	matrix.f95 \
	loads.f95

FMOD = constants.f95 \
	options.f95 \
	calls.f95

##########################################


FOBJ = $(FSRC:.f95=.o)

MODS = $(FMOD:.f95=.mod)

FCMP = gfortran
FOPT = -g -Wall -Wextra -fbounds-check

LIBS = -llapack -lblas


##########################################


a.out : $(FOBJ)
	$(FCMP) $(FOPT) -o $@ $^ $(LIBS)


##########################################


%.o : %.f95
	$(FCMP) $(FOPT) -c $<


##########################################

.PHONY : clean

clean :
	rm -f a.out a.test.out *.o *.mod *.dat fort*
