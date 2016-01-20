# 
#	Evaluation of the DAs of the vector and axial currents
#       Main targets:
#            libs:  builds libColorDipole.so and libColorDipole.a
#            bins:  builds Test_Dipole, mostly for illustration

vpath %.f    ../src
vpath %.c    ../src

# include path for the newly built libraries
HEADERS=../headers/
LIBPATH=../libraries/

#default target
default:	bins

.PHONY: clean distclean
clean:
	rm -f *.o
	$(MAKE) -C dynamic clean
	$(MAKE) -C static distclean
	

distclean:	clean
		rm -f *~
		rm -f ../src/*~
		rm -f ../headers/*~
		rm -f $(LIBPATH)/*



# all targets:
all: libs TEST_DIPOLE
bins: TEST_DIPOLE

math: DAIVM_ML

#__________________________________________________________________________________________________________

# Default rule how to make object files from Fortran files
%.o: %.f
	$(FC) -c  $(FFLAGS) $< -I$(HEADERS) -o $@

# Default rule how to make object files from C files
%.o: %.c
	$(CC) -c  $(CFLAGS) $< -I$(HEADERS) -o $@

# Default rule how to make object files from C++ files
%.o: %.cpp
	$(CXX) -c  $(CFLAGS) $(CPPFLAGS) $(CXXFLAGS) $< -I$(HEADERS) -o $@
#__________________________________________________________________________________________________________


#==========================================================================================================
#                         DA--(Axial)vector distribution amplitudes
#__________________________________________________________________________________________________________

libs:
	$(MAKE) -C dynamic $(MAKEFLAGS) LIBDIPOLE FC=$(FC)
	$(MAKE) -C static $(MAKEFLAGS) LIBDIPOLE FC=$(FC)


TEST_DIPOLE: Test_Dipole.o | libs
# Static version:
	$(FC) $(FFLAGS) $^ $(LIBPATH)libColorDipole.a -o ../"Test_Dipole"
# Dynamic version, requires export LD_LIBRARY_PATH+=:`pwd`/libraries before run
# 	$(FC) $(FFLAGS) $^ -lColorDipole -L$(LIBPATH) -o ../"Test_Dipole"

#__________________________________________________________________________________________________________#


