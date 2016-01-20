# 
#       Main targets:
#            libs:  builds dynamic libColorDipole.so

vpath %.f    ../../src
vpath %.c    ../../src

HEADERS=../../headers/
LIBPATH=../../libraries/


.PHONY: clean distclean
clean:
	rm -f *.o

distclean:	clean
		rm -f *~

#default target
default:	libs

# demos
w:
# all targets:
libs: LIBDIPOLE

#__________________________________________________________________________________________________________

# Default rule how to make object files from Fortran files
%.o: %.f
	$(FC) -c -fPIC -I../../inc $(FFLAGS) $< -o $@

# Default rule how to make object files from C files
%.o: %.c
	$(CC) -c  -fPIC $(CFLAGS) $< -I$(HEADERS) -o $@

# Default rule how to make object files from C++ files
%.o: %.cpp
	$(CXX) -c  -fPIC $(CFLAGS) $(CPPFLAGS) $(CXXFLAGS) $< -I$(HEADERS) -o $@
#__________________________________________________________________________________________________________



LIBDIPOLE: bmats.o debug.o gluonEvolve.o mbits.o mchar.o obsolete.o ome.o pij_lo.o pij_nlo.o qcdalf.o qcdcvol.o qcdevol.o qcdfast.o qcdgrd.o qcdini.o qcdpdf.o qcdspl.o qcdstore.o qcdthrs.o qcdutil.o qcdwfun.o qcdwgtq.o qcdwgtu.o qpolint.o sigma_dipole.o store.o usrcvol.o usrerr.o usrevol.o usrfast.o usrgrd.o usrini.o usrpdf.o usrwgt.o utils.o vogelsang.o wgplg.o xa2hgp.o xpij2p.o xpns2p.o sigmaDipole.o
	$(FC) $(FFLAGS) -shared $^ -L$(LIBPATH) -o $(LIBPATH)"libColorDipole.so"






 
