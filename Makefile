CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -I ../amplitudelib_v2/ -O2
LDFLAGS = `/opt/local/bin/gsl-config --libs` libColorDipole/libraries/libColorDipole.a -lgfortran

include filelist.m

CXX = /opt/local/bin/g++-mp-5

all: diffraction

diffraction: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o diffraction

.cpp.o: src/subnucleon_config.hpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

tools: tools/sample_color_charges.o tools/fftwpp/fftw++.o tools/satscale.o tools/dipxs_intb.o tools/average_amplitude.o tools/basis_expansion.o tools/modified_ipsat.o tools/qs_npart.o
	#$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/sample_color_charges.cpp src/dipole.o src/wilsonline.o src/ipsat_proton.o src/vector.o ../amplitudelib_v2/tools/tools.o src/gdist_dglap.o src/subnucleon_config.o tools/fftwpp/fftw++.o -o tools/sample_color_charges -lfftw3
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/satscale.cpp src/dipole.o src/wilsonline.o src/ipsat_proton.o src/vector.o src/gdist_dglap.o src/ipglasma.o src/subnucleon_config.o ../amplitudelib_v2/tools/tools.o -o tools/satscale
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/dipxs_intb.cpp src/dipole.o src/wilsonline.o src/ipsat_proton.o src/vector.o src/gdist_dglap.o src/ipglasma.o src/subnucleon_config.o ../amplitudelib_v2/tools/tools.o -o tools/dipxs_intb
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/average_amplitude.o src/dipole.o src/wilsonline.o src/ipsat_proton.o src/vector.o src/gdist_dglap.o src/ipglasma.o src/subnucleon_config.o ../amplitudelib_v2/tools/tools.o -o tools/average_amplitude
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/basis_expansion.o src/dipole.o src/wilsonline.o src/ipsat_proton.o src/vector.o src/gdist_dglap.o src/ipglasma.o src/subnucleon_config.o ../amplitudelib_v2/tools/tools.o -o tools/basis_expansion
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/modified_ipsat.o  libColorDipole/libraries/libColorDipole.a -o tools/modified_ipsat
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tools/qs_npart.cpp libColorDipole/libraries/libColorDipole.a -o tools/qs_npart src/vector.o ../amplitudelib_v2/tools/tools.o



clean:
	rm -f $(OBJECTS)
	rm -f tools/*.o
	rm -f tools/fftwpp/*.o
	rm -f diffraction
