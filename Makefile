CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -I ../amplitudelib_v2/ -O2
LDFLAGS = `/opt/local/bin/gsl-config --libs`

include filelist.m

CXX = g++-mp-5

all: diffraction

diffraction: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o diffraction

.cpp.o: src/subnucleon_config.hpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

tools: tools/sample_color_charges.o tools/fftwpp/fftw++.o
	$(CXX) `/opt/local/bin/gsl-config --cflags` -O2 $(LDFLAGS) tools/sample_color_charges.cpp src/dipole.o src/wilsonline.o src/ipsat_proton.o src/vector.o src/gdist_dglap.o src/subnucleon_config.o tools/fftwpp/fftw++.o -o tools/sample_color_charges -lfftw3

clean:
	rm -f $(OBJECTS)
	rm -f tools/*.o
	rm -f tools/fftwpp/*.o
	rm -f diffraction
