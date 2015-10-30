CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -I /Users/heikki/code/amplitudelib_v2/ -O2
LDFLAGS = `/opt/local/bin/gsl-config --libs`

include filelist.m

CXX = g++ #/usr/local/bin/g++
FOR = /usr/local/bin/gfortran

all: diffraction

diffraction: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o diffraction

.cpp.o: src/subnucleon_config.hpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f diffraction
