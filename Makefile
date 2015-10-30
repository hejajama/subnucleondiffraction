CXXFLAGS = `/opt/local/bin/gsl-config --cflags` -I /Users/heikki/code/amplitudelib_v2/ -O2
LDFLAGS = `/opt/local/bin/gsl-config --libs`

SOURCES = src/dipole.cpp src/gauss_boost.cpp src/wave_function.cpp src/main.cpp src/diffraction.cpp src/smooth_ws_nuke.cpp src/gdist_dglap.cpp src/ipsat_nucleons.cpp src/vector.cpp /Users/heikki/code/amplitudelib_v2/tools/tools.cpp src/wilsonline.cpp src/ipglasma.cpp
OBJECTS=$(SOURCES:.cpp=.o)

CXX = g++ #/usr/local/bin/g++
FOR = /usr/local/bin/gfortran

all: diffraction

diffraction: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o diffraction

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f diffraction
