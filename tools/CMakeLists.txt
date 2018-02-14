
add_executable(
	average_amplitude	
	average_amplitude.cpp
	../src/dipole.cpp
	../src/vector.cpp
	../src/wilsonline.cpp
	../src/ipglasma.cpp
	../src/ipsat_proton.cpp
	../src/gdist_dglap.cpp
	../src/subnucleon_config.cpp
	)

add_executable(
	f2fit
	f2fit.cpp
	../src/dipole.cpp
	../src/diffraction.cpp
	../src/wilsonline.cpp
	../src/ipsat_proton.cpp
	../src/vector.cpp
	../src/virtual_photon.cpp
	../src/ipglasma.cpp
	../src/wave_function.cpp
	../src/gdist_dglap.cpp
	../src/subnucleon_config.cpp
)

include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	average_amplitude
	PRIVATE
	mzipsat
	amplitude
	ColorDipole
	GSL::gsl
	GSL::gslcblas
)

target_link_libraries(
	f2fit      
	PRIVATE
	mzipsat
	amplitude
	ColorDipole
	GSL::gsl
	GSL::gslcblas
)
