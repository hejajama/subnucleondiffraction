configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../src/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

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
	amplitude_bline
	amplitude_bline.cpp
	../src/dipole.cpp
	../src/vector.cpp
	../src/wilsonline.cpp
	../src/ipglasma.cpp
	../src/ipsat_proton.cpp
	../src/gdist_dglap.cpp
	../src/subnucleon_config.cpp
	)

add_executable(
	b_r_angle 	
	b_r_angle.cpp
	../src/dipole.cpp
	../src/vector.cpp
	../src/wilsonline.cpp
	../src/ipglasma.cpp
	../src/ipsat_proton.cpp
	../src/gdist_dglap.cpp
	../src/subnucleon_config.cpp
	)

add_executable(
	b_r_grid  	
	tabulate_b_r_angle.cpp
	../src/dipole.cpp
	../src/vector.cpp
	../src/wilsonline.cpp
	../src/ipglasma.cpp
	../src/ipsat_proton.cpp
	../src/gdist_dglap.cpp
	../src/subnucleon_config.cpp
	)

#add_executable(
#	f2fit
#	f2fit.cpp
#	../src/dipole.cpp
#	../src/diffraction.cpp
#	../src/wilsonline.cpp
#	../src/ipsat_proton.cpp
#	../src/vector.cpp
#	#../src/virtual_photon.cpp
#	../src/ipglasma.cpp
##	../src/wave_function.cpp
##	../src/gdist_dglap.cpp
#	../src/subnucleon_config.cpp
#	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
#)

include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	average_amplitude
	PRIVATE
	mzipsat
	amplitude
	#dipolelib
	#ColorDipole
	GSL::gsl
	GSL::gslcblas
)
target_link_libraries(
	amplitude_bline
	PRIVATE
	mzipsat
	amplitude
	#dipolelib
	#ColorDipole
	GSL::gsl
	GSL::gslcblas
)
target_link_libraries(
	b_r_angle
	PRIVATE
	mzipsat
	amplitude
	#dipolelib
	#ColorDipole
	GSL::gsl
	GSL::gslcblas
)
target_link_libraries(
	b_r_grid 
	PRIVATE
	mzipsat
	#dipolelib
	amplitude
	#ColorDipole
	GSL::gsl
	GSL::gslcblas
)

#target_link_libraries(
#	f2fit      
#	PRIVATE
#	mzipsat
#	amplitude
#	dipolelib
#	#ColorDipole
#	GSL::gsl
#	GSL::gslcblas
#)
