configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	mzipsat
	STATIC
	mz_ipsat/dipoleamplitude.cpp
	mz_ipsat/dglap_cpp/AlphaStrong.cpp
	mz_ipsat/dglap_cpp/EvolutionLO_nocoupling.cpp
	)

add_library(
	ddis
	STATIC
	nucleons.cpp
	ipsat_proton.cpp
	dipole.cpp
	ipglasma.cpp
	gdist_dglap.cpp
	wilsonline.cpp
	virtual_photon.cpp
	dvcs_photon.cpp
	wave_function.cpp
	qcd.cpp
	subnucleon_config.cpp
	vector.cpp
	smooth_ws_nuke.cpp
	interpolation.cpp
)

add_executable(
	subnucleondiffraction	
	main.cpp
	diffraction.cpp
	gauss_boost.cpp
	gaus_lc.cpp
	nrqcd_wf.cpp
	dis.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	)
add_executable(
	test
	test.cpp
	test_interpolation.cpp
	diffraction.cpp
	gauss_boost.cpp
	gaus_lc.cpp
	nrqcd_wf.cpp
	dis.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)	
include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	subnucleondiffraction
	PRIVATE
	ddis
	#dipolelib # LCPT dipole
	#ColorDipole  #IPSAT12 fit
	mzipsat
	GSL::gsl
	GSL::gslcblas
)
target_link_libraries(
	test
	PRIVATE
	ddis
	#dipolelib # LCPT dipole
	#ColorDipole  #IPSAT12 fit
	mzipsat
	GSL::gsl
	GSL::gslcblas
)



