configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_executable(
	subnucleondiffraction	
	dipole.cpp
	gauss_boost.cpp
	wave_function.cpp
	main.cpp
	diffraction.cpp
	gdist_dglap.cpp
	vector.cpp
	wilsonline.cpp
	ipglasma.cpp
	subnucleon_config.cpp
	ipsat_proton.cpp
	nucleons.cpp
	virtual_photon.cpp
	dis.cpp
	smooth_ws_nuke.cpp
	gaus_lc.cpp
	mz_ipsat/dipoleamplitude.cpp
	mz_ipsat/dglap_cpp/AlphaStrong.cpp
	mz_ipsat/dglap_cpp/EvolutionLO.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
	)

include_directories(.)

find_package(GSL REQUIRED)    #
target_link_libraries(
	subnucleondiffraction
	PRIVATE
	amplitude
	ColorDipole
	GSL::gsl
	GSL::gslcblas
)


