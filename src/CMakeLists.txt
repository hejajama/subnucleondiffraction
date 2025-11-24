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
find_package(OpenMP)

# Build bundled Cuba library from vendored sources when requested
set(CUBA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Cuba-4.2.2")
set(CUBA_SRC_DIR "${CUBA_ROOT}/src")

# Build each algorithm in its own object library with precise include order
add_library(cuba_vegas_objs OBJECT ${CUBA_SRC_DIR}/vegas/Vegas.c)
target_include_directories(cuba_vegas_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/vegas"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_vegas_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_vegas_objs PROPERTIES C_STANDARD 99)

add_library(cuba_suave_objs OBJECT ${CUBA_SRC_DIR}/suave/Suave.c)
target_include_directories(cuba_suave_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/suave"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_suave_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_suave_objs PROPERTIES C_STANDARD 99)

add_library(cuba_divonne_objs OBJECT ${CUBA_SRC_DIR}/divonne/Divonne.c)
target_include_directories(cuba_divonne_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/divonne"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_divonne_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_divonne_objs PROPERTIES C_STANDARD 99)

add_library(cuba_cuhre_objs OBJECT ${CUBA_SRC_DIR}/cuhre/Cuhre.c)
target_include_directories(cuba_cuhre_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/cuhre"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_cuhre_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_cuhre_objs PROPERTIES C_STANDARD 99)

# Bundle into a single static library
add_library(cuba STATIC
    $<TARGET_OBJECTS:cuba_vegas_objs>
    $<TARGET_OBJECTS:cuba_suave_objs>
    $<TARGET_OBJECTS:cuba_divonne_objs>
    $<TARGET_OBJECTS:cuba_cuhre_objs>
    ${CUBA_SRC_DIR}/common/Data.c
    ${CUBA_SRC_DIR}/common/Global.c
)
target_include_directories(cuba PUBLIC "${CUBA_ROOT}")
target_compile_definitions(cuba PRIVATE NOUNDERSCORE=1)
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
target_link_libraries(subnucleondiffraction PRIVATE cuba)
if(OpenMP_CXX_FOUND)
	target_link_libraries(subnucleondiffraction PRIVATE OpenMP::OpenMP_CXX)
endif()
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
target_link_libraries(test PRIVATE cuba)
if(OpenMP_CXX_FOUND)
	target_link_libraries(test PRIVATE OpenMP::OpenMP_CXX)
endif()



