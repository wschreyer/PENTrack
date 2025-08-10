cmake_minimum_required(VERSION 3.15)
project(PENTrack)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

add_library(alglib SHARED alglib-3.15.0/cpp/src/alglibinternal.cpp
                          alglib-3.15.0/cpp/src/integration.cpp
                          alglib-3.15.0/cpp/src/optimization.cpp
                          alglib-3.15.0/cpp/src/alglibmisc.cpp
                          alglib-3.15.0/cpp/src/interpolation.cpp
                          alglib-3.15.0/cpp/src/solvers.cpp
                          alglib-3.15.0/cpp/src/ap.cpp
                          alglib-3.15.0/cpp/src/linalg.cpp
                          alglib-3.15.0/cpp/src/specialfunctions.cpp)
target_include_directories(alglib PUBLIC alglib-3.15.0/cpp/src)
						
add_library(tricubic SHARED tricubic-1.0/src/libtricubic/libtricubic.cpp tricubic-1.0/src/libtricubic/tricubic_utils.cpp)
target_include_directories(tricubic PUBLIC tricubic-1.0/src/libtricubic)

find_package(CGAL 5.3 REQUIRED HINTS cgal COMPONENTS Core)
find_package(Boost 1.73.0 CONFIG REQUIRED system filesystem iostreams OPTIONAL_COMPONENTS unit_test_framework)

add_library(PENTrack_src OBJECT src/globals.cpp src/geometry.cpp src/mc.cpp src/field.cpp src/edmfields.cpp src/tracking.cpp src/logger.cpp
                                src/field_2d.cpp src/field_3d.cpp src/fields.cpp src/harmonicfields.cpp src/conductor.cpp src/particle.cpp src/neutron.cpp src/microroughness.cpp
                                src/electron.cpp src/proton.cpp src/mercury.cpp src/xenon.cpp src/source.cpp src/config.cpp src/analyticFields.cpp src/mesh.cpp)
target_include_directories(PENTrack_src PUBLIC exprtk include)
target_link_libraries(PENTrack_src PUBLIC alglib tricubic Boost::boost Boost::system Boost::filesystem Boost::iostreams CGAL::CGAL CGAL::CGAL_Core)


find_package(ROOT COMPONENTS Core Tree)
if (ROOT_FOUND)
	target_compile_definitions(PENTrack_src PUBLIC USEROOT=1)
	target_link_libraries(PENTrack_src PUBLIC ROOT::Core ROOT::Tree)
else()
	message(WARNING "Could not find ROOT, you won't be able to log to ROOT files")
endif()

find_package(HDF5 COMPONENTS C HL)
if (HDF5_FOUND)
	target_compile_definitions(PENTrack_src PUBLIC USEHDF5=1)
	target_link_libraries(PENTrack_src PUBLIC HDF5::HDF5 hdf5::hdf5_hl)
else()
	message(WARNING "Could not find HDF5 library, you won't be able to log to HDF5 files")
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
	target_compile_options(PENTrack_src PRIVATE -Wall)
endif()


add_executable(PENTrack src/main.cpp)
target_link_libraries(PENTrack PENTrack_src)


if (TARGET Boost::unit_test_framework)
	message(STATUS "Found Boost::unit_test_framwork. Tests will be built")
	enable_testing()
	add_executable(runTests test/test.cpp test/fieldTests.cpp test/scatteringTest.cpp test/geometryTest.cpp)
	target_link_libraries(runTests PUBLIC PENTrack_src Boost::unit_test_framework)
	target_compile_definitions(runTests PRIVATE "BOOST_TEST_DYN_LINK=1")
	add_test(NAME Boost_unit_tests COMMAND runTests)
endif()
