cmake_minimum_required (VERSION 3.1)
project (PENTrack CXX)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("alglib-3.15.0/cpp/src")
add_library(alglib OBJECT alglib-3.15.0/cpp/src/alglibinternal.cpp
                          alglib-3.15.0/cpp/src/integration.cpp
                          alglib-3.15.0/cpp/src/optimization.cpp
                          alglib-3.15.0/cpp/src/alglibmisc.cpp
                          alglib-3.15.0/cpp/src/interpolation.cpp
                          alglib-3.15.0/cpp/src/solvers.cpp
                          alglib-3.15.0/cpp/src/ap.cpp
                          alglib-3.15.0/cpp/src/linalg.cpp
                          alglib-3.15.0/cpp/src/specialfunctions.cpp)
						
include_directories("tricubic-1.0/src/libtricubic")
add_library(libtricubic OBJECT tricubic-1.0/src/libtricubic/libtricubic.cpp tricubic-1.0/src/libtricubic/tricubic_utils.cpp)

include_directories("exprtk")
include_directories("include")

set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Tell CGAL not to override compiler flags")
find_package(CGAL HINTS "cgal")
if (CGAL_FOUND)
	message(STATUS "CGAL version: " ${CGAL_VERSION})
	if (CGAL_VERSION VERSION_LESS 4.14)
		message(SEND_ERROR "CGAL versions older than 4.14 are not supported!")
	endif()
	include_directories(${CGAL_INCLUDE_DIRS})
	include(${CGAL_USE_FILE})
else()
	message(SEND_ERROR "CGAL library not found! For installation instructions see README.md.")
endif()

set(Boost_REALPATH ON)
if (BUILD_TESTS)
	message(STATUS "Tests will be built, requires Boost 1.59.0 or newer with unit test framework")
	find_package(Boost 1.59.0 REQUIRED system filesystem unit_test_framework)
else()
	find_package(Boost 1.53.0 REQUIRED system filesystem)
endif()
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
endif()

find_package(ROOT)
if (ROOT_FOUND)
	message(STATUS "Found ROOT, you can use the ROOTlog option")
    include(${ROOT_USE_FILE})
else()
	message(STATUS "Cound not find ROOT, you won't be able to use the ROOTlog option")
endif()

				
add_library(PENTrack_src OBJECT src/globals.cpp src/trianglemesh.cpp src/geometry.cpp src/mc.cpp src/field.cpp src/edmfields.cpp src/tracking.cpp src/logger.cpp src/stepper.cpp
                        		src/field_2d.cpp src/field_3d.cpp src/fields.cpp src/harmonicfields.cpp src/conductor.cpp src/particle.cpp src/neutron.cpp src/microroughness.cpp
                        		src/electron.cpp src/proton.cpp src/mercury.cpp src/xenon.cpp src/source.cpp src/config.cpp src/analyticFields.cpp)

if (ROOT_FOUND)
	target_compile_definitions(PENTrack_src PUBLIC USEROOT=1)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
	target_compile_options(PENTrack_src PUBLIC -Wall)
endif()


add_executable(PENTrack src/main.cpp $<TARGET_OBJECTS:PENTrack_src> $<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:libtricubic>)
target_link_libraries (PENTrack ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${ROOT_LIBRARIES})


if (BUILD_TESTS)
	enable_testing()
	add_executable(runTests test/test.cpp test/fieldTests.cpp $<TARGET_OBJECTS:PENTrack_src> $<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:libtricubic>)
	target_link_libraries(runTests ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${ROOT_LIBRARIES})
	target_compile_definitions(runTests PRIVATE "BOOST_TEST_DYN_LINK=1")
	add_test(COMMAND runTests)
endif()
