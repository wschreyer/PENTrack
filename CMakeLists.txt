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
find_package(Boost 1.53.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
else()
	message(SEND_ERROR "Boost 1.53.0 or newer not found!")
endif()

include_directories("exprtk")
include_directories("include")
add_executable(PENTrack src/main.cpp src/globals.cpp src/trianglemesh.cpp src/geometry.cpp src/mc.cpp src/edmfields.cpp 
                        src/field_2d.cpp src/field_3d.cpp src/fields.cpp src/harmonicfields.cpp src/conductor.cpp src/particle.cpp src/neutron.cpp src/microroughness.cpp
                        src/electron.cpp src/proton.cpp src/mercury.cpp src/xenon.cpp src/source.cpp src/config.cpp src/analyticFields.cpp
                        $<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:libtricubic>)

if (CMAKE_COMPILER_IS_GNUCXX)
	target_compile_options(PENTrack PUBLIC -Wall)
endif()
	
				
target_link_libraries (PENTrack ${Boost_LIBRARIES} ${CGAL_LIBRARIES})
