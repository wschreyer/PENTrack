cmake_minimum_required (VERSION 2.8.11)
project (PENTrack CXX)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

include_directories("alglib-3.10.0/cpp/src")
add_library(alglib OBJECT alglib-3.10.0/cpp/src/alglibinternal.cpp 
						alglib-3.10.0/cpp/src/dataanalysis.cpp alglib-3.10.0/cpp/src/integration.cpp 
						alglib-3.10.0/cpp/src/optimization.cpp alglib-3.10.0/cpp/src/statistics.cpp 
						alglib-3.10.0/cpp/src/alglibmisc.cpp alglib-3.10.0/cpp/src/diffequations.cpp 
						alglib-3.10.0/cpp/src/interpolation.cpp alglib-3.10.0/cpp/src/solvers.cpp 
						alglib-3.10.0/cpp/src/ap.cpp alglib-3.10.0/cpp/src/fasttransforms.cpp 
						alglib-3.10.0/cpp/src/linalg.cpp alglib-3.10.0/cpp/src/specialfunctions.cpp)						
						
include_directories("muparser-2.2.5/include")
add_library(muparser OBJECT muparser-2.2.5/src/muParser.cpp muparser-2.2.5/src/muParserBase.cpp 
						muparser-2.2.5/src/muParserBytecode.cpp muparser-2.2.5/src/muParserCallback.cpp 
						muparser-2.2.5/src/muParserError.cpp muparser-2.2.5/src/muParserTokenReader.cpp)						
						
include_directories("tricubic-1.0/src/libtricubic")
add_library(libtricubic OBJECT tricubic-1.0/src/libtricubic/libtricubic.cpp tricubic-1.0/src/libtricubic/tricubic_utils.cpp)

set(Boost_REALPATH ON)
find_package(Boost 1.53.0 COMPONENTS system)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else()
	message(SEND_ERROR "Boost 1.53.0 or newer not found!")
endif()

find_package(CGAL PATHS "cgal/build")
if (CGAL_CONFIG_LOADED)
	message(STATUS "CGAL version: " ${CGAL_VERSION})
	if (CGAL_VERSION VERSION_LESS 4.1)
		message(SEND_ERROR "CGAL versions older than 4.1 are not supported!")
	endif()
	include_directories(${CGAL_INCLUDE_DIRS})
	include(${CGAL_USE_FILE})
else()
	message(SEND_ERROR "CGAL library not found! For installation instructions see README.md.")
endif()

find_library(RT_LIB rt)

include_directories("include")
add_executable(PENTrack src/main.cpp src/globals.cpp src/trianglemesh.cpp src/geometry.cpp src/mc.cpp src/edmfields.cpp 
						src/field_2d.cpp src/field_3d.cpp src/fields.cpp src/conductor.cpp src/particle.cpp src/neutron.cpp src/microroughness.cpp
						src/electron.cpp src/proton.cpp src/mercury.cpp src/xenon.cpp src/source.cpp 
						$<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:muparser> $<TARGET_OBJECTS:libtricubic>)		

if (CMAKE_COMPILER_IS_GNUCXX)
	target_compile_options(PENTrack PUBLIC -Wall PUBLIC -Wno-unused-local-typedefs)
	if (Boost_VERSION GREATER 105500)
		set_property(TARGET PENTrack PROPERTY CXX_STANDARD_REQUIRED ON)
		set_property(TARGET PENTrack PROPERTY CXX_STANDARD 11)
		target_compile_options(PENTrack PUBLIC -Wno-deprecated-declarations)
		message(STATUS "Boost 1.56.0+ requires C++11. I hope it is supported by your compiler.")
	endif()
endif()
	
				
target_link_libraries (PENTrack ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${RT_LIB})
