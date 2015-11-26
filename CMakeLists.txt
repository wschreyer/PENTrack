cmake_minimum_required (VERSION 2.6)
project (PENTrack)
set(CMAKE_BUILD_TYPE Release)

include_directories("alglib-3.9.0/cpp/src")
add_library(alglib OBJECT alglib-3.9.0/cpp/src/alglibinternal.cpp 
						alglib-3.9.0/cpp/src/dataanalysis.cpp alglib-3.9.0/cpp/src/integration.cpp 
						alglib-3.9.0/cpp/src/optimization.cpp alglib-3.9.0/cpp/src/statistics.cpp 
						alglib-3.9.0/cpp/src/alglibmisc.cpp alglib-3.9.0/cpp/src/diffequations.cpp 
						alglib-3.9.0/cpp/src/interpolation.cpp alglib-3.9.0/cpp/src/solvers.cpp 
						alglib-3.9.0/cpp/src/ap.cpp alglib-3.9.0/cpp/src/fasttransforms.cpp 
						alglib-3.9.0/cpp/src/linalg.cpp alglib-3.9.0/cpp/src/specialfunctions.cpp)						
						
include_directories("muparser_v2_2_4/include")
add_library(muparser OBJECT muparser_v2_2_4/src/muParser.cpp muparser_v2_2_4/src/muParserBase.cpp 
						muparser_v2_2_4/src/muParserBytecode.cpp muparser_v2_2_4/src/muParserCallback.cpp 
						muparser_v2_2_4/src/muParserError.cpp muparser_v2_2_4/src/muParserTokenReader.cpp)						
						
include_directories("libtricubic")
add_library(libtricubic OBJECT libtricubic/libtricubic.cpp libtricubic/tricubic_utils.cpp)

find_package(Boost 1.53.0)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
else()
	message(SEND_ERROR "Boost 1.53.0 or newer not found!")
endif()

find_package(CGAL)
if (CGAL_CONFIG_LOADED)
	if (CGAL_VERSION VERSION_LESS 4.1)
		message(SEND_ERROR "CGAL version 4.1 or newer not found!")
	endif()
	include_directories(${CGAL_INCLUDE_DIRS})
	link_directories(${CGAL_LIBRARIES_DIR})
	message(STATUS "CGAL version: " ${CGAL_VERSION})
else()
	message(SEND_ERROR "CGAL library not found!")
endif()
						


add_executable(PENTrack main.cpp globals.cpp trianglemesh.cpp geometry.cpp mc.cpp bruteforce.cpp 
						field_2d.cpp field_3d.cpp fields.cpp conductor.cpp particle.cpp neutron.cpp 
						electron.cpp proton.cpp ndist.cpp source.cpp
						$<TARGET_OBJECTS:alglib> $<TARGET_OBJECTS:muparser> $<TARGET_OBJECTS:libtricubic>)		

if (CMAKE_COMPILER_IS_GNUCXX)
	set_target_properties(PENTrack PROPERTIES COMPILE_FLAGS "-frounding-math -Wall -Wno-unused-local-typedefs")
	if (Boost_VERSION GREATER 105500)
		set_target_properties(PENTrack PROPERTIES COMPILE_FLAGS "-std=c++11 -frounding-math -Wall -Wno-unused-local-typedefs -Wno-deprecated-declarations")
		message(STATUS "Boost > 1.55.0 requires C++11. I hope it is supported by your compiler.")
	endif()
endif()
	
				
target_link_libraries (PENTrack rt boost_system CGAL)
