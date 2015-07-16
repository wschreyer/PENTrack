SRC = main.cpp globals.cpp trianglemesh.cpp geometry.cpp mc.cpp bruteforce.cpp \
		field_2d.cpp field_3d.cpp fields.cpp conductor.cpp particle.cpp neutron.cpp electron.cpp proton.cpp ndist.cpp source.cpp
OBJ = $(SRC:.cpp=.o)

TRICUBICSRC = libtricubic/libtricubic.cpp libtricubic/tricubic_utils.cpp
TRICUBICOBJ = $(TRICUBICSRC:.cpp=.o)
ALGLIBSRC = alglib-3.9.0/cpp/src/alglibinternal.cpp alglib-3.9.0/cpp/src/dataanalysis.cpp alglib-3.9.0/cpp/src/integration.cpp alglib-3.9.0/cpp/src/optimization.cpp alglib-3.9.0/cpp/src/statistics.cpp alglib-3.9.0/cpp/src/alglibmisc.cpp alglib-3.9.0/cpp/src/diffequations.cpp alglib-3.9.0/cpp/src/interpolation.cpp alglib-3.9.0/cpp/src/solvers.cpp alglib-3.9.0/cpp/src/ap.cpp alglib-3.9.0/cpp/src/fasttransforms.cpp alglib-3.9.0/cpp/src/linalg.cpp alglib-3.9.0/cpp/src/specialfunctions.cpp
ALGLIBOBJ = $(ALGLIBSRC:.cpp=.o)
MUPARSERSRC = muparser_v2_2_4/src/muParser.cpp muparser_v2_2_4/src/muParserBase.cpp muparser_v2_2_4/src/muParserBytecode.cpp muparser_v2_2_4/src/muParserCallback.cpp muparser_v2_2_4/src/muParserError.cpp muparser_v2_2_4/src/muParserTokenReader.cpp
MUPARSEROBJ = $(MUPARSERSRC:.cpp=.o)

BOOST_INCLUDE = #-I$(HOME)/boost_1_55_0/include # point gcc's -I option to Boost include directory if you have compiled Boost manually without installing it
BOOST_LIB = #-L$(HOME)/boost_1_55_0/stage/lib # point gcc's -L option to Boost lib directory if you have compiled Boost manually without installing it
BOOST_SHAREDLIB = #-Wl,-rpath=$(HOME)/boost_1_55_0/stage/lib # point gcc's -Wl,-rpath= option to Boost shared library if you have compiled Boost manually without installing it

CGAL_INCLUDE = #-I$(HOME)/CGAL-4.6/include # point gcc's -I option to CGAL include directory if you have compiled CGAL manually without installing it
CGAL_LIB = #-L$(HOME)/CGAL-4.6/lib # point gcc's -L option to CGAL lib directory if you have compiled CGAL manually without installing it
CGAL_SHAREDLIB = #-Wl,-rpath=$(HOME)/CGAL-4.6/lib # point gcc's -Wl,-rpath= option to CGAL shared library if you have compiled CGAL manually without installing it

CC=g++
LDFLAGS=-lrt -lboost_system $(BOOST_LIB) $(CGAL_LIB) -lCGAL
RM=rm
EXE=PENTrack

.PHONY: all
all: $(OBJ) $(TRICUBICOBJ) $(ALGLIBOBJ) $(MUPARSEROBJ)
	$(CC) -o $(EXE) $(OBJ) $(TRICUBICOBJ) $(ALGLIBOBJ) $(MUPARSEROBJ) $(CFLAGS) $(LDFLAGS)
	
$(OBJ): CFLAGS = -O3 -frounding-math -Wall -Ilibtricubic -Ialglib-3.9.0/cpp/src -Imuparser_v2_2_4/include $(BOOST_INCLUDE) $(BOOST_SHAREDLIB) $(CGAL_INCLUDE) $(CGAL_SHAREDLIB) #-O2: optimize, -Wno-*: suppress warnings from external libraries

$(TRICUBICOBJ): CFLAGS = -O3 -Wall -Ilibtricubic

$(ALGLIBOBJ): CFLAGS = -O3 -Ialglib-3.9.0/cpp/src
	
$(MUPARSEROBJ): CFLAGS = -O3 -Wall -Wno-switch -Imuparser_v2_2_4/include
	
%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	$(RM) $(EXE) $(OBJ) $(ALGLIBOBJ) $(MUPARSEROBJ)
