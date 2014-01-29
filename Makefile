SRC = main.c kdtree.c
OBJ=$(SRC:.c=.o)

CGAL_INCLUDE= #-I$(HOME)/CGAL-4.2/include # point gcc's -I option to CGAL include directory if you have compiled CGAL manually without installing it
CGAL_LIB= #-L$(HOME)/CGAL-4.2/lib # point gcc's -L option to CGAL lib directory if you have compiled CGAL manually without installing it
CGAL_SHAREDLIB= #-Wl,-rpath=$(HOME)/CGAL-4.2/lib # point gcc's -Wl,-rpath= option to CGAL shared library if you have compiled CGAL manually without installing it

MUPARSER_INCLUDE= #-I$(HOME)/muparser_v2_2_3/include/ # point gcc's -I option to muparser include directory if you have compiled muparser manually without installing it
MUPARSER_LIB= #-L$(HOME)/muparser_v2_2_3/lib/ # point gcc's -L option to muparser lib directory if you have compiled muparser manually without installing it
MUPARSER_SHAREDLIB= #-Wl,-rpath=$(HOME)/muparser_v2_2_3/lib/ # point gcc's -Wl,-rpath= option to muparser shared library if you have compiled muparser manually without installing it

CC=g++
CFLAGS=-O2 -frounding-math -Wall -Wno-reorder -Wno-parentheses -Wno-strict-aliasing $(CGAL_INCLUDE) $(CGAL_SHAREDLIB) $(MUPARSER_INCLUDE) $(MUPARSER_SHAREDLIB) #-O2: optimize, -Wno-*: suppress warnings from external libraries
LDFLAGS=-lrt -lboost_system $(CGAL_LIB) -lCGAL $(MUPARSER_LIB) -lmuparser
RM=rm
EXE=PENTrack

ADDSRC = libtricubic/libtricubic.cpp libtricubic/tricubic_utils.cpp 
ADDOBJ = $(ADDSRC:.cpp=.o)

.PHONY: all
all: $(SRC) $(ADDOBJ)
	$(CC) -o $(EXE) $(SRC) $(ADDOBJ) $(CFLAGS) $(LDFLAGS)
	
%.o:%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	$(RM) $(EXE) $(ADDOBJ)
