SRC = main.c kdtree.c libtricubic/*.cpp
OBJ=$(SRC:.c=.o)

CGAL_INCLUDE= #$(HOME)/CGAL-4.2/include	# change path to CGAL include directory if you have compiled CGAL manually without installing it
CGAL_LIB= #$(HOME)/CGAL-4.2/lib			# change path to CGAL libraries if you have compiled CGAL manually without installing it

CC=g++
CFLAGS=-O2 -frounding-math -I$(CGAL_INCLUDE) -Wl,-rpath=$(CGAL_LIB) #-Wall #-pedantic #-g # -O3: optimize -g: debug switch
LDFLAGS=-lrt -lboost_system -L$(CGAL_LIB) -lCGAL
RM=rm
EXE=PENTrack

.PHONY: all
all:
	$(CC) $(SRC) -o $(EXE) $(CFLAGS) $(LDFLAGS)
	
.PHONY: clean
clean:
	$(RM) $(EXE)
