SRC = main.c kdtree.c libtricubic/*.cpp
OBJ=$(SRC:.c=.o)

CGAL_INCLUDE= #-I$(HOME)/CGAL-4.2/include	# point gcc's -I option to CGAL include directory if you have compiled CGAL manually without installing it
CGAL_LIB= #-L$(HOME)/CGAL-4.2/lib			# point gcc's -L option to CGAL lib directory if you have compiled CGAL manually without installing it
CGAL_SHAREDLIB= #-Wl,-rpath=$(HOME)/CGAL-4.2/lib			# point gcc's -Wl,-rpath= option to CGAL shared library if you have compiled CGAL manually without installing it

CC=g++
CFLAGS=-O2 -frounding-math $(CGAL_INCLUDE) $(CGAL_SHAREDLIB) #-Wall #-pedantic #-g # -O3: optimize -g: debug switch
LDFLAGS=-lrt -lboost_system $(CGAL_LIB) -lCGAL
RM=rm
EXE=PENTrack

.PHONY: all
all:
	$(CC) $(SRC) -o $(EXE) $(CFLAGS) $(LDFLAGS)
	
.PHONY: clean
clean:
	$(RM) $(EXE)
