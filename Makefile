SRC = main.c kdtree.c libtricubic/*.cpp
OBJ=$(SRC:.c=.o)


CC=g++
CFLAGS=-O2 -frounding-math #-Wall #-pedantic #-g # -O3: optimize -g: debug switch
LDFLAGS=-lrt -lboost_system -lCGAL
RM=rm
EXE=PENTrack

.PHONY: all
all:
	$(CC) $(SRC) -o $(EXE) $(CFLAGS) $(LDFLAGS)
	
.PHONY: clean
clean:
	$(RM) $(EXE)
