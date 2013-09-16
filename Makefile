SRC = main.c 2dinterpolfeld.c racetrack.c globals.c bruteforce.c adiabacity.c ndist.c mersenne/mt.c nr/stepperdopr853.c kdtree.c libtricubic/*.cpp #BForbes.c
OBJ=$(SRC:.c=.o)


CC=g++
CFLAGS=-O3 -Wall #`root-config --cflags` # -O3: optimize -g: debug switch -std=c99 -pedantic
LDFLAGS=#-lc `root-config --libs` # -lm #/lib/mingw/libmingwex.a /lib/mingw/libmsvcrt.a
RM=rm
EXE=PENTrack

.PHONY: all
all:
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC) -o $(EXE)
	
.PHONY: clean
clean:
	$(RM) $(EXE)
