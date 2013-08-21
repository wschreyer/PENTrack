SRC = main.c 2dinterpolfeld.c racetrack.c globals.c bruteforce.c adiabacity.c ndist.c mersenne/mt.c nr/stepperdopr853.c kdtree.c libtricubic/*.cpp #BForbes.c
OBJ=$(SRC:.c=.o)


CC=g++
CFLAGS=-O3 -Wall #`root-config --cflags` # -O3: optimize -g: debug switch -std=c99 -pedantic
LDFLAGS=#-lc `root-config --libs` # -lm #/lib/mingw/libmingwex.a /lib/mingw/libmsvcrt.a
RM=rm
EXE=Track2013

#%.o: %.c   		 # combined w/ next line will compile recently changed .c files
#	$(CC) $(CFLAGS) -o $@ -c $<

#.PHONY : all     # .PHONY ignores files named all

#all: $(EXE)      # all is dependent on $(EXE) to be complete

#$(EXE): $(OBJ)   # $(EXE) is dependent on all of the files in $(OBJ) to exist
#	$(CC) $(OBJ) $(LDFLAGS) -o $@

#.PHONY : clean   # .PHONY ignores files named clean
#clean:
#	-$(RM) $(OBJ)    # '-' causes errors not to exit the process

.PHONY: all
all:
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC) -o $(EXE)
	
.PHONY: clean
clean:
	$(RM) $(EXE)
