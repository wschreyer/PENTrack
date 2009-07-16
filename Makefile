SRC = bahnen2005.c bruteforce.c cel.c 2dinterpolfeld.c nrutil.c adiabacity.c integrator.c main.c mc.c neutron.c mersenne/mt.c maplefield.c BForbes.c geometry.c kdtree.c
OBJ=$(SRC:.c=.o)


CC=g++
CFLAGS=-Wall  `root-config --cflags`  -o3   -g  -Wno-deprecated    # -o3: optimize -g: debug switch -std=c99 -pedantic
LDFLAGS=-lm -lc `root-config --libs` #/lib/mingw/libmingwex.a /lib/mingw/libmsvcrt.a
RM=rm
EXE=Track2009

%.o: %.c   		 # combined w/ next line will compile recently changed .c files
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : all     # .PHONY ignores files named all

all: $(EXE)      # all is dependent on $(EXE) to be complete

$(EXE): $(OBJ)   # $(EXE) is dependent on all of the files in $(OBJ) to exist
	$(CC) $(OBJ) $(LDFLAGS) -o $@

.PHONY : clean   # .PHONY ignores files named clean
clean:
	-$(RM) $(OBJ) core    # '-' causes errors not to exit the process


