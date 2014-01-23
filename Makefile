SRC = main.c kdtree.c
OBJ=$(SRC:.c=.o)

MUPARSER_INCLUDE= #-I$(HOME)/muparser_v2_2_3/include/ # point gcc's -I option to muparser include directory if you have compiled muparser manually without installing it
MUPARSER_LIB= #-L$(HOME)/muparser_v2_2_3/lib/ # point gcc's -L option to muparser lib directory if you have compiled muparser manually without installing it
MUPARSER_SHAREDLIB= #-Wl,-rpath=$(HOME)/muparser_v2_2_3/lib/ # point gcc's -Wl,-rpath= option to muparser shared library if you have compiled muparser manually without installing it

CC=g++
CFLAGS=-O2 $(MUPARSER_INCLUDE) $(MUPARSER_SHAREDLIB) #-Wall #-pedantic #-g # -O3: optimize -g: debug switch
LDFLAGS=-lrt $(MUPARSER_LIB) -lmuparser
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
