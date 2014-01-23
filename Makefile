SRC = main.c kdtree.c
OBJ=$(SRC:.c=.o)

MUPARSER_INCLUDE=$(HOME)/muparser_v2_2_3/include/ # uncomment and adapt if you have built muparser manually without installing it
MUPARSER_LIB=$(HOME)/muparser_v2_2_3/lib/

CC=g++
CFLAGS=-O2 -I$(MUPARSER_INCLUDE) -Wl,-rpath=$(MUPARSER_LIB) #-Wall #-pedantic #-g # -O3: optimize -g: debug switch
LDFLAGS=-lrt -L$(MUPARSER_LIB) -lmuparser
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
