SRC = main.c kdtree.c
OBJ=$(SRC:.c=.o)

CC=g++
CFLAGS=-O2 -Imuparser_v2_2_3/include #-Wall #-pedantic #-g # -O3: optimize -g: debug switch
LDFLAGS=-lrt
RM=rm
EXE=PENTrack

ADDSRC = libtricubic/libtricubic.cpp libtricubic/tricubic_utils.cpp muparser_v2_2_3/src/muParser.cpp muparser_v2_2_3/src/muParserBase.cpp muparser_v2_2_3/src/muParserBytecode.cpp muparser_v2_2_3/src/muParserCallback.cpp muparser_v2_2_3/src/muParserError.cpp muparser_v2_2_3/src/muParserTokenReader.cpp 
ADDOBJ = $(ADDSRC:.cpp=.o)

.PHONY: all
all: $(SRC) $(ADDOBJ)
	$(CC) -o $(EXE) $(SRC) $(ADDOBJ) $(CFLAGS) $(LDFLAGS)
	
%.o:%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	$(RM) $(EXE) $(ADDOBJ)
