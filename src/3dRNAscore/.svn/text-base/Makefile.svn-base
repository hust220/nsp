CC      := g++
CFLAGS  := 
LD      := ld
LDFLAGS := 
SOURCE  := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp, %.o, $(SOURCE))

all:3dRNAscore
	cp 3dRNAscore ${H3DRNA}/bin

3dRNAscore:$(OBJS) ../lib3dRNA.a
	g++ $(OBJS) ../lib3dRNA.a -o 3dRNAscore

%.o:%.cpp
	@echo Compiling $< ...
	$(CC) -c $(CFLAGS) $< -o $*.o
 
clean:
	@rm -f *.o 3dRNAscore

