CC      := g++
CFLAGS  := -std=c++14 -g
LD      := ld
LDFLAGS := 
SOURCE  := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp, %.o, $(SOURCE))

all:$(OBJS)

%.o:%.cpp
	@echo Compiling $< ...
	$(CC) -c $(CFLAGS) $< -o $*.o
 
clean:
	@rm -f *.o 

