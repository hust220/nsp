CXX      := g++
CXXFLAGS  := -std=c++14 -g
LD      := ld
LDFLAGS := 
SOURCE  := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp, %.o, $(SOURCE))

all: 3dRNA
	cp 3dRNA ${H3DRNA}/bin/3dRNA

3dRNA:$(OBJS) ../lib3dRNA.a
	@echo 3dRNA ...
	$(CXX) $(CXXFLAGS) $(OBJS) ../lib3dRNA.a -o 3dRNA

%.o:%.cpp %.h
	@echo Compiling $< ...
	$(CXX) -c $(CXXFLAGS) $< -o $*.o
 
clean:
	@rm -f *.o 3dRNA

