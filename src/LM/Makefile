CXX      := g++
CXXFLAGS  := -std=c++14 -g
LD      := ld
LDFLAGS := 

all:lm5.o
	@echo link ...
	$(CXX) $(CXXFLAGS) lm5.o ../lib3dRNA.a -o lm5
	cp lm5 ${H3DRNA}/bin
 
lm5.o:lm5.cpp
	$(CXX) $(CXXFLAGS) -c lm5.cpp

clean:
	@rm -f lm5.o lm5

