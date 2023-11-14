CXX := g++
CXXFLAGS := -I./common -O3 -std=c++17 -Wall -lpthread

ifeq (/usr/bin/g++-10,$(wildcard /usr/bin/g++-10*))
    CXX=g++-10
endif

all: pi

pi: pi.cpp
	$(CXX) pi.cpp $(CXXFLAGS) -o ./pi.out


clean:
	rm -f *.out 