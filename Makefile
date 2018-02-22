.PHONY: all clean test

CXX = g++
CXXFLAGS = -Wall -Wextra -fsanitize=undefined -std=c++1z -O3 -Iinclude -march=native
LDFLAGS =

OBJS =

TESTSRC = $(wildcard test/*.cpp)
TESTOUT = $(TESTSRC:.cpp=.out)

all: test

clean:
	rm src/*.o test/*.out ; true

test: $(TESTOUT)

%.out: %.cpp
	$(CXX) $(CXXFLAGS) -Os -o $(OBJS) $@ $<
	echo "Running $@"
	$@
