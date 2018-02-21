.PHONY: all clean test

CXX = g++
CXXFLAGS = -Wall -Wextra -fsanitize=undefined -std=c++14 -O3 -Iinclude
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
