# Makefile for StoutSmearing project

CXX = g++
CXXFLAGS = -std=c++20 -O2 -Wall
LDFLAGS = -llapacke -llapack -lblas -lboost_program_options

SRCS = main.cpp gauge_field.cpp linalg.cpp stout_smearing.cpp su3maximization.cpp
OBJS = $(SRCS:.cpp=.o)

EXEC = main

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJS) $(EXEC)
