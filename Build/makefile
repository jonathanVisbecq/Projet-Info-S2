OBJS = main.o p_adic.o scripts.o
INCLUDES = ../Include

CXX = g++
DEBUG = -g
LFLAGS = -Wall $(DEBUG)
CXXFLAGS = -Wall -c $(DEBUG) -std=c++0x -I $(INCLUDES)

vpath %.hpp ../Include
vpath %.cpp ../Source

all : projetPages.exe

projetPages.exe: $(OBJS)
	$(CXX) $(LFLAGS) $(OBJS) -o projetPages.exe


main.o:
p_adic.o: p_adic.hpp
scripts.o: scripts.hpp


clean:
	\rm *.o *~ projetPages.exe
