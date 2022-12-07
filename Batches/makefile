OBJS = Particle.o main.o
EXE = CoaleBatches

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
# LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++ -std=c++17
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)
%.o: %.cxx
	$(COMPILE)  $< 

.PHONY : clean
clean :
	-rm -f $(OBJS) $(EXE) *.gch
