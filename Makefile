CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -Wno-conversion-null -Wno-deprecated-declarations -I. -I/home/miguelalcaniz/Escritorio/PACS/pacs-examples/Examples/include -I/home/miguelalcaniz/Escritorio/PACS/pacs-examples/Examples/include/muparserx -I/home/miguelalcaniz/Escritorio/PACS/JSON/json/single_include/nlohmann -I/home/miguelalcaniz/Escritorio/PACS/pacs-examples/Examples/src/muParserInterface
EXEC      = main

LDFLAGS ?= -L/home/miguelalcaniz/Escritorio/PACS/pacs-examples/Examples/lib
LIBS  ?= -lmuparser -lmuparserx

all: $(EXEC)

%.o: %.cpp NewtonSolver.hpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $<

$(EXEC): %: %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< $(LIBS) -o $@

clean:
	$(RM) *.o

distclean: clean
	$(RM) *~
