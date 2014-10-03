CXX=g++
CXXARGS=-Wall -O3 -g
all: alignAcc

alignAcc: alignAcc.o StringRef.o miscFunctions.o
	$(CXX) $(CXXARGS) -o alignAcc alignAcc.o StringRef.o miscFunctions.o

alignAcc.o: alignAcc.cpp
	$(CXX) $(CXXARGS) -c -o alignAcc.o alignAcc.cpp

StringRef.o: StringRef.cpp StringRef.h
	$(CXX) $(CXXARGS) -c -o StringRef.o StringRef.cpp

miscFunctions.o: miscFunctions.cpp miscFunctions.h
	$(CXX) $(CXXARGS) -c -o miscFunctions.o miscFunctions.cpp

clean:
	rm alignAcc