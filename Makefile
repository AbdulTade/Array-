CXX=g++.exe
CXXFLAGS=-O3 -fpermissive -std=c++17
DEBUGGER=gdb.exe

array.exe: array.cpp
	$(CXX) $(CXXFLAGS) array.cpp -o ./build/array.exe

clean:
	rm ./build/array.exe

debug:
	$(DEBUGGER) ./build/array.exe