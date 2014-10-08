CXX := clang++
CFLAGS := -std=c++11 -O2 -Wall -Ieigen
LDFLAGS := 

as2: as2.cpp
	$(CXX) $^ -o $@ $(CFLAGS) $(LDFLAGS)
