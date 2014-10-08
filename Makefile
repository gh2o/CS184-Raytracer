CFLAGS := -std=c++11 -Os -Wall
LDFLAGS := 

as2: as2.cpp
	g++ $^ -o $@ $(CFLAGS) $(LDFLAGS)
