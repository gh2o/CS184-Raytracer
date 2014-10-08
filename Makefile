CFLAGS := -std=c++11 -O2 -Wall -Ieigen -Wno-unused-local-typedefs
LDFLAGS := 

as2: as2.cpp
	g++ $^ -o $@ $(CFLAGS) $(LDFLAGS)
