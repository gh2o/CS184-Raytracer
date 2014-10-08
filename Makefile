CXX := clang++
LIBS := build/libpng/.libs/libpng16.a
CXXFLAGS := -std=c++11 -O2 -Wall -Ibuild/libpng -Ieigen -Ilibpng
LDFLAGS := -lz

as2: as2.cpp $(LIBS)
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

build/libpng/.libs/libpng16.a:
	rm -rf build/libpng
	mkdir -p build/libpng
	cd build/libpng && ../../libpng/configure --enable-static --disable-shared
	make -C build/libpng -j 4 libpng16.la
