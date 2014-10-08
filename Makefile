CXX := clang++
LIBS := build/libpng/.libs/libpng16.a
CFLAGS := -std=c++11 -O2 -Wall -Ieigen -Ilibpng
LDFLAGS :=

as2: as2.cpp $(LIBS)
	$(CXX) $^ -o $@ $(CFLAGS) $(LDFLAGS) $(LIBS)

build/libpng/.libs/libpng16.a:
	rm -rf build/libpng
	mkdir -p build/libpng
	cd build/libpng && ../../libpng/configure --enable-static --disable-shared
	make -C build/libpng libpng16.la
