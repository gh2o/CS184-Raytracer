CXX := clang++
LIBS := build/libpng/.libs/libpng16.a
CXXFLAGS := -std=c++11 -O2 -Wall -Ibuild/libpng -Ieigen -Ilibpng
LDFLAGS := -lz

SRCS := $(wildcard src/*.cpp)
OBJS := $(patsubst src/%.cpp,build/%.o,$(SRCS))
DEPS := $(patsubst src/%.cpp,build/%.d,$(SRCS))

as2: $(LIBS) $(OBJS)
	$(CXX) $(OBJS) $(LIBS) $(CXXFLAGS) $(LDFLAGS) -o $@

build/%.o: src/%.cpp build/%.d
	$(CXX) $< $(CXXFLAGS) -c -o $@

build/%.d: src/%.cpp Makefile
	@mkdir -p build
	@$(CXX) $< $(CXXFLAGS) -MM -MP -MT $@ -MF $@

build/libpng/.libs/libpng16.a:
	rm -rf build/libpng
	mkdir -p build/libpng
	cd build/libpng && ../../libpng/configure --enable-static --disable-shared
	make -C build/libpng -j 4 libpng16.la

-include $(DEPS)
