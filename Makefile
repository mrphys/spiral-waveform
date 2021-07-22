CXX := g++

SOURCES = $(wildcard src/*.cpp)

CFLAGS = -std=c++98 -fPIC -O2

TARGET_SHARED = build/libspiral_waveform.so
TARGET_STATIC = build/libspiral_waveform.a
TARGET_OBJECT = build/spiral_waveform.o

all: lib lib-static

lib: $(TARGET_SHARED)

$(TARGET_SHARED): $(SOURCES)
	mkdir -p build/
	$(CXX) $(CFLAGS) -o $@ $^ -shared
	
lib-static: $(TARGET_STATIC)

$(TARGET_STATIC): $(SOURCES)
	mkdir -p build/
	$(CXX) -c $(CFLAGS) -o $(TARGET_OBJECT) $^
	ar rvs $@ $(TARGET_OBJECT)

.PHONY: clean
clean:
	rm -rf build/
