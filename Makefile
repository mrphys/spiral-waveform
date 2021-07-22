CXX := g++

HEADERS = $(wildcard src/*.h)
SOURCES = $(wildcard src/*.cpp)

CFLAGS = -std=c++98 -fPIC -O2

TARGET_SHARED = lib/libspiral_waveform.so
TARGET_STATIC = lib/libspiral_waveform.a
TARGET_INCLUDE = include

all: lib lib-static headers

lib: $(TARGET_SHARED)

$(TARGET_SHARED): $(SOURCES)
	mkdir -p lib/
	$(CXX) $(CFLAGS) -o $@ $^ -shared
	
lib-static: $(TARGET_STATIC)

$(TARGET_STATIC): $(SOURCES)
	mkdir -p build/
	mkdir -p lib/
	$(CXX) -c $(CFLAGS) -o build/spiral_waveform.o $^
	ar rvs $@ build/spiral_waveform.o

headers: $(TARGET_INCLUDE) $(HEADERS)

$(TARGET_INCLUDE): $(HEADERS)
	mkdir -p $(TARGET_INCLUDE)
	cp -t $@ $^

.PHONY: clean
clean:
	rm -rf build/
	rm -rf lib/
	rm -rf include/
