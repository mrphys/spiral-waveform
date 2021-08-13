CXX := g++

HEADERS = $(wildcard src/*.h)
SOURCES = $(wildcard src/*.cpp)

CFLAGS = -std=c++98 -fPIC -O3 -march=x86-64 -mtune=generic

TARGET_SHARED = build/libspiral_waveform.so
TARGET_STATIC = build/libspiral_waveform.a
TARGET_INCLUDE = include

INSTALL_PREFIX = /usr/local

all: lib lib-static

lib: $(TARGET_SHARED)

$(TARGET_SHARED): $(SOURCES)
	mkdir -p build/
	$(CXX) $(CFLAGS) -o $@ $^ -shared
	
lib-static: $(TARGET_STATIC)

$(TARGET_STATIC): $(SOURCES)
	mkdir -p build/
	$(CXX) -c $(CFLAGS) -o build/spiral_waveform.o $^
	ar rvs $@ build/spiral_waveform.o

# headers: $(TARGET_INCLUDE) $(HEADERS)

# $(TARGET_INCLUDE): $(HEADERS)
# 	mkdir -p $(TARGET_INCLUDE)
# 	cp -t $@ $^

install: $(TARGET_SHARED) $(TARGET_STATIC)
	cp src/spiral_waveform.h $(INSTALL_PREFIX)/include/
	cp build/libspiral_waveform.so $(INSTALL_PREFIX)/lib
	cp build/libspiral_waveform.a $(INSTALL_PREFIX)/lib

.PHONY: clean
clean:
	rm -rf build/
