CXX := g++

HEADERS = $(wildcard src/*.h)
SOURCES = $(wildcard src/*.cpp) \
          $(wildcard ../RewinderWaveform/*.cpp)

CFLAGS = -std=c++98 -fPIC -O3 -march=x86-64 -mtune=generic
CFLAGS += -I../RewinderWaveform

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

OBJECTS = $(SOURCES:%.cpp=%.o)

$(TARGET_STATIC): $(OBJECTS)
	ar rvs $@ $(OBJECTS)

%.o: %.cpp
	$(CXX) -c $(CFLAGS) -I../RewinderWaveform -o $@ $<

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
