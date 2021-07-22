CXX := g++

SOURCES = $(wildcard src/*.cpp)

CFLAGS = -std=c++98 -fPIC -O2
LDFLAGS = -shared

TARGET = build/spiral_waveform.so

all: lib

lib: $(TARGET)

$(TARGET): $(SOURCES)
	mkdir -p build/
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf build/
