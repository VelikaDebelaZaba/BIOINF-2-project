
---
Pretpostavke:
- `.cpp` datoteke su u `src/`
- `.hpp` datoteke su u `include/`
- build ide u `build/`
- binarka ide u `bin/`

```makefile
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra
INCLUDES := -Iinclude

SRC := $(shell find src -name "*.cpp")
OBJ := $(SRC:src/%.cpp=build/%.o)

TARGET := bin/launcher

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build bin

.PHONY: all clean
