CXX ?= g++
CXXFLAGS ?= -g -std=c++17 -Wall -Werror -I../include $(EXTRA_FLAGS)

ifeq ($(OS),Windows_NT)
    DETECTED_OS := Windows
else
    DETECTED_OS := $(shell uname -s)
endif

ifeq ($(DETECTED_OS),Darwin)
	CXXFLAGS += -mmacosx-version-min=13.3
endif

test_doubledouble: test_doubledouble.cpp checkit.h ../include/doubledouble.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -rf test_doubledouble
