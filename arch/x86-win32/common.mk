CROSS=i686-w64-mingw32-
CC=$(CROSS)gcc
CXX=$(CROSS)g++
LD=$(CROSS)gcc
AR=$(CROSS)ar

LDFLAGS=-static-libgcc -static-libstdc++
