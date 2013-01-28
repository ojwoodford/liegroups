CPPFLAGS= -I../
CXXFLAGS=-g -Wall -O2 -finline-functions

OBJS=   se2.o \
	se2_io.o \
	sim2.o \
	sim2_io.o \
	so3.o \
	so3_io.o \
	se3.o \
	se3_io.o \
	sl3.o \
	sl3_io.o \
	matrix.o

all: libliegroups.a
	+ make -C test/

libliegroups.a: $(OBJS) .deps
	$(AR) rvs $@ $^

.deps: *.cpp *.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MG -MP *.cpp > .deps

clean:
	$(RM) $(OBJS) libliegroups.a


-include .deps