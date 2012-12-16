CPPFLAGS= -I../
CXXFLAGS=-g -Wall -O2 -finline-functions

OBJS=   se2.o \
	se2_io.o \
	sim2.o \
	sim2_io.o \
	matrix.o


libliegroups.a: $(OBJS) .deps
	$(AR) rvs $@ $^

.deps: *.cpp *.hpp
	rm -f .deps
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MG -MP *.cpp > .deps

clean:
	$(RM) $(OBJS) libliegroups.a


-include .deps