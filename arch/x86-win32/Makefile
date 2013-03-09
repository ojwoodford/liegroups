include common.mk

CPPFLAGS= -I../../../
CXXFLAGS=-g -Wall -O2 -finline-functions

include ../../common.mk

%.o: ../../%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

all: libliegroups.a
	+ make -C test/ all

libliegroups.a: $(OBJS)
	$(AR) rvs $@ $^

.deps: ../../*.cpp ../../*.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MG -MP $^ > .deps

clean:
	$(RM) $(OBJS) libliegroups.a
	+ make -C test/ clean

-include .deps
