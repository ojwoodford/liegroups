include ../common.mk

CPPFLAGS=-I../../../../
CXXFLAGS=-g -Wall -O2 -finline-functions

LDFLAGS+= -L../
LOADLIBES= -lliegroups

PROGS= main expm extrinsics moments se3_moments se3_gen

%.o: ../../../test/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

%: %.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(LOADLIBES)


vpath %.a ../

all: $(PROGS)

$(PROGS): -lliegroups

clean:
	$(RM) *.o $(PROGS)


