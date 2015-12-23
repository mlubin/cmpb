JL_SHARE = $(shell julia -e 'print(joinpath(JULIA_HOME,Base.DATAROOTDIR,"julia"))')
CFLAGS   += -g $(shell $(JL_SHARE)/julia-config.jl --cflags) -fPIC -Iinclude
CXXFLAGS += $(shell $(JL_SHARE)/julia-config.jl --cflags)
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDLIBS   += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)

DEPS = include/cmpb.h

all: libcmpb.so

src/cmpb.o: src/cmpb.c
	$(CC) -c ${CFLAGS} -o $@ $^

libcmpb.so: src/cmpb.o
	$(CC) ${LDFLAGS} -shared -o $@ $^ ${LDLIBS}

test_cmpb: test/test_cmpb.c libcmpb.so
	$(CXX) ${CFLAGS} ${LDFLAGS} -o $@ test/test_cmpb.c libcmpb.so

clean:
	rm -f libcmpb.so src/cmpb.o test_cmpb
