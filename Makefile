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

test/test_cmpb.o: test/test_cmpb.c
	$(CC) -c ${CFLAGS} -o $@ test/test_cmpb.c

test_cmpb: test/test_cmpb.o libcmpb.so
	$(CXX) -fPIC -std=c++11 ${LDFLAGS} test/test_cmpb.o src/cmpb.o ${LDLIBS} -o $@

clean:
	rm -f libcmpb.so src/cmpb.o test_cmpb test/test_cmpb.o
