JL_SHARE = $(shell /home/mlubin/foo/julia-bb73f3489d/bin/julia -e 'print(joinpath(JULIA_HOME,Base.DATAROOTDIR,"julia"))')
#CFLAGS   += -g $(shell $(JL_SHARE)/julia-config.jl --cflags) -fPIC -Iinclude
CFLAGS   += -g -DJULIA_INIT_DIR=\"/home/mlubin/foo/julia-bb73f3489d/lib\" -I/home/mlubin/foo/julia-bb73f3489d/include/julia -fPIC -Iinclude
CXXFLAGS += $(shell $(JL_SHARE)/julia-config.jl --cflags)
#LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDFLAGS  += -L/home/mlubin/foo/julia-bb73f3489d/lib/julia
#LDLIBS   += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)
LDLIBS   += -Wl,-rpath,/home/mlubin/foo/julia-bb73f3489d/lib/julia -ljulia

DEPS = include/cmpb.h

all: libcmpb.so

src/cmpb.o: src/cmpb.c
	$(CC) -c ${CFLAGS} -o $@ $^

libcmpb.so: src/cmpb.o
	$(CC) ${LDFLAGS} -shared -o $@ $^ ${LDLIBS}

test_cmpb: test/test_cmpb.c libcmpb.so
	$(CC) ${CFLAGS} ${LDFLAGS} -o $@ test/test_cmpb.c -Wl,-rpath,$(shell pwd) libcmpb.so

clean:
	rm -f libcmpb.so src/cmpb.o test_cmpb
