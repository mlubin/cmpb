JL_SHARE = $(shell julia -e 'print(joinpath(JULIA_HOME,Base.DATAROOTDIR,"julia"))')
CFLAGS   += -g $(shell $(JL_SHARE)/julia-config.jl --cflags) -fPIC -Iinclude
CXXFLAGS += $(shell $(JL_SHARE)/julia-config.jl --cflags)
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDLIBS   += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)

DEPS = include/cmpb.h
OS := $(shell uname)

ifeq ($(OS),Darwin)
	LIBCMPB := libcmpb.dylib
else
	LIBCMPB := libcmpb.so
endif

all: $(LIBCMPB)

src/cmpb.o: src/cmpb.c
	$(CC) -c ${CFLAGS} -o $@ $^

$(LIBCMPB): src/cmpb.o
	$(CC) ${LDFLAGS} -shared -o $@ $^ ${LDLIBS}

test/test_cmpb.o: test/test_cmpb.c
	$(CC) -c ${CFLAGS} -o $@ test/test_cmpb.c

test_cmpb: test/test_cmpb.o $(LIBCMPB)
	$(CC) ${LDFLAGS} test/test_cmpb.o src/cmpb.o ${LDLIBS} -o $@

clean:
	rm -f $(LIBCMPB) src/cmpb.o test_cmpb test/test_cmpb.o
