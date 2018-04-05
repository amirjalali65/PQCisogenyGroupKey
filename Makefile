####  Makefile for compilation on Linux  ####

OPT=-O3     # Optimization option by default

CC=gcc
ifeq "$(CC)" "gcc"
    COMPILER=gcc
else ifeq "$(CC)" "clang"
    COMPILER=clang
endif

ARCHITECTURE=_AMD64_
ifeq "$(ARCH)" "x64"
    ARCHITECTURE=_AMD64_
else ifeq "$(ARCH)" "x86"
    ARCHITECTURE=_X86_
else ifeq "$(ARCH)" "ARM"
    ARCHITECTURE=_ARM_
#    ARM_SETTING=-lrt
else ifeq "$(ARCH)" "ARM64"
    ARCHITECTURE=_ARM64_
#    ARM_SETTING=-lrt
endif

ADDITIONAL_SETTINGS=
ifeq "$(SET)" "EXTENDED"
    ADDITIONAL_SETTINGS=-fwrapv -fomit-frame-pointer -march=native
endif

USE_OPT_LEVEL=_OPTIMIZED_GENERIC_

AR=ar rcs
RANLIB=ranlib

CFLAGS=$(OPT) -static $(ADDITIONAL_SETTINGS) -D $(ARCHITECTURE) -D __LINUX__ -D $(USE_OPT_LEVEL)
LDFLAGS=-lm
EXTRA_OBJECTS_747=objs747/fp_generic.o
OBJECTS_747=objs747/P747.o $(EXTRA_OBJECTS_747) objs/random.o 

all: lib747 tests 

objs747/%.o: %.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

objs747/fp_generic.o: generic/fp_generic.c
	$(CC) -c $(CFLAGS) generic/fp_generic.c -o objs747/fp_generic.o

objs/random.o: random.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) random.c -o objs/random.o


lib747: $(OBJECTS_747)
	rm -rf sigk
	mkdir sigk
	$(AR) sigk/libsigk.a $^
	$(RANLIB) sigk/libsigk.a

tests: lib747
	$(CC) $(CFLAGS) -L./sigk tests/test_SIGKp747.c tests/test_extras.c -lsigk $(LDFLAGS) -o sigk/test_SIGK_747 $(ARM_SETTING)

check: tests

.PHONY: clean

clean:
	rm -rf *.req objs747 objs sigk

