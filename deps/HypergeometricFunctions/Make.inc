LIB = hypergeometricfunctions
LIBDIR = .
LIBFLAGS = -shared -fPIC

ifeq ($(HF_DETECT_GCC), 1)
    _GCCMAJOR := $(shell gfortran -dumpversion | cut -d'.' -f1)
    _GCCMINOR := $(shell gfortran -dumpversion | cut -d'.' -f2)
    ifeq ($(_GCCMAJOR),4)
      _GCCVER := $(_GCCMAJOR).$(_GCCMINOR)
    else
      _GCCVER := $(_GCCMAJOR)
    endif
    CC := gcc-$(_GCCVER)
    COV := gcov-$(_GCCVER)
endif

ifndef target
    target := $(shell $(CC) -dumpmachine)
endif

ifneq (, $(findstring darwin, $(target)))
    UNAME := Darwin
else ifneq (, $(findstring cygwin, $(target)))
    UNAME := Windows
else ifneq (, $(findstring mingw, $(target)))
    UNAME := Windows
else ifneq (, $(findstring msys, $(target)))
    UNAME := Windows
else ifneq (, $(findstring linux, $(target)))
    UNAME := Linux
else ifneq (, $(findstring freebsd, $(target)))
    UNAME := FreeBSD
else
    UNAME := Unknown
endif

simd := $(shell $(CC) -march=native -E -v - </dev/null 2>&1 | grep cc1)

ifneq (, $(filter -msse +sse, $(simd)))
    MSSE := -msse
endif
ifneq (, $(filter -msse2 +sse2, $(simd)))
    MSSE2 := -msse2
endif
ifneq (, $(filter -mavx +avx, $(simd)))
    MAVX := -mavx
endif
ifneq (, $(filter -mfma +fma, $(simd)))
    ifeq (, $(filter -mfma4 +fma4, $(simd)))
        MFMA := -mfma
    endif
endif
ifneq (, $(filter -mavx512f +avx512f -march=knl -march=knm -march=skylake-avx512 -march=cannonlake -march=icelake -march=cascadelake -march=cooperlake -march=tigerlake -march=sapphirerapids, $(simd)))
    ifneq ($(UNAME), Windows)
        MAVX512F := -mavx512f
    endif
endif

ifeq ($(UNAME), Darwin)
    SLIB = dylib
else ifeq ($(UNAME), Windows)
    SLIB = dll
else
    SLIB = so
endif

SRC = src/hypergeometricfunctions.c

machine := $(shell $(CC) -dumpmachine | cut -d'-' -f1)

AFLAGS += -O3 -fPIC -std=gnu11 -I./src

ifndef CFLAGS
    CFLAGS = -O3
    ifneq (, $(findstring 86, $(machine)))
        CFLAGS += -march=native
    endif
    ifneq (, $(findstring aarch64, $(machine)))
        ifeq (, $(findstring clang, $(shell $(CC) --version)))
            CFLAGS += -mcpu=native
        endif
    endif
endif
CFLAGS += -std=gnu11 -I./src

ifeq ($(HF_COVERAGE), 1)
    CFLAGS += -coverage
endif

ifdef HF_PREFIX
    CFLAGS += -I$(HF_PREFIX)/include
    ifeq ($(UNAME), Windows)
        LDFLAGS += -L$(HF_PREFIX)/bin
    else
        LDFLAGS += -L$(HF_PREFIX)/lib
    endif
else
    ifeq ($(UNAME), Darwin)
        CFLAGS += -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
        LDFLAGS += -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
        CFLAGS += -I/usr/local/opt/mpfr/include
        LDFLAGS += -L/usr/local/opt/mpfr/lib
        CFLAGS += -I/usr/local/opt/gmp/include
        LDFLAGS += -L/usr/local/opt/gmp/lib
    else ifeq ($(UNAME), Windows)
        CFLAGS += -Ic:\tools\vcpkg\installed\x64-windows\include\openblas
        CFLAGS += -Ic:\tools\vcpkg\installed\x64-windows\include
        LDFLAGS += -Lc:\tools\vcpkg\installed\x64-windows\bin
    endif
endif


LDLIBS += -lm -lmpfr

ifdef HF_QUADMATH
    AFLAGS += -DHF_QUADMATH
    AFLAGS += -lquadmath
    CFLAGS += -DHF_QUADMATH
    LDLIBS += -lquadmath
endif
