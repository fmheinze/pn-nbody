# Portable Makefile for pn-nbody
#
# OpenMP control:
#   make                    # auto-detect OpenMP (default)
#   make USE_OPENMP=1       # require OpenMP
#   make USE_OPENMP=0       # force a serial build
#
# Compiler examples:
#   make CC=gcc
#   make CC=clang
#   make CC=gcc-15
#   make CC="$(brew --prefix llvm)/bin/clang"   # macOS/Homebrew LLVM

SHELL := /bin/sh

# Compiler
CC ?= cc

# Directories
SRC_DIR   := src
BUILD_DIR := build
EXE_DIR   := exe

# Target
TARGET := $(EXE_DIR)/pn-nbody

# Sources / objects
SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# User-overridable flags.
# Project-required and OpenMP flags are appended separately, so a command such as
#   make CFLAGS="-O0 -g" USE_OPENMP=1
# still receives warning/dependency/OpenMP flags.
CPPFLAGS ?=
CFLAGS   ?= -O3 -march=native -mtune=native -fno-math-errno -DNDEBUG
LDFLAGS  ?=
LDLIBS   ?=

PROJECT_CFLAGS := -Wall -Wextra -MMD -MP
PROJECT_LDLIBS := -lm

# -----------------------------------------------------------------------------
# OpenMP detection
# -----------------------------------------------------------------------------
# USE_OPENMP=auto (default): enable OpenMP if a compile+link probe succeeds.
# USE_OPENMP=1:              require OpenMP; fail if it cannot be configured.
# USE_OPENMP=0:              disable OpenMP explicitly.
USE_OPENMP ?= auto

UNAME_S    := $(shell uname -s 2>/dev/null)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)

OMP_MODE     := disabled
OMP_CPPFLAGS :=
OMP_CFLAGS   :=
OMP_LDFLAGS  :=
OMP_LDLIBS   :=

# Validate USE_OPENMP before probing.
ifneq ($(USE_OPENMP),auto)
ifneq ($(USE_OPENMP),0)
ifneq ($(USE_OPENMP),1)
$(error USE_OPENMP must be auto, 0, or 1)
endif
endif
endif

# Only probe when OpenMP has not been explicitly disabled.
ifneq ($(USE_OPENMP),0)

# Native compiler OpenMP support.  The test uses '-include omp.h' rather than a
# literal preprocessor line inside $(shell ...); this is compatible with the
# older GNU Make shipped by macOS as well as current GNU Make releases.
OMP_NATIVE_OK := $(shell tmp=/tmp/pn-nbody-omp-native-$$$$; printf '%s\n' 'int main(void) { return omp_get_max_threads() < 1; }' | $(CC) -include omp.h -fopenmp -x c - -o $$tmp -fopenmp >/dev/null 2>&1; rc=$$?; rm -f $$tmp; if test $$rc -eq 0; then echo 1; else echo 0; fi)

ifeq ($(OMP_NATIVE_OK),1)
OMP_MODE    := native
OMP_CFLAGS  := -fopenmp
OMP_LDFLAGS := -fopenmp
else

# Apple Clang normally needs Homebrew's separate libomp package.
ifeq ($(UNAME_S),Darwin)
LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
ifneq ($(strip $(LIBOMP_PREFIX)),)
OMP_BREW_OK := $(shell tmp=/tmp/pn-nbody-omp-brew-$$$$; printf '%s\n' 'int main(void) { return omp_get_max_threads() < 1; }' | $(CC) -I$(LIBOMP_PREFIX)/include -Xpreprocessor -fopenmp -include omp.h -x c - -o $$tmp -L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib -lomp >/dev/null 2>&1; rc=$$?; rm -f $$tmp; if test $$rc -eq 0; then echo 1; else echo 0; fi)
ifeq ($(OMP_BREW_OK),1)
OMP_MODE     := homebrew-libomp
OMP_CPPFLAGS := -I$(LIBOMP_PREFIX)/include
OMP_CFLAGS   := -Xpreprocessor -fopenmp
OMP_LDFLAGS  := -L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib
OMP_LDLIBS   := -lomp
endif
endif
else

# Optional Clang/libomp fallback on Linux and other Unix-like systems when
# pkg-config exposes libomp but plain '-fopenmp' does not work.
LIBOMP_PKG_CFLAGS := $(shell pkg-config --cflags libomp 2>/dev/null)
LIBOMP_PKG_LIBS   := $(shell pkg-config --libs libomp 2>/dev/null)
ifneq ($(strip $(LIBOMP_PKG_LIBS)),)
OMP_PKG_OK := $(shell tmp=/tmp/pn-nbody-omp-pkg-$$$$; printf '%s\n' 'int main(void) { return omp_get_max_threads() < 1; }' | $(CC) -fopenmp=libomp $(LIBOMP_PKG_CFLAGS) -include omp.h -x c - -o $$tmp $(LIBOMP_PKG_LIBS) >/dev/null 2>&1; rc=$$?; rm -f $$tmp; if test $$rc -eq 0; then echo 1; else echo 0; fi)
ifeq ($(OMP_PKG_OK),1)
OMP_MODE     := pkg-config-libomp
OMP_CPPFLAGS := $(LIBOMP_PKG_CFLAGS)
OMP_CFLAGS   := -fopenmp=libomp
OMP_LDLIBS   := $(LIBOMP_PKG_LIBS)
endif
endif
endif
endif

# If OpenMP was explicitly requested, a failed probe is an error.  In auto mode
# the build simply remains serial.
ifeq ($(USE_OPENMP),1)
ifeq ($(OMP_MODE),disabled)
ifeq ($(UNAME_S),Darwin)
$(error OpenMP was requested but no working setup was found. With Apple Clang, install Homebrew libomp with 'brew install libomp', or select an OpenMP-capable compiler)
else
$(error OpenMP was requested but no working setup was found. Install an OpenMP-capable compiler/runtime, or build with USE_OPENMP=0)
endif
endif
endif

endif

# Effective flags
ALL_CPPFLAGS := $(CPPFLAGS) $(OMP_CPPFLAGS)
ALL_CFLAGS   := $(CFLAGS) $(PROJECT_CFLAGS) $(OMP_CFLAGS)
ALL_LDFLAGS  := $(LDFLAGS) $(OMP_LDFLAGS)
ALL_LDLIBS   := $(LDLIBS) $(PROJECT_LDLIBS) $(OMP_LDLIBS)

# Build information
$(info Compiler: $(CC))
$(info Compiler version: $(CC_VERSION))
$(info OpenMP: $(OMP_MODE))

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS) | $(EXE_DIR)
	$(CC) $(OBJS) -o $@ $(ALL_LDFLAGS) $(ALL_LDLIBS)

# Compile
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

# Directories
$(BUILD_DIR) $(EXE_DIR):
	mkdir -p $@

# Show the resolved configuration without compiling.
show-config:
	@echo "CC          = $(CC)"
	@echo "OpenMP      = $(OMP_MODE)"
	@echo "CPPFLAGS    = $(ALL_CPPFLAGS)"
	@echo "CFLAGS      = $(ALL_CFLAGS)"
	@echo "LDFLAGS     = $(ALL_LDFLAGS)"
	@echo "LDLIBS      = $(ALL_LDLIBS)"

clean:
	rm -rf $(BUILD_DIR) $(EXE_DIR)

-include $(DEPS)

.PHONY: all clean show-config
