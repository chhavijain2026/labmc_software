# Portable Makefile for labmc
#
# Usage examples:
#   make                 # MPI/parallel build (default)
#   make MPI=0           # serial build
#   make CXX=g++ MPI=0   # serial build with a specified compiler
#   make CXX=clang++ MPI=0
#   make DEBUG=1         # debug build
#   make clean
#
# Notes:
# - Modified (using AI) the original Makefile so no platform-specific compiler 
#   paths are hard-coded.
# - Users can override compiler and flags from the command line.
# - On Windows, this Makefile is intended for WSL, MSYS2/MinGW, or similar
#   Unix-like build environments.

# -----------------------------------------------------------------------------
# Build options
# -----------------------------------------------------------------------------

MPI ?= 1
DEBUG ?= 0

# The executable can be renamed with, for example, make TARGET=my_program.
# We intentionally ignore environment variables named TARGET, because some
# systems define TARGET for unrelated purposes.
ifeq ($(origin TARGET), command line)
    TARGET_BASE := $(TARGET)
else
    TARGET_BASE := labmc
endif

# Add .exe automatically when building under native Windows make environments.
ifeq ($(OS),Windows_NT)
    EXEEXT := .exe
else
    EXEEXT :=
endif

TARGET := $(TARGET_BASE)$(EXEEXT)

# -----------------------------------------------------------------------------
# Compiler selection
# -----------------------------------------------------------------------------

ifeq ($(MPI),1)
# for parallel
    ifneq ($(origin CXX), command line)
        CXX := mpic++
    endif
    SRC_LABMC = main_parallel.cc parameter.cc flowlaw_new.cc state.cc statmodel.cc \
                sm_mnbrak.cc sm_brent.cc gaussj.cc ran2.cc gasdev.cc util.cc error.cc
else
# for serial
    ifneq ($(origin CXX), command line)
        CXX := g++
    endif
    SRC_LABMC = main_serial.cc parameter.cc flowlaw_new.cc state.cc statmodel.cc \
                sm_mnbrak.cc sm_brent.cc gaussj.cc ran2.cc gasdev.cc util.cc error.cc
endif

# -----------------------------------------------------------------------------
# Flags
# -----------------------------------------------------------------------------

# Use a conservative C++ standard. Change to -std=c++11 or newer if needed.
CXXSTD ?= -std=c++11 #c++98

# User-overridable include paths and libraries, e.g.
#   make INCLUDES="-I/path/to/include" LDLIBS="-L/path/to/lib -lname"
INCLUDES ?= -I.
LDLIBS ?= -lm

CXXFLAGS ?= $(CXXSTD) -Wall -Wextra

ifeq ($(DEBUG),1)
    CXXFLAGS += -O0 -g
else
    CXXFLAGS += -ansi -O3
endif

OBJ_LABMC = $(SRC_LABMC:.cc=.o)

# -----------------------------------------------------------------------------
# Build rules
# -----------------------------------------------------------------------------

all: $(TARGET)

$(TARGET): $(OBJ_LABMC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ_LABMC) labmc labmc.exe

print-config:
	@echo "MPI       = $(MPI)"
	@echo "CXX       = $(CXX)"
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "INCLUDES  = $(INCLUDES)"
	@echo "LDLIBS    = $(LDLIBS)"
	@echo "TARGET    = $(TARGET)"
	@echo "SOURCES   = $(SRC_LABMC)"

.PHONY: all clean print-config
