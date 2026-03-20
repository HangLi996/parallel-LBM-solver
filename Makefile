CC := mpicxx
# CC := mpiicpc  # For Intel MPI with Intel compiler (requires Intel compiler module)
# CC := mpic++  # For OpenMPI
# CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR    := src
SPECDIR   := specs
BUILDDIR  := build
TARGETDIR := bin
TARGET    := main

# Compile all .cpp files in the src directory and put the object files in the build dir
SRCEXT  := cpp

# All except main.cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT) -not -wholename src/main.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# All spec .cpp files
SPECS   := $(shell find $(SPECDIR) -type f -name *.$(SRCEXT))
SPEC_OBJECTS := $(patsubst $(SPECDIR)/%,$(BUILDDIR)/$(SPECDIR)/%,$(SPECS:.$(SRCEXT)=.o))

#OpenMP flags
OMPFLAGS := -fopenmp

# MPI flags
# CFLAGS  := -std=c++14 -g -MMD $(OMPFLAGS)
CFLAGS  := -std=c++14 -O3 -march=native -fopt-info-vec-all=vec_report.txt -g -MMD $(OMPFLAGS)
# For domain decomposition strategies, use:
# -DDECOMP_VERTICAL (default, no flag needed)
# -DDECOMP_HORIZONTAL for horizontal decomposition
# -DDECOMP_2D for 2D block decomposition
# Example: CFLAGS += -DDECOMP_HORIZONTAL

LIB     := -lmpi $(OMPFLAGS)
INC     := -I include

all: lbm

$(TARGET): $(OBJECTS)
	@echo "Linking..."
	$(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/specs/%.o: $(SPECDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -I lib/bandit -c -o $@ $<

lbm: $(OBJECTS) build/main.o
	@mkdir -p $(TARGETDIR)
	$(CC) $^ -o $(TARGETDIR)/$(TARGET) $(LIB)

# Build with horizontal decomposition
lbm_horizontal: CFLAGS += -DDECOMP_HORIZONTAL
lbm_horizontal: lbm

# Build with 2D block decomposition
lbm_2d: CFLAGS += -DDECOMP_2D
lbm_2d: lbm

# Tests
spec: $(SPEC_OBJECTS) $(OBJECTS)
	@echo "Linking spec runner..."
	$(CC) $(CFLAGS) $^ -o $(TARGETDIR)/spec $(LIB)

clean:
	@echo "Cleaning...";
	$(RM) -r $(BUILDDIR) $(TARGETDIR)/$(TARGET)

.PHONY: clean lbm lbm_horizontal lbm_2d

