PROJNAME = gauss_eliminate
VIEWER = less
reverse = $(if $(wordlist 2,2,$(1)),$(call reverse,$(wordlist 2,$(words $(1)),$(1))) $(firstword $(1)),$(1))

IDIR = include
ODIR = build
SRCDIR = src

BUILD_ENV_BASE = $(ODIR)
BUILD_ENV =  $(BUILD_ENV_BASE)

DEPF = gauss_eliminate.h 
DEPS = $(patsubst %,$(IDIR)/%,$(DEPF))

OBJF =  main.o
OBJF += compute_gold.o
OBJF += compute_threaded.o
OBJ = $(patsubst %,$(ODIR)/%,$(OBJF))

DRIVER = $(SRCDIR)/main.c

NUM_THREADS=8
MATRIX_SIZE=1024

CC = gcc
CCOPTS = -g -Wall -Wextra -lm -fopenmp -Wno-unused-variable -Wno-missing-field-initializers -I$(IDIR)

.PHONY: all view build env clean rebuild preprocessor test
.DEFAULT: build

build : | env all

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	- $(CC) -c -o $@ $< $(CCOPTS) -DARG_THREADS=$(NUM_THREADS) -DARG_MATRIX_SIZE=$(MATRIX_SIZE)

all : $(OBJ)
	- $(CC) -o $(PROJNAME) $^ $(CCOPTS)


env: 
	-@ mkdir -p $(BUILD_ENV)

clean :
	-@ \rm -f $(OBJ) $(SRCDIR)/*~ $(T_SRCDIR)/*~ core $(INCDIR)/*~ $(T_INCDIR)/*~ $(PROJNAME)
	-@ \rmdir $(call reverse, $(BUILD_ENV)) 

rebuild : | clean build

view :
	-@ $(VIEWER) $(DRIVER) $(DEPS)

run :
	-@ ./$(PROJNAME)

preprocessor :
	-@ $(CC) -c $(DRIVER) -E $(CCOPTS)

test :
	-@ make rebuild NUM_THREADS=$(NUM_THREADS) MATRIX_SIZE=$(MATRIX_SIZE)
	-@ ./$(PROJNAME)


