COMPILER=gcc-6
CFLAGS=-std=c99 -fopenmp -lm -o
ARGS=1000000
all:
	$(COMPILER) $(CFLAGS) histogram histogram.c 
	./histogram $(ARGS)
