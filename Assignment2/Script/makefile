COMPILERH=gcc-6
COMPILER=gcc
CFLAGS=-std=c99 -fopenmp -lm -o
ARGS=100000000
THREADS=-D NUM_THREADS=
main:
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram $(ARGS)

## Test case on xunil
test: test2 test4 test8 test16

test2:
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)2
	./histogram 1000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)2
	./histogram 10000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)2
	./histogram 100000000

test4:
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)4
	./histogram 1000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)4
	./histogram 10000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)4
	./histogram 100000000

test8:
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)8
	./histogram 1000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)8
	./histogram 10000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)8
	./histogram 100000000

test16:
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram 1000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram 10000000
	$(COMPILER) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram 100000000

## Test case on my personal machine using gcc-6 with OpenMP
testH: test2H test4H test8H test16H

test2H:
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)2
	./histogram 1000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)2
	./histogram 10000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)2
	./histogram 100000000

test4H:
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)4
	./histogram 1000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)4
	./histogram 10000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)4
	./histogram 100000000

test8H:
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)8
	./histogram 1000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)8
	./histogram 10000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)8
	./histogram 100000000

test16H:
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram 1000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram 10000000
	$(COMPILERH) $(CFLAGS) histogram histogram.c $(THREADS)16
	./histogram 100000000
