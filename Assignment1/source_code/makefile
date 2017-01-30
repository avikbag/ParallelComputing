CCT=gcc-6
CC=gcc
CFLAGS=-w -fopenmp -std=c99 -lm -O3 -o
THREADS=-D NUM_THREADS=
MATRIX=-D MATRIX_SIZE=

main:
	$(CC) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)8192

# This is for my personal machine running 
# gcc version 6 with openMP support
test_p:
	$(CCT) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c $(THREADS)16 $(MATRIX)8192

## This is for running tests on the various threads vs matrices sizes
timer_test: timer_p2 timer_p4 timer_p8 timer_p16

timer_p2:
	$(CCT) $(CFLAGS) gauss_eliminate21 gauss_eliminate.c compute_gold.c $(THREADS)2 $(MATRIX)1024
	$(CCT) $(CFLAGS) gauss_eliminate22 gauss_eliminate.c compute_gold.c $(THREADS)2 $(MATRIX)2048
	$(CCT) $(CFLAGS) gauss_eliminate23 gauss_eliminate.c compute_gold.c $(THREADS)2 $(MATRIX)4096
	$(CCT) $(CFLAGS) gauss_eliminate24 gauss_eliminate.c compute_gold.c $(THREADS)2 $(MATRIX)8192
	./gauss_eliminate21 && ./gauss_eliminate22 && ./gauss_eliminate23 && ./gauss_eliminate24

timer_p4:
	$(CCT) $(CFLAGS) gauss_eliminate41 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)1024
	$(CCT) $(CFLAGS) gauss_eliminate42 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)2048
	$(CCT) $(CFLAGS) gauss_eliminate43 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)4096
	$(CCT) $(CFLAGS) gauss_eliminate44 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)8192
	./gauss_eliminate41 && ./gauss_eliminate42 && ./gauss_eliminate43 && ./gauss_eliminate44 

timer_p8:
	$(CCT) $(CFLAGS) gauss_eliminate81 gauss_eliminate.c compute_gold.c $(THREADS)8 $(MATRIX)1024
	$(CCT) $(CFLAGS) gauss_eliminate82 gauss_eliminate.c compute_gold.c $(THREADS)8 $(MATRIX)2048
	$(CCT) $(CFLAGS) gauss_eliminate83 gauss_eliminate.c compute_gold.c $(THREADS)8 $(MATRIX)4096
	$(CCT) $(CFLAGS) gauss_eliminate84 gauss_eliminate.c compute_gold.c $(THREADS)8 $(MATRIX)8192
	./gauss_eliminate81 && ./gauss_eliminate82 && ./gauss_eliminate83 && ./gauss_eliminate84

timer_p16:
	$(CCT) $(CFLAGS) gauss_eliminate161 gauss_eliminate.c compute_gold.c $(THREADS)16 $(MATRIX)1024
	$(CCT) $(CFLAGS) gauss_eliminate162 gauss_eliminate.c compute_gold.c $(THREADS)16 $(MATRIX)2048
	$(CCT) $(CFLAGS) gauss_eliminate163 gauss_eliminate.c compute_gold.c $(THREADS)16 $(MATRIX)4096
	$(CCT) $(CFLAGS) gauss_eliminate164 gauss_eliminate.c compute_gold.c $(THREADS)16 $(MATRIX)8192
	./gauss_eliminate161 && ./gauss_eliminate162 && ./gauss_eliminate163 && ./gauss_eliminate164

clean:
	rm gauss_eliminate21 gauss_eliminate22 gauss_eliminate23 gauss_eliminate24 gauss_eliminate41 gauss_eliminate42 gauss_eliminate43 gauss_eliminate44 gauss_eliminate81 gauss_eliminate82 gauss_eliminate83 gauss_eliminate84 gauss_eliminate161 gauss_eliminate162 gauss_eliminate163 gauss_eliminate164 gauss_eliminate
