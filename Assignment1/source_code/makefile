CCT=gcc-6
CC=gcc
CFLAGS=-w -fopenmp -std=c99 -lm -o
THREADS=-D NUM_THREADS=
MATRIX=-D MATRIX_SIZE=

main:
	$(CC) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)1024

timer:
	$(CC) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c $(THREADS) $(MATRIX)

# This is for my personal machine running 
# gcc version 6 with openMP supprt
test_p:
	$(CCT) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)1024

timer_p4:
	$(CCT) $(CFLAGS) gauss_eliminate1 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)1024
	$(CCT) $(CFLAGS) gauss_eliminate2 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)2048
	$(CCT) $(CFLAGS) gauss_eliminate3 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)4096
	$(CCT) $(CFLAGS) gauss_eliminate4 gauss_eliminate.c compute_gold.c $(THREADS)4 $(MATRIX)8192
	./gauss_eliminate1
	./gauss_eliminate2
	./gauss_eliminate3
	./gauss_eliminate4




