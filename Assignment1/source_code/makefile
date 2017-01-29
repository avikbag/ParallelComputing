CCT=gcc-6
CC=gcc
CFLAGS=-w -fopenmp -std=c99 -lm -o

main:
	$(CC) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c


# This is for my personal machine running 
# gcc version 6 with openMP supprt
test:
	$(CCT) $(CFLAGS) gauss_eliminate gauss_eliminate.c compute_gold.c

