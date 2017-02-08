
THREADS=-D NUM_THREADS=
MATRIX=-D MATRIX_SIZE=
main:
	gcc-6 -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -lm $(MATRIX)2048 $(THREADS)8
	./gauss
