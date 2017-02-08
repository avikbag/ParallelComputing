
THREADS=-D NUM_THREADS=
MATRIX=-D MATRIX_SIZE=
main:
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)2048 $(THREADS)8
	./gauss
