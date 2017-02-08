
THREADS=-D NUM_THREADS=
MATRIX=-D MATRIX_SIZE=

main:
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)2048 $(THREADS)8
	./gauss

test: test2 test4 test8 test16

test2:
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)1024 $(THREADS)2
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)2048 $(THREADS)2
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)4096 $(THREADS)2
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)8192 $(THREADS)2
	./gauss

test4:
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)1024 $(THREADS)4
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)2048 $(THREADS)4
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)4096 $(THREADS)4
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)8192 $(THREADS)4
	./gauss

test8:
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)1024 $(THREADS)8
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)2048 $(THREADS)8
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)4096 $(THREADS)8
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)8192 $(THREADS)8
	./gauss

test16:
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)1024 $(THREADS)16
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)2048 $(THREADS)16
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)4096 $(THREADS)16
	./gauss
	gcc -o gauss gauss_eliminate.c compute_gold.c -lpthread -std=c99 -w -lm $(MATRIX)8192 $(THREADS)16
	./gauss
