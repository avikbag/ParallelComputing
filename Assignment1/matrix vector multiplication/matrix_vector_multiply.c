/* Matrix vector multiplication b = Ax using OpenMP. A is nxn matrix and x is a 1xn vector.
 * Author: Naga Kandasamy
 * Date created: 4/25/2011
 * Date last modified: 10/07/2014
 * Compile as follows: gcc -fopenmp matrix_vector_multiply.c -o matrix_vector_multiply -std=c99
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 8
#define NUM_ROWS 8192
#define NUM_COLUMNS 8192

#define THREADING_THRESHOLD 256

// Function prototypes
void compute_gold(float *, float *, float *, int, int);
void compute_using_openmp_v1(float *, float *, float *, int, int);
void compute_using_openmp_v2(float *, float *, float *, int, int);
void compute_using_openmp_v3(float *, float *, float *, int, int);


int 
main(int argc, char **argv)
{
		  
    int num_rows = NUM_ROWS;
    int num_cols = NUM_COLUMNS;

	
    /* Allocate memory. */
    float *matrix_a = (float *)malloc(sizeof(float) * num_rows * num_cols); // Matrix A
    float *vector_x = (float *)malloc(sizeof(float) * num_cols); // Vector x
    float *vector_b = (float *)malloc(sizeof(float) * num_rows); // Vector b
	
    int i, j;
    srand(time(NULL));
    for(i = 0; i < num_rows; i++){       /* Initialize matix A with random nnumbers betwwen [.5, -.5]. */
        for(j = 0; j < num_cols; j++){
            matrix_a[i*num_cols + j] = ((float)rand()/(float)RAND_MAX) - 0.5;
        }
    }
		
    for(i = 0; i < num_cols; i++){  /* Initialize the X vector. */
            vector_x[i] =  ((float)rand()/(float)RAND_MAX) - 0.5;
    }

    /* Multiply using the single-threaded solution. */
    struct timeval start, stop;	
    printf("Vector matrix multiplication using the single-threaded version. \n");
    gettimeofday(&start, NULL);
    
    compute_gold(matrix_a, vector_x, vector_b, num_rows, num_cols); 
    
    gettimeofday(&stop, NULL);
    printf("Execution time = %fs. \n", (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)/(float)1000000));
    printf("\n");
	
    printf("Vector matrix multiplication using OpenMP. Version 1. \n");
    gettimeofday(&start, NULL);
    
    compute_using_openmp_v1(matrix_a, vector_x, vector_b, num_rows, num_cols);
    
    gettimeofday(&stop, NULL);
    printf("Execution time = %fs. \n", (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)/(float)1000000));
    printf("\n");

    printf("Vector matrix multiplication using OpenMP. Version 2. \n");
    gettimeofday(&start, NULL);
	
    compute_using_openmp_v2(matrix_a, vector_x, vector_b, num_rows, num_cols);

    gettimeofday(&stop, NULL);
    printf("Execution time = %fs. \n", (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)/(float)1000000));
    printf("\n");

    printf("Vector matrix multiplication using OpenMP. Version 3. \n");
    gettimeofday(&start, NULL);
	
    compute_using_openmp_v3(matrix_a, vector_x, vector_b, num_rows, num_cols);
	
    gettimeofday(&stop, NULL);
    printf("Execution time = %fs. \n", (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec)/(float)1000000));
    printf("\n");

    /* Free memory here. */ 
    free((void *)matrix_a);
    free((void *)vector_b);
    free((void *)vector_x);
}

/* This function computes the reference soution using a single thread. */
void 
compute_gold(float *matrix_a, float *vector_x, float *vector_b, int num_rows, int num_cols)
{
    int i, j;
    for(i = 0; i < num_rows; i++){
        vector_b[i] = 0;
        for(j = 0; j < num_cols; j++){
            vector_b[i] += matrix_a[i*num_cols + j] * vector_x[j];
        }
    }
}

/* This function computes the solution using OpenMP. */
void 
compute_using_openmp_v1(float *matrix_a, float *vector_x, float *vector_b, int num_rows, int num_cols)
{
    int i, j;

	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel for default(none) shared(matrix_a, vector_b, vector_x, num_rows, num_cols) private(i, j)
	 for(i = 0; i < num_rows; i++){
         vector_b[i] = 0;
         for(j = 0; j < num_cols; j++){
             vector_b[i] += matrix_a[i*num_cols + j] * vector_x[j];
         }
     }
}

/* This function computes the solution using OpenMP. It implements a thresholding scheme for deciding whether to go multi-threaded or not. 
 * This code exhibits some slowdown due to false sharing since the partial sums are being constantly accumulated directly into the 
 * vector_b variable. */
void 
compute_using_openmp_v2(float *matrix_a, float *vector_x, float *vector_b, int num_rows, int num_cols)
{
	int i, j;

	omp_set_num_threads(NUM_THREADS); // Set the number of threads
	#pragma omp parallel for if(num_rows*num_cols > THREADING_THRESHOLD) default(none) shared(matrix_a, vector_b, vector_x, num_rows, num_cols) private(i, j)
	 for(i = 0; i < num_rows; i++){
         vector_b[i] = 0;
         for(j = 0; j < num_cols; j++){
             vector_b[i] += matrix_a[i*num_cols + j] * vector_x[j];
         }
     }
}

/* This function computes the solution using OpenMP. It implements a thresholding scheme for deciding whether to go multi-threaded or not. */
void 
compute_using_openmp_v3(float *matrix_a, float *vector_x, float *vector_b, int num_rows, int num_cols)
{
	int i, j;
	double sum;

	omp_set_num_threads(NUM_THREADS); // Set the number of threads
	#pragma omp parallel for if(num_rows*num_cols > THREADING_THRESHOLD) default(none) shared(matrix_a, vector_b, vector_x, num_rows, num_cols) private(i, j, sum)
	 for(i = 0; i < num_rows; i++){
         sum = 0;
         for(j = 0; j < num_cols; j++){
             sum += matrix_a[i*num_cols + j] * vector_x[j];
         }
         vector_b[i] = (float)sum;
     }
}



