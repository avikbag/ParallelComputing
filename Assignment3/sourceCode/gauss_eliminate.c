/* Gaussian elimination code.
 * Author: Naga Kandasamy
 * Date created: 02/07/2014
 * Date of last update: 01/30/2017
 * Compile as follows: gcc -o gauss_eliminate gauss_eliminate.c compute_gold.c -lpthread -std=c99 -lm
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "gauss_eliminate.h"

#define MIN_NUMBER 2
#define MAX_NUMBER 50
#define NUM_THREADS 4

/* Struct for passing to parallel function*/
struct Args{
  unsigned int start; 
  unsigned int stop;
  float *src;
  unsigned int num_elements;
};

void *gauss_parallel(void *arg);
void *gauss_reduction(void *arg);
void *gauss_elimination(void *arg);

/* Function prototypes. */
extern int compute_gold (float *, unsigned int);
Matrix allocate_matrix (int num_rows, int num_columns, int init);
void gauss_eliminate_using_pthreads (Matrix);
int perform_simple_check (const Matrix);
void print_matrix (const Matrix);
float get_random_number (int, int);
int check_results (float *, float *, unsigned int, float);


int
main (int argc, char **argv)
{
  /* Check command line arguments. */
  if (argc > 1)
    {
      printf ("Error. This program accepts no arguments. \n");
      exit (0);
    }

  /* Matrices for the program. */
  Matrix A;			// The input matrix
  Matrix U_reference;		// The upper triangular matrix computed by the reference code
  Matrix U_mt;			// The upper triangular matrix computed by the pthread code

  /* Initialize the random number generator with a seed value. */
  srand (time (NULL));

  /* Allocate memory and initialize the matrices. */
  A = allocate_matrix (MATRIX_SIZE, MATRIX_SIZE, 1);	// Allocate and populate a random square matrix
  U_reference = allocate_matrix (MATRIX_SIZE, MATRIX_SIZE, 0);	// Allocate space for the reference result
  U_mt = allocate_matrix (MATRIX_SIZE, MATRIX_SIZE, 0);	// Allocate space for the multi-threaded result

  /* Copy the contents of the A matrix into the U matrices. */
  for (int i = 0; i < A.num_rows; i++)
  {
    for (int j = 0; j < A.num_rows; j++)
	  {
	    U_reference.elements[A.num_rows * i + j] = A.elements[A.num_rows * i + j];
	    U_mt.elements[A.num_rows * i + j] = A.elements[A.num_rows * i + j];
	  }
  }

  printf ("\nPerforming gaussian elimination using the reference code. \n");
  struct timeval start, stop;
  gettimeofday (&start, NULL);
  int status = compute_gold (U_reference.elements, A.num_rows);
  gettimeofday (&stop, NULL);
  printf ("CPU run time = %0.2f s. \n",
	  (float) (stop.tv_sec - start.tv_sec +
		   (stop.tv_usec - start.tv_usec) / (float) 1000000));

  if (status == 0)
    {
      printf("Failed to convert given matrix to upper triangular. Try again. Exiting. \n");
      exit (0);
    }
  status = perform_simple_check (U_reference);	// Check that the principal diagonal elements are 1 
  if (status == 0)
    {
      printf ("The upper triangular matrix is incorrect. Exiting. \n");
      exit (0);
    }
  printf ("Single-threaded Gaussian elimination was successful. \n");

  /* Perform the Gaussian elimination using pthreads. The resulting upper triangular matrix should be returned in U_mt */
  printf ("\nPerforming gaussian elimination using the pthread implementation. \n");
  gettimeofday (&start, NULL);
  gauss_eliminate_using_pthreads (U_mt);
  gettimeofday (&stop, NULL);
  printf ("CPU run time = %0.2f s. \n",
	  (float) (stop.tv_sec - start.tv_sec +
		   (stop.tv_usec - start.tv_usec) / (float) 1000000));

  status = perform_simple_check (U_mt);	// Check that the principal diagonal elements are 1 
  if (status == 0)
  {
    printf ("The upper triangular matrix is incorrect. Exiting. \n");
    exit (0);
  }
  printf ("Multi-threaded Gaussian elimination was successful. \n");

  /* check if the pthread result is equivalent to the expected solution within a specified tolerance. */
  int size = MATRIX_SIZE * MATRIX_SIZE;
  int res = check_results (U_reference.elements, U_mt.elements, size, 0.001f);
  printf ("\nTest %s\n", (1 == res) ? "PASSED" : "FAILED");

  /* Free memory allocated for the matrices. */
  free (A.elements);
  free (U_reference.elements);
  free (U_mt.elements);

  return 0;
}


/* Write code to perform gaussian elimination using pthreads. */
void
gauss_eliminate_using_pthreads (Matrix A)
{
  struct Args *args = malloc(sizeof(struct Args)); 
  unsigned int start, stop;
  unsigned int step;
  
  pthread_t *thread_handles = malloc(NUM_THREADS * sizeof(pthread_t));

  step = MATRIX_SIZE / NUM_THREADS;
  for(int i = 0; i < NUM_THREADS; i++)
  {
    start = i * step;
    stop = start + step;
    
    args->start = start;
    args->stop = stop;
    args->src = A.elements;
    args->num_elements = A.num_rows;
    
    pthread_create(&thread_handles[i], NULL, gauss_parallel, (void *)args);
    
  }
}

void* gauss_parallel(void *arg)
{
  unsigned int i, j, k;
  struct Args *current = (struct Args*)arg;
  unsigned int start, num_elements;
  
  float *U = (*current).src;
  start = (*current).start;
  num_elements = (*current).num_elements;

  for (k = start; k < num_elements; k++)
  {
    // TODO: pthread(gauss_reduction)
    // Set the principal diagonal entry in U to be 1
    U[num_elements * k + k] = 1;	 
    // TODO: pthread(gauss_elimination)
    
  }
}

// TODO: turn this parallel
void* gauss_reduction(void *arg)
{
  for (j = (k + 1); j < num_elements; j++)
  {			
    // Reduce the current row
    if (U[num_elements * k + k] == 0)
      {
        printf("Numerical instability detected. The principal diagonal element is zero. \n");
        /*break;*/
      }

      // Division step
      U[num_elements * k + j] = (float) (U[num_elements * k + j] / U[num_elements * k + k]);
  }
}
// TODO: turn this parallel
void* gauss_elimination(void *arg)
{

  for (i = (k + 1); i < num_elements; i++)
  {
    for (j = (k + 1); j < num_elements; j++)
      // Elimination step
      U[num_elements * i + j] = U[num_elements * i + j] - (U[num_elements * i + k] * U[num_elements * k + j]);	
    U[num_elements * i + k] = 0;
  }
}

/* Function checks if the results generated by the single threaded and multi threaded versions match. */
int
check_results (float *A, float *B, unsigned int size, float tolerance)
{
  for (int i = 0; i < size; i++)
    if (fabsf (A[i] - B[i]) > tolerance)
      return 0;
  return 1;
}


/* Allocate a matrix of dimensions height*width
 * If init == 0, initialize to all zeroes.  
 * If init == 1, perform random initialization. 
*/
Matrix
allocate_matrix (int num_rows, int num_columns, int init)
{
  Matrix M;
  M.num_columns = M.pitch = num_columns;
  M.num_rows = num_rows;
  int size = M.num_rows * M.num_columns;

  M.elements = (float *) malloc (size * sizeof (float));
  for (unsigned int i = 0; i < size; i++)
    {
      if (init == 0)
	M.elements[i] = 0;
      else
	M.elements[i] = get_random_number (MIN_NUMBER, MAX_NUMBER);
    }
  return M;
}


/* Returns a random floating-point number between the specified min and max values. */ 
float
get_random_number (int min, int max)
{
  return (float)
    floor ((double)
	   (min + (max - min + 1) * ((float) rand () / (float) RAND_MAX)));
}

/* Performs a simple check on the upper triangular matrix. Checks to see if the principal diagonal elements are 1. */
int
perform_simple_check (const Matrix M)
{
  for (unsigned int i = 0; i < M.num_rows; i++)
    if ((fabs (M.elements[M.num_rows * i + i] - 1.0)) > 0.001)
      return 0;
  return 1;
}
