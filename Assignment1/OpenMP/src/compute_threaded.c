#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include "gauss_eliminate.h"

#ifdef ARG_THREADS
#define NUM_THREADS ARG_THREADS
#else
#define NUM_THREADS 8           /* Number of threads. */
#endif

void gauss_eliminate_using_openmp(const Matrix A, Matrix U)
{
    unsigned int i, j, k;

    unsigned int rows = A.num_rows;
    unsigned int cols = A.num_columns;

    omp_set_num_threads(NUM_THREADS);

    /* Copy the contents of the A matrix into the U matrix. */
    for (i = 0; i < rows; ++i)
    {
        #pragma omp parallel for default(none) private(j) shared(U, i, rows, cols)
        for(j = 0; j < cols; ++j)
        {
            U.elements[rows * i + j] = A.elements[rows * i + j];
        }
    }

    unsigned int row;
    unsigned int col;
    unsigned int row_2;

    /* Perform Gaussian elimination in place on the U matrix. */
    for (row = 0; row < rows; ++row)
    {
        if (U.elements[rows * row + row] == 0)
        {
            printf("Numerical instability detected. The principal diagonal element is zero. \n");
            return;
        }

        /* Reduce the current row. */
        #pragma omp parallel for default(none) private(col) shared(U, row, rows, cols)
        for (col = (row + 1); col < cols; ++col)
        {
            /* Division step. */
            U.elements[rows * row + col] = (float)(U.elements[rows * row + col] / U.elements[rows * row + row]);
        }

        /* Set the principal diagonal entry in U to be 1. */
        U.elements[rows * row + row] = 1;

        #pragma omp parallel for default(none) private(row_2, col) shared(U, row, rows, cols)
        for (row_2 = (row + 1); row_2 < rows; ++row_2)
        {
            for (col = (row + 1); col < cols; ++col)
            {
                /* Elimination step. */
                U.elements[rows * row_2 + col] = U.elements[rows * row_2 + col]
                                               - ( U.elements[rows * row_2 + row]
                                                 * U.elements[rows * row + col]
                                                 );
            }

            U.elements[rows * row_2 + row] = 0;
        }
    }
}
