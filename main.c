#include <mpi.h>
#include "omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOLERANCE 0.00001

long usecs();
void initialize(double **A, int rows, int cols);
int calc_serial(double **A, int rows, int cols, int iters, double tolerance);
int calc_serial_v1(double **A, int rows, int cols, int iters, double tolerance);
int calc_omp(double **A, int rows, int cols, int iters, double tolerance, int num_threads);
double verify(double **A, double **B, int rows, int cols);

int main(int argc, char** argv) {
    int rank = 0;
    int comm_size = 0;
    int tag = 0;
    int retVal = 0;
    int choice = 0;
    int i, j;
    if (argc > 1) {
        choice = atoi(argv[1]);
    }
    int iters = 1;
    if (argc > 2) {
        iters = atoi(argv[2]);
    }
    int rows = 100, cols = 100;
    if (argc > 3) {
        rows = atoi(argv[3]);
        cols = rows;
    }

    int num_threads = 2;
    if (argc > 4) {
        num_threads = atoi(argv[4]);
    }

    double **A = malloc((rows+1)*sizeof(double*));
    double **orig = malloc((rows+1)*sizeof(double*));
    for (i = 0; i < rows+1; ++i) {
        A[i] = malloc((cols+1) * sizeof(double));
        orig[i] = malloc((cols+1) * sizeof(double));
    }

    initialize(A, rows, cols);
    initialize(orig, rows, cols);

    long startTime = 0;
    long endTime = 0;
    long diffTime = 0;

    calc_serial_v1(orig, rows, cols, iters, TOLERANCE);
    
    startTime = usecs();
    if (choice == 0) {
        calc_serial(A, rows, cols, iters, TOLERANCE);
    } else if (choice == 1) {
        calc_omp(A, rows, cols, iters, TOLERANCE, num_threads);
    }
    endTime = usecs();
    diffTime = endTime-startTime;
    
    double err = verify(orig, A, rows, cols);
    
    printf("Time = %ld us, with error %f", diffTime, err);
    
    return 0;
}

long usecs() {
    double time = MPI_Wtime();
    return (long)(time*1000000);
}

void initialize(double **A, int rows, int cols) {
    int i,j;

    for (j = 0; j < cols+1;j++){
        A[0][j]=1.0;
    }

    for (i = 1; i < rows+1; i++){
        A[i][0]=1.0;
        for (j = 1; j < cols+1;j++) 
            A[i][j]=0.0;
    }
}

double verify(double **A, double **B, int rows, int cols) {
    double error = 0.0;

    int i,j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            error += fabs(A[i][j] - B[i][j]);
        }
    } 
    return error;
}

int calc_serial(double **A, int rows, int cols, int iters, double tolerance) {
    int convergence = 0;
    double diff, tmp;
    int i,j;
    int for_iters;


    for (for_iters = 0; for_iters < iters; for_iters++) 
    { 
        diff = 0.0;

        for (i = 1;i < rows;i++)
        {
            for (j = 1; j < cols; j++)
            {
                tmp = A[i][j];
                A[i][j] = 0.2 * (A[i][j] + A[i][j-1] + A[i-1][j] + A[i][j+1] + A[i+1][j]);
                diff += fabs(A[i][j] - tmp);
            }
        }

        if (diff/((double)rows*(double)cols) < tolerance)
            convergence=1;

    }

    //printf("\nConv = %d", convergence);
    return convergence;
}

int calc_serial_v1(double **A, int rows, int cols, int iters, double tolerance) {
    int convergence = 0;
    double diff, tmp;
    int i,j;
    int for_iters;

    double** B = (double**)malloc((rows+1)*sizeof(double*));
    for (i = 0; i < rows+1; ++i) {
        B[i] = (double*)malloc((cols+1) * sizeof(double));
    }

    for (for_iters = 0; for_iters < iters; for_iters++) 
    { 
        for (i = 0; i < rows+1; ++i) {
            for (j = 0; j < cols+1; ++j) {
                B[i][j] = A[i][j];
            }
        }

        diff = 0.0;

        for (i = 1;i < rows;i++)
        {
            for (j = 1; j < cols; j++)
            {
                A[i][j] = 0.2 * (B[i][j] + B[i][j-1] + B[i-1][j] + B[i][j+1] + B[i+1][j]);
                diff += fabs(A[i][j] - B[i][j]);
            }
        }

        if (diff/((double)rows*(double)cols) < tolerance)
            convergence=1;

    }

    //printf("\nConv = %d", convergence);
    return convergence;
}

int calc_omp(double **A, int rows, int cols, int iters, double tolerance, int num_threads) {
    int convergence = 0;
    double diff;
    int i,j;
    int for_iters;

    double** B = (double**)malloc((rows+1)*sizeof(double*));
    for (i = 0; i < rows+1; ++i) {
        B[i] = (double*)malloc((cols+1) * sizeof(double));
    }

    for (for_iters = 0; for_iters < iters; for_iters++)
    {
#pragma omp parallel for private(j) num_threads(num_threads)
        for (i = 0; i < rows+1; ++i) {
            for (j = 0; j < cols+1; ++j) {
                B[i][j] = A[i][j];
            }
        }

        diff = 0.0;

#pragma omp parallel for reduction(+:diff) private(j) num_threads(num_threads)
        for (i = 1;i < rows;i++)
        {
            for (j = 1; j < cols; j++)
            {
                A[i][j] = 0.2 * (B[i][j] + B[i][j-1] + B[i-1][j] + B[i][j+1] + B[i+1][j]);
                diff += fabs(A[i][j] - B[i][j]);
            }
        }

        if (diff/((double)rows*(double)cols) < tolerance)
            convergence=1;
    } 

    //printf("\nConvP = %d", convergence);
    return convergence;
}
