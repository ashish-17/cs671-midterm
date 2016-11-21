#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define Tolerance 0.00001
#define TRUE 1
#define FALSE 0

#define N 5000

double ** A;

void initialize (double **A, int n)
{
	int i,j;

	for (j=0;j<n+1;j++){
		A[0][j]=1.0;
	}
	for (i=1;i<n+1;i++){
		A[i][0]=1.0;
		for (j=1;j<n+1;j++) A[i][j]=0.0;
	}

}

void solve(double **A, int n)
{
	int convergence=FALSE;
	double diff, tmp;
	int i,j, iters=0;
	int for_iters;


	for (for_iters=1;for_iters<21;for_iters++) 
	{ 
		diff = 0.0;

		for (i=1;i<n;i++)
		{
			for (j=1;j<n;j++)
			{
				tmp = A[i][j];
				A[i][j] = 0.2*(A[i][j] + A[i][j-1] + A[i-1][j] + A[i][j+1] + A[i+1][j]);
				diff += fabs(A[i][j] - tmp);
			}
		}
		iters++;

		if (diff/((double)N*(double)N) < Tolerance)
			convergence=TRUE;

	} /*for*/
}


long usecs (void)
{
	struct timeval t;

	gettimeofday(&t,NULL);
	return t.tv_sec*1000000+t.tv_usec;
}

__global__ void calc_gpu(int* d_arr) {
	d_arr[threadIdx.x]++;
}

int main(int argc, char * argv[])
{
	int i;
	long t_start,t_end;
	double time;

	dim3 grid(10);
	dim3 block(1);
	
	int h_arr[] = {1,2,3,4,5,6,7,8,9,10};
	int *d_arr;
	cudaMalloc(&d_arr, 10*sizeof(int));
	cudaMemcpy(d_arr, h_arr, 10*sizeof(int), cudaMemcpyHostToDevice);
 
	calc_gpu<<<1,10>>>(d_arr);
	
	cudaMemcpy(h_arr, d_arr, 10*sizeof(int), cudaMemcpyDeviceToHost);
	for (i=0; i < 10; ++i) {
		printf(" %d", h_arr[i]);
	}

	A = (double**)malloc((N+2) * sizeof(double *));
	for (i=0; i<N+2; i++) {
		A[i] = (double*)malloc((N+2) * sizeof(double)); 
	}

	initialize(A, N);

	t_start=usecs();
	solve(A, N);
	t_end=usecs();

	time = ((double)(t_end-t_start))/1000000;
	printf("Computation time = %f\n", time);

}