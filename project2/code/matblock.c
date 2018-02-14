#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "c_timer.h"
#include "string.h"
#include "emmintrin.h"
#define LOOP 4

double *A, *B, *C, *D;
int size[] = {64, 128, 256, 512, 1024, 2048};

static double BLOCK_A_64[64*64] __attribute__((aligned(16)));
static double BLOCK_A_128[128*128] __attribute__((aligned(16)));
static double BLOCK_A_256[256*256] __attribute__((aligned(16)));

static double BLOCK_B_64[64*64] __attribute__((aligned(16)));
static double BLOCK_B_128[128*128] __attribute__((aligned(16)));
static double BLOCK_B_256[256*256] __attribute__((aligned(16)));

static double BLOCK_C_64[64*64] __attribute__((aligned(16)));
static double BLOCK_C_128[128*128] __attribute__((aligned(16)));
static double BLOCK_C_256[256*256] __attribute__((aligned(16)));


/*dynamically allocating Matrices as 1D linear array and then initializing with random floating values*/
void init(int N)
{
	A = (double*) calloc(N, sizeof(double));
	if(A==NULL){
		printf("Error A");
		exit(1);
	}
	B = (double*) calloc(N, sizeof(double));
	if(B==NULL){
		printf("Error B");
		exit(1);
	}
	C = (double*) calloc(N, sizeof(double));
	if(C==NULL){
		printf("Error C");
		exit(1);
	}
	D = (double*) calloc(N, sizeof(double));
	if(D==NULL){
		printf("Error D");
		exit(1);
	}
	for(int i=0;i<N;i++){
		A[i] = (double)rand()*100/RAND_MAX;
		B[i] = (double)rand()*100/RAND_MAX;
		C[i] = (double)rand()*100/RAND_MAX;
        D[i] = C[i];
	}
}

/*printing Matrices to check or debug*/
void print(int N, int AB)
{
	int n = sqrt(N);
	if(AB==1){
		printf("A=[");
		for(int i=0;i<N;i++){
			if( i%n == 0 ) printf("[");
			printf("%.2lf", A[i]);
			if( (i+1)%n == 0 ) printf("],\n");
			else printf(",");
		}
		printf("]\n");
	
		printf("B=[");
		for(int i=0;i<N;i++){
			if( i%n == 0 ) printf("[");			
			printf("%.2lf", B[i]);
			if( (i+1)%n == 0 ) printf("],");
			else printf(",");
		}
		printf("]\n");
	}
	printf("C=[");
	for(int i=0;i<N;i++){
		if( i%n == 0 ) printf("[");
		printf("%.2lf", C[i]);
		if( (i+1)%n == 0 ) printf("],\n");
		else printf(",");
	}
	printf("]\n");

	// printf("D=[");
	// for(int i=0;i<N;i++)
	// 	printf("%.2lf,", D[i]);
	// printf("]\n");	
	for(int i=0;i<N;i++)
		C[i] = D[i];    
}

/*Matrix multiplication using IJK combination to check error*/
void dgemmIJK(int N)
{
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			double cij = D[i*N+j];
			for(int k=0; k<N; k++)
				cij += A[i*N+k] * B[k*N+j];
			D[i*N+j] = cij;
		}
		
	//checking error	
	double error = 0.0;
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			error += fabs(D[i*N+j] - C[i*N+j]);
	printf("Error=%lf, ", error);		
		
}

/*blocked matrix multiplication*/
void dgemmTiled(int n, int b)
{
	int N = n/b;
	int bb = sizeof(double)*b*b;
	int db = sizeof(double)*b;
	static double *Caux, *Aaux, *Baux;
	for(int i=0; i<N; i++){
		int ibn = i*b*n;
		for(int j=0; j<N; j++){
			int jb = j*b;
			if(b==64)
				Caux = BLOCK_C_64;
			else if(b==128)
				Caux = BLOCK_C_128;
			else if(b==256)
				Caux = BLOCK_C_256;
			else
				posix_memalign((void **)&Caux, 16, bb);//bb=sizeof(double)*b*b
			
			if(b==64)
				Aaux = BLOCK_A_64;
			else if(b==128)
				Aaux = BLOCK_A_128;
			else if(b==256)
				Aaux = BLOCK_A_256;
			else
				posix_memalign((void **)&Aaux, 16, bb);

			if(b==64)
				Baux = BLOCK_B_64;
			else if(b==128)
				Baux = BLOCK_B_128;
			else if(b==256)
				Baux = BLOCK_B_256;
			else
				posix_memalign((void **)&Baux, 16, bb);

			//load C to fast memory
			for(int s=0; s<b; s+=LOOP){
				memcpy(Caux+s*b, C+ibn+jb+s*n, db);
				memcpy(Caux+(s+1)*b, C+ibn+jb+(s+1)*n, db);
				memcpy(Caux+(s+2)*b, C+ibn+jb+(s+2)*n, db);
				memcpy(Caux+(s+3)*b, C+ibn+jb+(s+3)*n, db);
			}

			for(int k=0; k<N; k++){
				int kb = k*b;
				int kbn = k*b*n;
				//load A to fast memory
				for(int s=0; s<b; s+=LOOP){
					memcpy(Aaux+s*b, A+ibn+kb+s*n, db);
					memcpy(Aaux+(s+1)*b, A+ibn+kb+(s+1)*n, db);
					memcpy(Aaux+(s+2)*b, A+ibn+kb+(s+2)*n, db);
					memcpy(Aaux+(s+3)*b, A+ibn+kb+(s+3)*n, db);
				}
				//load B to fast memory					
				for(int s=0; s<b; s+=LOOP){
					memcpy(Baux+s*b, B+kbn+jb+s*n, db);
					memcpy(Baux+(s+1)*b, B+kbn+jb+(s+1)*n, db);
					memcpy(Baux+(s+2)*b, B+kbn+jb+(s+2)*n, db);
					memcpy(Baux+(s+3)*b, B+kbn+jb+(s+3)*n, db);
				}
				
				//IKJ
				for(int ii=0; ii<b; ii++)
					for(int kk=0; kk<b; kk++){
					register double r = Aaux[ii*b+kk];
					for(int jj=0; jj<b; jj+=LOOP){
						Caux[ii*b+jj] += r*Baux[kk*b+jj];
						Caux[ii*b+jj+1] += r*Baux[kk*b+jj+1];
						Caux[ii*b+jj+2] += r*Baux[kk*b+jj+2];
						Caux[ii*b+jj+3] += r*Baux[kk*b+jj+3];
					}
				}			
			}
			
			for(int s=0; s<b; s+=LOOP){
				memcpy(C+ibn+jb+s*n, Caux+s*b, db);
				memcpy(C+ibn+jb+(s+1)*n, Caux+(s+1)*b, db);
				memcpy(C+ibn+jb+(s+2)*n, Caux+(s+2)*b, db);
				memcpy(C+ibn+jb+(s+3)*n, Caux+(s+3)*b, db);
			}
		}
	}
	
}

void compute()
{
	srand(time(NULL));

	double btime, etime;
	int len = sizeof(size)/sizeof(size[0]);
	for(int i=0; i<len; i++)
	{
		int N = size[i];
		init(N*N);
		//print(N*N,1);
		for(int block=4; block<=N/2; block*=2){
			printf("%d, %d, ", N, block);
			double btime = get_cur_time();
			dgemmTiled(N, block);
			double etime = get_cur_time();
			double gflop = (2.0*N*N*N)/((etime-btime)*1000000000.0);
			//dgemmIJK(N);			
			printf("exec=%lf, gflop=%lf \n", (etime-btime), gflop);
		}

		free(A);
		free(B);
		free(C);
		free(D);
	}
}

int main()
{
	compute();
	return 0;
}
