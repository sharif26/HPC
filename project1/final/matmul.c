#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include <papi.h>
#include "papi_timer.h"
#define NUM 5

double *A, *B, *C;
int size[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

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

	for(int i=0;i<N;i++){
		A[i] = (double)rand()*100/RAND_MAX;
		B[i] = (double)rand()*100/RAND_MAX;
		C[i] = (double)rand()*100/RAND_MAX;
	}
}

/*Matrix multiplication using IJK combination*/
void dgemmIJK(int N)
{
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
		{
			double cij = C[i*N+j];
			for(int k=0; k<N; k++)
				cij += A[i*N+k] * B[k*N+j];
			C[i*N+j] = cij;
		}
}

/*Matrix multiplication using IKJ combination*/
void dgemmIKJ(int N)
{
	for(int i=0; i<N; i++){
		for(int k=0; k<N; k++){
			for(int j=0; j<N; j++){
				C[i*N+j] += A[i*N+k] * B[k*N+j];
			}
		}
	}
}

/*Matrix multiplication using JIK combination*/
void dgemmJIK(int N)
{
	for(int j=0; j<N; j++){
		for(int i=0; i<N; i++){
			double cji = C[j*N+i];
			for(int k=0; k<N; k++){
				cji += A[j*N+k] * B[k*N+i];
			}
			C[j*N+i] = cji;
		}
	}
}

/*Matrix multiplication using JKI combination*/
void dgemmJKI(int N)
{
	for(int j=0; j<N; j++){
		for(int k=0; k<N; k++){
			for(int i=0; i<N; i++){
				C[j*N+i] += A[j*N+k] * B[k*N+i];
			}
		}
	}
}

/*Matrix multiplication using KIJ combination*/
void dgemmKIJ(int N)
{
	for(int k=0; k<N; k++){
		for(int i=0; i<N; i++){
			double cki = C[k*N+i];
			for(int j=0; j<N; j++){
				cki += A[k*N+j] * B[j*N+i];
			}
			C[k*N+i] = cki;
		}
	}
}

/*Matrix multiplication using KJI combination*/
void dgemmKJI(int N)
{
	for(int k=0; k<N; k++){
		for(int j=0; j<N; j++){
			for(int i=0; i<N; i++){
				C[k*N+i] += A[k*N+j] * B[j*N+i];
			}
		}

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
			if( (i+1)%n == 0 ) printf("],");
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
		if( (i+1)%n == 0 ) printf("],");
		else printf(",");
	}
	printf("]\n");

}

/*function used to handle PAPI library's error situation (got from PAPI tutorial)*/
void test_fail(char *file, int line, char *call, int retval){
    printf("%s\tFAILED\nLine # %d\n", file, line);
    if ( retval == PAPI_ESYS ) {
        char buf[128];
        memset( buf, '\0', sizeof(buf) );
        sprintf(buf, "System error in %s:", call );
        perror(buf);
    }
    else if ( retval > 0 ) {
        printf("Error calculating: %s\n", call );
    }
    else {
        printf("Error in %s: %s\n", call, PAPI_strerror(retval) );
    }
    printf("\n");
    exit(1);
}

/*Calculate execution time using C timer function (was not used to draw plots)*/
/*void calcTime( void(*f)(int), char* com )
{
	int len = sizeof(size)/sizeof(size[0]);
	double btime, etime;
	for(int i=0; i<len; i++)
	{
		int n = size[i];
		double sum = 0.0, avg = 0.0;
		init(n*n);
		for(int j=0; j<NUM; j++){
			btime = get_cur_time();
			//dgemmIJK(n);
			(*f)(n);
			etime = get_cur_time();
			sum += (etime-btime);
			//printf("Runtime of %s for n=%d is %lf \n", com, n, (etime-btime));
		}
		avg = sum/NUM;
		printf("Average runtime of %s for n=%d is %lf \n", com, n, avg);
	}
}*/

/*Calculate Gflop/s and execution time using PAPI library's time function*/
void calcPapiTime(void(*f)(int), char* com)
{
	int len = sizeof(size)/sizeof(size[0]);
	double btime, etime, gflop;
	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
		printf("PAPI can not be initialized!");
		exit(1);	
	}
	for(int i=0; i<len; i++)
	{
		int n = size[i];
		double sum = 0.0, avg = 0.0;
		init(n*n);
		for(int j=0; j<NUM; j++){
			btime = get_cur_time();
			(*f)(n);
			etime = get_cur_time();
			sum += (etime-btime);
		}
		avg = sum/NUM;
		gflop = (2.0*n*n*n)/(avg*1000000000.0);
		printf("%s, %d, %lf, %lf \n", com, n, avg, gflop);
	free(A);
	free(B);
	free(C);
	}
}

/*Calculate cache miss rate using PAPI library's high level events*/
void calcPapiEvents(void(*f)(int), char* com)
{
	int len = sizeof(size)/sizeof(size[0]);
	int Events[] = {PAPI_L1_TCM, PAPI_L2_TCM, PAPI_L1_TCH, PAPI_L1_TCA};
	int elen = sizeof(Events)/sizeof(Events[0]);
	long long values[4];
	int retval;

	for(int i=0; i<len; i++){
		int n = size[i];
		init(n*n);

		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		
		(*f)(n);

		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);

		printf("%s, %d, %d, %d, %d, %d, %lf\n", com, n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]); 	
		free(A);
		free(B);
		free(C);
	}
}

int main()
{
	srand(time(NULL));

	printf("Order, N, Execution Time, Gflop/s\n");
	calcPapiTime(dgemmIJK, "IJK");
	calcPapiTime(dgemmIKJ, "IKJ");
	calcPapiTime(dgemmJIK, "JIK");
	calcPapiTime(dgemmJKI, "JKI");
	calcPapiTime(dgemmKIJ, "KIJ");
	calcPapiTime(dgemmKJI, "KJI");
	
	printf("Order, N, L1_TCM, L2_TCM, L1_TCH, L1_TCA, L1_Miss_Rate\n");
	calcPapiEvents(dgemmIJK, "IJK");
	calcPapiEvents(dgemmIKJ, "IKJ");
	calcPapiEvents(dgemmJIK, "JIK");
	calcPapiEvents(dgemmJKI, "JKI");
	calcPapiEvents(dgemmKIJ, "KIJ");
	calcPapiEvents(dgemmKJI, "KJI");	

	return 0;
}
