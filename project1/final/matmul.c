// ConApp.cpp : Defines the entry point for the console application.
//
/*
#include "stdafx.h"
int _tmain(int argc, _TCHAR* argv[])
{
	return 0;
}
*/
//#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
//#include "c_timer.h"
#include <papi.h>
#include "papi_timer.h"
#define NUM 5

double *A, *B, *C, *D;
int size[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
//int size[] = {10, 20, 30, 40, 50};//, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

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

void dgemmJKI(int N)
{
	for(int j=0; j<N; j++){
		for(int k=0; k<N; k++){
			for(int i=0; i<N; i++){
				C[j*N+i] += A[j*N+k] * B[k*N+i];
			}
			//C[j*N+i] = cji;
		}
	}
}

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

	// printf("D=[");
	// for(int i=0;i<N;i++)
	// 	printf("%.2lf,", D[i]);
	// printf("]\n");	
	for(int i=0;i<N;i++)
		C[i] = D[i];	
	
}

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

/*
void calcTime( void(*f)(int), char* com )
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
}
*/

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
			//dgemmIJK(n);
			(*f)(n);
			etime = get_cur_time();
			sum += (etime-btime);
			//printf("Runtime of %s for n=%d is %lf \n", com, n, (etime-btime));
		}
		avg = sum/NUM;
		gflop = (2.0*n*n*n)/(avg*1000000000.0);
		//printf("Average runtime of %s for n=%d is %lf \n", com, n, avg);
		printf("%s, %d, %lf, %lf \n", com, n, avg, gflop);
	}
}

void calcPapiEvents(void(*f)(int), char* com)
{
	int len = sizeof(size)/sizeof(size[0]);
	int Events[] = {PAPI_L1_TCM, PAPI_L2_TCM, PAPI_L1_TCH, PAPI_L1_TCA};
	int elen = sizeof(Events)/sizeof(Events[0]);
	long long values[4];	//should be dynamcally allocated
	int retval;

	printf("Order, N, L1_TCM, L2_TCM, L1_TCH, L1_TCA, L1_Miss_Rate\n");
	for(int i=0; i<len; i++){
		int n = size[i];
		init(n*n);

		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		
		(*f)(n);
		//dgemmIJK(n);

		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);

		printf("%s, %d, %d, %d, %d, %d, %lf\n", com, n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]); 	
	    //printf("L1_TCM=%d, L2_TCM=%d of %s for n=%d\n", values[0], values[1], com, n);
	    //printf("L1_TCM=%d, L2_TCM=%d, L1_TCH=%d, PAPI_L1_TCA=%d of %s for n=%d\n", values[0], values[1], values[2], values[3], com, n);
	    //printf("%s, %d, %d, %d, %d, %d\n", com, n, values[0], values[1], values[2], values[3]);

/*		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		dgemmIKJ(n);
		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);
	    printf("%d, IKJ, %d, %d, %d, %d, %.2f\n", n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]);

		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		dgemmJIK(n);
		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);
	    //printf("JIK, %d, %d, %d, %d, %d\n", n, values[0], values[1], values[2], values[3]);
		printf("%d, JIK, %d, %d, %d, %d, %.2f\n", n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]);

		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		dgemmJKI(n);
		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);
	    //printf("JKI, %d, %d, %d, %d, %d\n", n, values[0], values[1], values[2], values[3]);
		printf("%d, JKI, %d, %d, %d, %d, %.2f\n", n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]);

		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		dgemmKIJ(n);
		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);
	    //printf("KIJ, %d, %d, %d, %d, %d\n", n, values[0], values[1], values[2], values[3]);
		printf("%d, KIJ, %d, %d, %d, %d, %.2f\n", n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]);

		if ((retval=PAPI_start_counters(Events, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_start_counters", retval);
		dgemmKJI(n);
		if ((retval=PAPI_stop_counters(values, elen)) != PAPI_OK)
			test_fail(__FILE__, __LINE__, "PAPI_stop_counters", retval);
	    //printf("KJI, %d, %d, %d, %d, %d\n", n, values[0], values[1], values[2], values[3]);
		printf("%d, KJI, %d, %d, %d, %d, %.2f\n", n, values[0], values[1], values[2], values[3], values[0]*100.0/values[3]);
*/
		free(A);
		free(B);
		free(C);
		free(D);
	}
	//printf("%s\tPASSED\n", __FILE__);
	//PAPI_shutdown();
}

int main()
{
	int retval;
	srand(time(NULL));
	if ((retval=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT){
		test_fail(__FILE__, __LINE__, "PAPI_library_init", retval);
	}

/*	clock_t start_t, end_t;
	double diff_t;
	const int N = 3;
	start_t = clock();	
	double btime = get_cur_time();

	init(N*N);	
	print(N*N,1);
	dgemmIJK(N);
	print(N*N,0);
	dgemmIKJ(N);
	print(N*N,0);
	dgemmJIK(N);
	print(N*N,0);
	dgemmJKI(N);
	print(N*N,0);
	printf("Hello World\n");

	end_t = clock();
	double etime = get_cur_time();
	diff_t = (double)(end_t-start_t)/CLOCKS_PER_SEC;
	printf("Execution time of IJK = %f\n", diff_t);
	printf("Execution time = %f\n", etime-btime);

	calcTime(dgemmIJK, "IJK");
	calcTime(dgemmIKJ, "IKJ");
	calcTime(dgemmJIK, "JIK");
	calcTime(dgemmJKI, "JKI");
*/
	calcPapiTime(dgemmIJK, "IJK");
	calcPapiTime(dgemmIKJ, "IKJ");
	calcPapiTime(dgemmJIK, "JIK");
	calcPapiTime(dgemmJKI, "JKI");
	calcPapiTime(dgemmKIJ, "KIJ");
	calcPapiTime(dgemmKJI, "KJI");

	// calcPapiEvents(dgemmIJK, "IJK");
	// calcPapiEvents(dgemmIKJ, "IKJ");
	// calcPapiEvents(dgemmJIK, "JIK");
	// calcPapiEvents(dgemmJKI, "JKI");
	// calcPapiEvents(dgemmKIJ, "KIJ");
	// calcPapiEvents(dgemmKJI, "KJI");	

	// free(A);
	// free(B);
	// free(C);
	// free(D);

	PAPI_shutdown();
	return 0;
}
