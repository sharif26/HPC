#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "time.h"
#include "math.h"
#include "c_timer.h"
#include <mpi.h>
#define TOL 1e-4

#define DGEMM dgemm_
extern void DGEMM (char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
//const char* dgemm_desc = "Reference dgemm.";
/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are N-by-N matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 * This function wraps a call to the BLAS-3 routine DGEMM, via the standard FORTRAN interface - hence the reference semantics. */
void square_dgemm (int N, double* A, double* B, double* C)
{
  char TRANSA = 'N';
  char TRANSB = 'N';
  int M = N;
  int K = N;
  double ALPHA = 1.;
  double BETA = 1.;
  int LDA = N;
  int LDB = N;
  int LDC = N;
  DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}

int myrank;
int nprocs;

void print_matrix(const int rows, const int cols, const double *matr)
{
	int N = rows * cols;
	int n = sqrt(N);
		printf("[");
		for(int i=0;i<N;i++){
			if( i%n == 0 ) printf("[");
			printf("%.2lf", matr[i]);
			if( i==N-1 ) printf("]]\n");
			else if( (i+1)%n == 0 ) printf("],\n");
			else printf(",");
		}
		//printf("]\n");

}

double validate(const int n, const double *Csumma, double *Cnaive) {
    double eps = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i*n + j;
            Cnaive[idx] = fabs(Cnaive[idx] - Csumma[idx]);
            if ( eps < Cnaive[idx] ){
                eps = Cnaive[idx];
            }
        }
    }
    return eps;
}

void SUMMA(MPI_Comm comm_cart, const int n, const int nb, double *A_loc, double *B_loc, double *C_loc) {
    int coords[2];
    MPI_Cart_coords(comm_cart, myrank, 2, coords);
    int my_col = coords[0];
    int my_row = coords[1];

    MPI_Comm row_comm;
    MPI_Comm col_comm;
    int remain_dims[2];

    // create row comms for A
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(comm_cart, remain_dims, &row_comm);

    // create col comms for B
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(comm_cart, remain_dims, &col_comm);

    double *A_loc_save = (double *) calloc(nb*nb, sizeof(double));
    double *B_loc_save = (double *) calloc(nb*nb, sizeof(double));
    double *C_loc_tmp = (double *) calloc(nb*nb, sizeof(double));

    memcpy(A_loc_save, A_loc, nb*nb*sizeof(double));
    memcpy(B_loc_save, B_loc, nb*nb*sizeof(double));

    //memset(C_loc, 0, nb*nb*sizeof(double));
    int nblks = n / nb;		//this makes nblks to nproc ??

	for (int b=0; b<nblks; b++) {
       int root_col = b;
       int root_row = b;

       if (my_col == root_col) {
           memcpy(A_loc, A_loc_save, nb*nb*sizeof(double));
       }
       // broadcast A_loc from root_col within row_comm
		MPI_Bcast(A_loc, nb*nb, MPI_DOUBLE, root_col, row_comm);
		//}
       if (my_row == root_row) {
           memcpy(B_loc, B_loc_save, nb*nb*sizeof(double));
       }
       // broadcast B_loc from root_row within col_comm
		MPI_Bcast(B_loc, nb*nb, MPI_DOUBLE, root_row, col_comm);
		//}
	   //here we should use Cray LibSci’s dgemm function
	   //matmul_naive(nb, A_loc, B_loc, C_loc_tmp);
	   //matmul_naive(nb, A_loc, B_loc, C_loc_tmp);
	   square_dgemm (nb, A_loc, B_loc, C_loc);

	   //plus_matrix(nb, C_loc, C_loc_tmp, C_loc);
	}

    free(A_loc_save);
    free(B_loc_save);
    free(C_loc_tmp);
}

void scatterMatrix(int n)
{
	//MPI_Init(&argc, &argv);
    int p, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    //char i;
	int ROWS = n, COLS = n;

    //char a[ROWS*COLS];
    double A[n*n];
    double B[n*n];
    double C[n*n];
    double D[n*n];
    const int NPROWS=sqrt(p);//2;  /* number of rows in _decomposition_ */
    const int NPCOLS=sqrt(p);//3;  /* number of cols in _decomposition_ */
    const int BLOCKROWS = ROWS/NPROWS;  /* number of rows in _block_ */
    const int BLOCKCOLS = COLS/NPCOLS; /* number of cols in _block_ */

    if (p != NPROWS*NPCOLS) {
        fprintf(stderr,"Error: number of PEs %d != %d x %d\n", p, NPROWS, NPCOLS);
        MPI_Finalize();
        exit(-1);
    }
    double A_loc[BLOCKROWS*BLOCKCOLS];
    double B_loc[BLOCKROWS*BLOCKCOLS];
    double C_loc[BLOCKROWS*BLOCKCOLS];

	if (myrank == 0) {
        for (int i=0; i<n*n; i++) {
            A[i] = (double)i;
            B[i] = (double)(i+i);
            C[i] = (double)i;
            D[i] = (double)i;
        }
		//this is to validate the result
		square_dgemm (n, A, B, D);
    }
	
    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_DOUBLE, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(double), &blocktype);
    MPI_Type_commit(&blocktype);
    int disps[NPROWS*NPCOLS];
    int counts[NPROWS*NPCOLS];
    for (int ii=0; ii<NPROWS; ii++) {
        for (int jj=0; jj<NPCOLS; jj++) {
            disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;
            counts [ii*NPCOLS+jj] = 1;
        }
    }

    MPI_Scatterv(A, counts, disps, blocktype, A_loc, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(B, counts, disps, blocktype, B_loc, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(C, counts, disps, blocktype, C_loc, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int ndims = 2;
    const int dims[2] = {NPROWS, NPCOLS};
    const int periods[2] = {0, 0};
    int reorder = 0;
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);

    double tstart, tend;
    MPI_Barrier(MPI_COMM_WORLD);
    tstart = MPI_Wtime();
    SUMMA(comm_cart, n, BLOCKROWS, A_loc, B_loc, C_loc);
    tend = MPI_Wtime();

    double etime = tend - tstart;
    double max_etime = 0.0;

    MPI_Reduce(&etime, &max_etime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        printf("n=%d, p=%d, exec=%lf\n", n, p, max_etime);
        //printf("Gflops %f\n", (2.0*n*n*n)/(max_etime*1000000000.0));
		//print_matrix(ROWS, COLS, A);
		//print_matrix(ROWS, COLS, B);
		//print_matrix(ROWS, COLS, C);
    }

    MPI_Gatherv(C_loc, BLOCKROWS*BLOCKCOLS, MPI_DOUBLE, C, counts, disps, blocktype, 0, MPI_COMM_WORLD);

	//print block matrix from all processes
/*    for (int proc=0; proc<p; proc++) {
        if (proc == myrank) {
            printf("Rank = %d\n", myrank);
            if (myrank == 0) {
                printf("Global matrix: \n");
                for (int ii=0; ii<ROWS; ii++) {
                    for (int jj=0; jj<COLS; jj++) {
                        printf("%5.2f ",C[ii*COLS+jj]);
                    }
                    printf("\n");
                }
            }
            printf("Local Matrix:\n");
            for (int ii=0; ii<BLOCKROWS; ii++) {
                for (int jj=0; jj<BLOCKCOLS; jj++) {
                    printf("%5.2f ",C_loc[ii*BLOCKCOLS+jj]);
                }
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/ 
	
	//validate the result
/*	if(myrank == 0){
		double eps = validate(n, C, D);
		if (eps > TOL) {
			fprintf(stderr, "ERROR: eps = %f\n", eps);
			//MPI_Abort(MPI_COMM_WORLD, 1);
		} else {
			printf("SUMMA: OK: eps = %f\n", eps);
		}
	}*/
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	srand(time(NULL));
	 if(argc!=2){
	 	printf("Usage: aprun -n <number of procs> ./summa <n>\n");
	 	exit(1);
	 }
	int n = atoi(argv[1]);
	scatterMatrix(n);
	MPI_Finalize();
	return 0;
}

