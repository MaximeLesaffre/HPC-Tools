#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mkl_lapacke.h"
#include <omp.h>

double *generate_matrix(int size, double *matrix)
{
    int i;
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", name);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

double *Transposee(double *Q, int size, double *tQ) {	
    int i,j;
    #pragma omp parallel for private (i,j) schedule(static)
    for ( i = 0; i < size ; i++) {
        for ( j = 0; j < size ; j++) {				
            tQ[i*size+j] = Q[j*size+i];
        }
    }
    return tQ;
}


double *MatriceMatrice(double *matrix1, double *matrix2, int size, double *res) {
    int i,j,k;
    #pragma omp parallel for private (i,j,k) schedule(static)
	for(int i=0;i<size;i++) {
		for(int j=0;j<size;j++) {		// Parcourt du tableau 2D
			res[i*size+j]=0;		// Initilisation à 0 pour chaque élément de la matrice finale
				for(int k=0;k<size;k++) {		// Parcourt du tableau 2D
					res[i*size+j]+=matrix1[i*size+k]*matrix2[k*size+j];	// Les élements de la matrice vide sont calculés par la multiplication des éléments des 2 autres matrices et la somme
				}
		}
	}
    return res;
}


int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}

void QR(double *A, int size, double *Q, double *R) {
    int i,j,k;

    

    for (int i=0; i<size; i++) {
        double doo = 0;
        #pragma omp parallel for reduction(+:doo)
        for (int j=0; j<size; j+=2) {
            doo = doo + pow(A[j*size+i],2);
            doo = doo + pow(A[(j+1)*size+i],2);
        }
        R[i*size+i] = sqrt(doo);
        for (int j=0; j<size; j++) {
            Q[j*size+i] = A[j*size+i]/R[i*size+i];
        }
        #pragma omp parallel for private (j,k) schedule(dynamic) reduction(+:doo)
        for (int k=i+1; k<size; k++) {
            double doo = 0;
            
            for (int j=0; j<size; j++) {
                doo = doo + A[j*size+k]*Q[j*size+i];
            }
            R[i*size+k] = doo;
            for (int j=0; j<size; j++) {
                A[j*size+k] = A[j*size+k] - R[i*size+k]*Q[j*size+i];
            }
        }
    }
}

int Computation(double *A, double *B, int size, double *Q, double *R, double *X, double *tQ, double *res) {

    QR(A, size, Q, R);
    double *Y = MatriceMatrice(Transposee(Q, size, tQ), B, size, res);
    int i,j,k;
    for ( i=size-1; i>-1;i--) {
        #pragma omp parallel for private(j,k) schedule(dynamic)
        for ( j=size-1; j>-1;j--) {
            double s = 0;
            #pragma omp parallel for simd schedule(simd: static) reduction(+:s)
            for ( k=i+1;k<size;k++) {
                s += R[i*size+k]*X[k*size+j];
            }
        X[i*size+j] = (Y[i*size+j]-s)/R[i*(size+1)];
        }
    }
    return 1;
}


//int my_dgesv(int size, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {
//    
//
//
//    //Replace with your implementation
//    LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
//    
//}


    int main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *A = (double *)malloc(sizeof(double) * size * size);
        double *B = (double *)malloc(sizeof(double) * size * size);   
        double *Q = (double *)malloc(sizeof(double) * size * size);   
        double *R = (double *)malloc(sizeof(double) * size * size);   
        double *X = (double *)malloc(sizeof(double) * size * size);   
        double *tQ = (double *)malloc(sizeof(double) * size * size);   
        double *res = (double *)malloc(sizeof(double) * size * size);   




        A = generate_matrix(size, A);      
        B = generate_matrix(size, B);



        

        //print_matrix("A", a, size);
        //print_matrix("B", b, size);

        // Using MKL to solve the system
        //MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        //MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        //info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        //printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        //tStart = clock();    
        //MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);       // balek 
        //my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb); // balek
        Computation(A, B, size, Q, R, X, tQ, res);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        //if (check_result(bref,B,size)==1)
        //    printf("Result is ok!\n");
        //else    
        //    printf("Result is wrong!\n");
        
        //print_matrix("X", X, size);
        //print_matrix("Xref", bref, size);

        free(A);
        free(B);
        free(Q);
        free(R);
        free(X);
        free(tQ);
        free(res);
        return 0;
    }