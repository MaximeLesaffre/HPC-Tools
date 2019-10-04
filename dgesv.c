#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mkl_lapacke.h"

double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
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

double *Transposee(double *Q, int size) {
    double *tQ = (double *)malloc(sizeof(double) * size * size);		
    for (int i = 0; i < size ; i++) {
        for (int j = 0; j < size ; j++) {				
            tQ[i*size+j] = Q[j*size+i];
        }
    }
    return tQ;
}

double *SousMatrice(double *flt_mat, int int_Ligneimat, int int_Colonnejmat, int size) {
    int int_Ligneires;
    int int_Colonnejres;
    double *newMat = (double *)malloc(sizeof(double) * (size-1) * (size-1));
    for (int i = 0 ; i < size ; i++) {
        for (int j = 0 ; j < size ; j++) {
            if ( i != int_Ligneimat && j != int_Colonnejmat) {
                int_Ligneires = (i > int_Ligneimat) ? i-1 : i;
                int_Colonnejres = (j > int_Colonnejmat) ? j-1 : j;
                newMat[int_Ligneires*(size-1)+int_Colonnejres] = flt_mat[i*size+j];
            }
        }
    }
    return newMat;
}



double CalculDeterminant(double *matrix, int size) {
    double coeff;
    double det=0;
    if (size == 1) {
      return matrix[0]; 		// Determinant égal à sa valeur si la matrice est de taille 1
    } 
    else {
        for (int i =0 ; i < size ; i++) {		// Sinon on calcule le déterminant via la sousmatrice
            coeff = (i %2 == 0) ? 1 : -1;
            det += coeff * matrix[i*size] * CalculDeterminant(SousMatrice(matrix, i, 0, size), size-1);
        }
        return det;
    }
}

double * Comatrice(double* flt_mat, int size) {
    double * coma;
    double coeff;
    coma =	(double *)malloc(sizeof(double) * size * size);	// On crée de l'espace
    for (int i = 0; i < size ; i++) {
        for (int j = 0; j < size ; j++) {
            coeff = ((i+j) % 2 == 0) ? 1 : -1;		// Calcul des coefficients 
            coma[i*size+j] = coeff * CalculDeterminant(SousMatrice (flt_mat, i, j, size), size - 1);		// Le résultat renvoie la commatrice
        }
    }
    return coma;
}

double *inverse(double *matrix, int size)
{
    double *inv = (double *)malloc(sizeof(double) * size * size);
    double *coma = (double *)malloc(sizeof(double) * size * size);
    double det;
    det = CalculDeterminant(matrix,size);
    coma = Comatrice(matrix,size);
    coma = Transposee(coma,size);

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            inv[i*size+j] = (1/det) * coma[i*size+j];
        }
    }

    return inv;

}

double *MatriceMatrice(double *matrix1, double *matrix2, int size) {
    double *res = (double *)malloc(sizeof(double) * size * size);	// Définition des entiers
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

    //double *Q_1 = (double *)malloc(sizeof(double) * size * size);
    //double *R_1 = (double *)malloc(sizeof(double) * size *size);
    //double *prod1 = (double *)malloc(sizeof(double) * size *size);

    for (int i=0; i<size; i++) {
        double doo = 0;
        for (int j=0; j<size; j++) {
            doo = doo + pow(A[j*size+i],2);
        }
        R[i*size+i] = sqrt(doo);
        for (int j=0; j<size; j++) {
            Q[j*size+i] = A[j*size+i]/R[i*size+i];
        }
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

double *Computation(double *A, double *B, int size) {

    double *Q = (double *)malloc(sizeof(double) * size * size);
    double *R = (double *)malloc(sizeof(double) * size * size);
    double *X = (double *)malloc(sizeof(double) * size *size);

    QR(A, size, Q, R);
    double *Y = MatriceMatrice(Transposee(Q, size), B, size);

    for (int i=size-1; i>-1;i--) {
        for (int j=size-1; j>-1;j--) {
            double s = 0;
            for (int k=i+1;k<size;k++) {
                s += R[i*size+k]*X[k*size+j];
            }
        X[i*size+j] = (Y[i*size+j]-s)/R[i*(size+1)];
        }
    }
    return X;
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
        double *aref;
        double *B = (double *)malloc(sizeof(double) * size * size);
        double *bref;

        A = generate_matrix(size);
        aref = generate_matrix(size);        
        B = generate_matrix(size);
        bref = generate_matrix(size);


        

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
        B = Computation(A, B, size);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        //if (check_result(bref,B,size)==1)
        //    printf("Result is ok!\n");
        //else    
        //    printf("Result is wrong!\n");
        
        //print_matrix("X", b, size);
        //print_matrix("Xref", bref, size);
        return 0;
    }
