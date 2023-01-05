/*
import all the libraries that are useful
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <time.h>

/*
define macros which are constant through the code
*/

#define MAX_THREADS 65536

/*
define your global variables that can access most functions
*/

int k, base, threads_num;

/*
size is the size of the matrix
Allocate the memory for matrix dynamically using new
Allocate n^2 memory locations inorder to store the elements of matrix
*/

int** Allocmatrix(int size){
    //obtain the pointer for the first element in the row to start iterating
    int **ptr = new int*[size];
    int *row = new int[size*size];

    for(int i = 0; i < size; i++){
        *(ptr+i) = row;
        row+= size;
    }
    return ptr;
}

//freeup the memory allocated to the matrix
void deAlloc_mat(int **ptr1){
    //deallocate the memory allocated for the elements in the matrix
    delete []*ptr1;
    //delete the row pointers as well
    delete []ptr1;
}


/*
Always divide the matrix into 4 parts for recursive divide and conquer
submat is the one-fourth part of the matrix
mat is the original matrix which should be divided
m,n are variables that identify the quarter in the matrix
*/

void divMatrix(int **submat, int size, int **mat, int m, int n){
    for(int i = 0; i < size; i++){
       submat[i] = &mat[m+i][n];
    }
}

/*
Compute matrix addition 
Input matrices A & B 
Result is stored in submat
size of matrix is stored in size
*/
void addMatrix(int **submat, int size, int **A, int **B){
    for (int i = 0; i < size ; i++)
		for (int j = 0; j < size; j++)
			submat[i][j] = A[i][j] + B[i][j];
}

/*
Compute matrix subtraction
Input matrices A & B 
Result is stored in submat
size of matrix is stored in size
*/
void subMatrix(int **submat, int size, int **A, int **B){
    for (int i = 0; i < size ; i++)
		for (int j = 0; j < size; j++)
			submat[i][j] = A[i][j] - B[i][j];
}


/*
compute the resultant matrix from the naive matrix multiplication of input matrices A and B
C is the resultant matrix
A and B are input matrices
size is the size of the matrix
*/

void basicmatrixmul(int size, int **A, int **B, int **C){
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            C[i][j] = 0;
            for(int k = 0; k < size; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
/*
Strassen algo employs a recursive divide and conquer with reduced multiplications
A,B are input matrices
C is the matrix obtained after computing the product
Size is the size of the matrix
*/
void Strassen(int size, int **A, int **B, int** C){

    //if size doesnot cross the threshold compute naive matrix multiplication
    if(((float)size) <= pow(2,k-base)){
        for(int i = 0; i<size; i++){
			for(int j=0; j<size; j++){
				C[i][j] =0;
				for(int k=0; k<size; k++)
					C[i][j] += A[i][k]*B[k][j];
			}
		}
    }
    //if the size crosses our threshold employ strassen
    else
    {   
        //use the split to divide the matrix into four parts
        int split = size/2;

        //Allocate memory for the strassen multiplication matrices
        int **M1 = Allocmatrix(split);
        int **M2 = Allocmatrix(split);
        int **M3 = Allocmatrix(split);
        int **M4 = Allocmatrix(split);
        int **M5 = Allocmatrix(split);
        int **M6 = Allocmatrix(split);
        int **M7 = Allocmatrix(split);

        //initialize partial matrices to store intermediate addition and subtraction results

        int **subM1 = Allocmatrix(split);
        int **subM2 = Allocmatrix(split);
        int **subM3 = Allocmatrix(split);
        int **subM4 = Allocmatrix(split);
        int **subM5 = Allocmatrix(split);
        int **subM6 = Allocmatrix(split);
        int **subM7 = Allocmatrix(split);
        int **subM8 = Allocmatrix(split);
        int **subM9 = Allocmatrix(split);
        int **subM10 = Allocmatrix(split); 

        //Allocate memory for elements of A
        int **A11 = new int*[split];
		int **A12 = new int*[split];
		int **A21 = new int*[split];
		int **A22 = new int*[split];

        //Allocate memory for elements of B
        int **B11 = new int*[split];
		int **B12 = new int*[split];
		int **B21 = new int*[split];
		int **B22 = new int*[split];

        //Allocate memory for elements of C
        int **C11 = new int*[split];
		int **C12 = new int*[split];
		int **C21 = new int*[split];
		int **C22 = new int*[split];

        // Split the matrix into 4 parts for each matrix 
        divMatrix(A11, split, A, 0, 0);
		divMatrix(A12, split, A, 0, split);
		divMatrix(A21, split, A, split, 0);
		divMatrix(A22, split, A, split, split);

        //Split the matrix into 4 parts for each matrix
        divMatrix(B11, split, B, 0, 0);
		divMatrix(B12, split, B, 0, split);
		divMatrix(B21, split, B, split, 0);
		divMatrix(B22, split, B, split, split);

        //Split the matrix into 4 parts for each matrix
        divMatrix(C11, split, C, 0, 0);
		divMatrix(C12, split, C, 0, split);
		divMatrix(C21, split, C, split, 0);
		divMatrix(C22, split, C, split, split);

        //Employing parallelism in computing M1-7

        // M1 = (A11 + A22) * (B11 + B22)
        #pragma omp task
		{
			addMatrix(subM1, split, A11, A22);
			addMatrix(subM2, split, B11, B22);
			Strassen(split, subM1, subM2, M1);
		}

        // M2 = (A21 + A22) * B11
        #pragma omp task
		{
			addMatrix(subM3, split, A21, A22);
			Strassen(split, subM3, B11, M2);
		}

        // M3 = A11 * (B12 - B22)
        #pragma omp task
		{
			subMatrix(subM4, split, B12, B22);
			Strassen(split, A11, subM4, M3);
		}

        //M4 = A22 * (B21 - B11)
        #pragma omp task
		{
			subMatrix(subM5, split, B21, B11);
			Strassen(split, A22, subM5, M4);
		}
        
        //M5 = (A11 + A12) * B22
        #pragma omp task
		{
			addMatrix(subM6, split, A11, A12);
			Strassen(split, subM6, B22, M5);
		}

        // M6 = (A21 - A11) * (B11 + B12)
        #pragma omp task
		{
			subMatrix(subM7, split, A21, A11);
			addMatrix(subM8, split, B11, B12);
			Strassen(split, subM7, subM8, M6);
		}

        //M7 = (A12 - A22) * (B21 + B22)
         #pragma omp task
		{
			subMatrix(subM9, split, A12, A22);
			addMatrix(subM10, split, B21, B22);
			Strassen(split, subM9, subM10, M7);
		}

        // wait for all threads to finish their execution

        #pragma omp taskwait

        //calculate the final result of strassen based on M1-M7

        for(int i=0; i<split; i++){
			for(int j=0; j<split; j++){
				C11[i][j] = M1[i][j]+M4[i][j]-M5[i][j]+M7[i][j];
				C12[i][j] = M3[i][j]+M5[i][j];
				C21[i][j] = M2[i][j] + M4[i][j];
				C22[i][j] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
			}
		} 

        // Deallocate memory and perform garbage collection

        deAlloc_mat(subM1);
        deAlloc_mat(subM2);
        deAlloc_mat(subM3);
        deAlloc_mat(subM4);
        deAlloc_mat(subM5);
        deAlloc_mat(subM6);
        deAlloc_mat(subM7);
        deAlloc_mat(subM8);
        deAlloc_mat(subM9);
        deAlloc_mat(subM10);
        deAlloc_mat(M1);
        deAlloc_mat(M2);
        deAlloc_mat(M3);
        deAlloc_mat(M4);
        deAlloc_mat(M5);
        deAlloc_mat(M6);
        deAlloc_mat(M7);

        //deallocate memory for the elements of A
        delete[] A11;
        delete[] A12;
        delete[] A21;
        delete[] A22;

        //deallocate memory for the elements of B
        delete[] B11;
        delete[] B12;
        delete[] B21;
        delete[] B22;

        //deallocate memory for the elements of C
        delete[] C11;
        delete[] C12;
        delete[] C21;
        delete[] C22;

    }
}


int main(int argc, char* argv[]){

    //use k to calculate the matrix size for n*n matrix i.e n = 2^k
    k = atoi(argv[1]);
    int size = pow(2,k);
    //obtain your kdash threshold for base case
    base = atoi(argv[2]);
    int power = atoi(argv[3]);
    threads_num = (1 << power);

    if(threads_num > MAX_THREADS){
        printf("\n Maximum limits of threads exceeded error : %d.", MAX_THREADS);
        exit(0);
    }

    //set number of threads
    omp_set_dynamic(0);
    omp_set_num_threads(threads_num);


    //Allocate the memory for the matrix
    int **A = Allocmatrix(size);
    int **B = Allocmatrix(size);
    int **C = Allocmatrix(size);
    int **result = Allocmatrix(size);

    //Assign some random values for defining the matrix
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            A[i][j] = rand()%100;
            B[i][j] = rand()%100;
        }
    }
    //cross verify the result with normal multiplication
    basicmatrixmul(size, A, B, result);

    //calculate run time
    struct timespec start, stop;
    double total_time;
    clock_gettime(CLOCK_REALTIME, &start);

    //Code for Strassen multiplication

    #pragma omp parallel
    {
        #pragma omp single
        {
            Strassen(size, A, B,C);
        }
    }
    //calculating execution time
    clock_gettime(CLOCK_REALTIME, &stop);
    total_time = (stop.tv_sec-start.tv_sec)+0.000000001*(stop.tv_nsec-start.tv_nsec);
    
    //checking the correctness of the result
    bool valid = true;
    for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if(result[i][j] != C[i][j]) valid = false;
    if(valid){
        printf("Correct: k=%d, Size of Matrix = %d X %d, base= %d, number of threads = %d, Execution time for strassen = %8.4f sec \n",  k,size,size, base,threads_num,total_time );          
	}
    else{
        printf("Incorrect: k=%d, Size of Matrix = %d X %d, base= %d, number of threads = %d, Execution time for strassen = %8.4f sec \n",  k,size,size, base,threads_num,total_time );          
    }

    //perform garbage collection

    deAlloc_mat(A);
    deAlloc_mat(B);
    deAlloc_mat(C);
    deAlloc_mat(result);

    return 0;

}