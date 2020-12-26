#include <iostream>
#include <vector>
#include <iomanip>
#include <cstring>
#include <cmath>

#include "matrix_func.h"

using namespace std;

int Check (Matrix A, int n) {
	int flag = 0;
	if (Determinant(A,n)<=0) {
		flag = 1;
		return flag;
	}
    	for (int i=1; i<n; i++) {
		Matrix B=Matrix_Initialization(i);
		for (int j=0;j<i;j++)
			for (int k=0;k<i;k++)
				B[j*i+k]=A[j*n+k];
		if (Determinant(B,i)<0) flag = 1;
	}
        return flag;
}

double Determinant (Matrix A, int n) {
	int det;
  	Matrix B = Matrix_Initialization(n-1);
  	det = 0;
  	int k = 1; //(-1) в степени i

  	if (n == 1) det = A[0];
  	if (n == 2) det = A[0] * A[3] - (A[1] * A[2]);
 	if (n>2) {
    		for (int i = 0; i<n; i++) {
			int j=0;
  			int i1, j1, di, dj;
  			di = 0;
  			for (i1 = 0; i1<n-1; i1++) { // проверка индекса строки
    				if (i1 == i) di = 1;
   	 			dj = 0;
    				for (j1 = 0; j1<n-1; j1++) { // проверка индекса столбца
      					if (j1 == j) dj = 1;
      					B[i1*n+j1] = A[(i1 + di)*n+(j1 + dj)];
    				}
  			}
      			det = det + k * A[i*n] * Determinant(B,n-1);
      			k = -k;
		}
	}
	return det;
}

Matrix Matrix_Initialization(int n){
    double *A = new double[n * n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i*n+j] = 0.0;
        }
    }
    return A;
}

Matrix Create_Matrix(int n, int a, int b, int seed){
    srand(seed);
    double *mat = new double[n * n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat[i*n+j] = a + (rand() % static_cast<int>(b - a + 1));
        }
    }
    return mat;
}

Matrix Matrix_Multiplication(Matrix A, Matrix B, int n) {
    Matrix C = Matrix_Initialization(n);
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                C[i+n*j] += A[k*n+i] * B[j*n+k];
            }
        }
    }
    return C;
}

Matrix Matrix_Transp(Matrix A, int n) {
	Matrix At = Matrix_Initialization(n);
    	for (int i = 0; i < n; ++i) {
        	for (int j = 0; j < n; ++j) {
            		At[i*n+j] = A[j*n+i];
        	}
    	}
	return At;
}

void Print_Matrix(Matrix A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j != 0) {
                cout << " ";
            }
            cout << fixed << setiosflags(ios::left) << setw(6)
                      << setfill(' ') << setprecision(2) << A[i*n+j];
        }
        cout << endl;
    }
}

Matrix Choleskij_Decomposition(Matrix A, int n) {
    Matrix L = Matrix_Initialization(n);
    for (int i = 0; i < n; i++) {
        double temp;
        //Сначала вычисляем значения элементов слева от диагонального элемента,
        //так как эти значения используются при вычислении диагонального элемента.
        for (int j = 0; j < i; j++) {
            temp = 0;
            for (int k = 0; k < j; k++) {
                temp += L[i*n+k] * L[j*n+k];
            }
            L[i*n+j] = (A[i*n+j] - temp) / L[j*n+j];
        }
        //Находим значение диагонального элемента
        temp = A[i*n+i];
        for (int k = 0; k < i; k++) {
            temp -= L[i*n+k] * L[i*n+k];
        }
        L[i*n+i] = sqrt(temp);
    }
    return L;
}

double Error_Rate(Matrix A, Matrix LU, int n) {
    Matrix dA = Matrix_Initialization(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dA[i*n+j] = A[i*n+j] - LU[i*n+j];
        }
    }

    double norm = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            norm += dA[i*n+j] * dA[i*n+j];
        }
    }
    return sqrt(norm);
}

Matrix manualInput(int n) {
    Matrix A = Matrix_Initialization(n);
    cout << "Enter square " << n << "-dimensional matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> A[i*n+j];
        }
    }
    return A;
}

Matrix randomInput(int n) {
    int a, b;
    cout << "Define range of values (a, b) for randomly generated square " << n <<"-dimensional matrix" << endl;
    cout << "a = ";
    cin >> a;
    cout << "b = ";
    cin >> b;
    return Create_Matrix(n, a, b);
}

