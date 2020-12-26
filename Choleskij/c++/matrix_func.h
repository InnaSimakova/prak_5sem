using namespace std;

typedef double* Matrix;

Matrix Matrix_Initialization(int n);
Matrix Create_Matrix(int n, int a, int b, int seed=42);
void Print_Matrix(Matrix A, int n);
double Determinant (Matrix A, int n);
Matrix Matrix_Multiplication(Matrix A, Matrix B, int n);
Matrix Matrix_Transp(Matrix A, int n);
double Error_Rate(Matrix A, Matrix LU, int n);
Matrix manualInput(int n);
Matrix randomInput(int n);
Matrix Choleskij_Decomposition(Matrix A, int n);
int Check (Matrix A, int n);
