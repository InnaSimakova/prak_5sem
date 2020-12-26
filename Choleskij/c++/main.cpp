#include <iostream>
#include <vector>
#include <iomanip>
#include <cstring>
#include <cctype>
#include <ctime>
#include "matrix_func.h"

using namespace std;

int main() {
	int n;
	cout << "Enter square matrix (n x n) dimension:" << endl << "n = ";
    	cin >> n;

    	char manual_mode;
    	Matrix A;
    	while (1) {
        	cout << "Manual input of matrix? (Y/n)" << endl;
        	cin >> manual_mode;
        	manual_mode = toupper(manual_mode);
        	if (manual_mode == 'Y') {
            		A = manualInput(n);
            		break;
        	} 
		else if (manual_mode == 'N') {
            		A = randomInput(n);
            		break;
        	} 
		else {
            		cout << "Answer not recognized, try again." << endl;
        	}
    	}
    	cout << "A:" << endl;
    	Print_Matrix(A,n);
    	cout<<endl;

	Matrix L,Lt,LLt;
	
    	if (Check(A,n) == 1) {
        	cout << "Matrix A is not suitible for our metod." << endl;
    	} 
    	else {
		int start = clock();
               	L = Choleskij_Decomposition(A, n);
        	int end = clock();
		
		cout << "L:" <<endl;
		Print_Matrix(L,n);
		cout << "Lt:" << endl;
		Lt = Matrix_Transp(L,n);
		Print_Matrix(Lt,n);
		LLt = Matrix_Multiplication(L, Lt, n);
		cout << "LLt:" << endl;
		Print_Matrix(LLt, n);
        
        	cout << "Time: " << (end - start) << " ms" << endl;
        	cout << "||LLt - A|| = " << Error_Rate(A, LLt, n) << endl;
    	}
    	return 0;
}
