#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Matrix.h"

using namespace std;

void readFile(Matrix& a, Matrix& b) {
    ifstream input;
    bool isFirstNum = true;
    int n = 0;
    int counter = 0;
    double value;
    input.open("q5data.txt");

    if (input.is_open()) {
        cout << "input opened" << endl;
        while (!input.eof()) {
            int num;
            if (isFirstNum) {
                input >> n;
                a.setRows(n);
                a.setCols(n);
                isFirstNum = false;
                cout << "N: " << n << endl;
            } else {
                // the next n^2 numbers will fill the A matrix
                if (counter < (n * n)) {
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            input >> value;
                            //cout << value << endl;
                            a.setValue(i, j, value);
                            //cout << a(i, j) << endl;
                            counter++;
                        }
                    }
                } else {
                    // the next n values will fill the B matrix
                    b.setRows(n);
                    b.setCols(1);
                    for (int i = 0; i < n; i++) {
                        input >> value;
                        b.setValue(i, 0, value);
                    }
                    input.close();
                }
            }
        }
    }
}

void readFile(Matrix& a) {
    ifstream input;
    double value;
    int n = 4;

    input.open("q8data.txt");

    if (input.is_open()) {
        while (!input.eof()) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    input >> value;
                    a.setValue(i, j, value);
                }
            }
        }
        input.close();
    }
}

void writeSolutionToFile(const Matrix& soln) {
    ofstream stream;
    stream.open("q8soln.txt");

    if (stream.is_open()) {
        stream << setprecision(64) << soln;
/*        int counter = 0;
        int rootN = sqrt(soln.getRows());
        cout << rootN << endl;

        cout << "Open" << endl;
        for (int i = 0; i < soln.getRows(); i++) {
            for (int j = 0; j < soln.getCols(); j++) {
                stream << soln(i, j) << " ";
                counter++;
                if (counter == rootN) {
                    counter = 0;
                    stream << endl;
                }
            }
        }*/
        stream.close();
    }
}

/*
int main() {
    //Matrix A, b, L, U, x, y;
    //readFile(A, b);
    //cout << "Matrix A:\n" << A << endl;
    //cout << "Matrix b:\n" << b << endl;
   */
/* L.setRows(A.getRows());
    L.setCols(A.getCols());
    U = L;
    A.decompose(L, U);
    //cout << "Matrix L:\n" << L << endl;
    //cout << "Matrix U:\n" << U << endl;
    x.setRows(b.getRows());
    x.setCols(b.getCols());
    y = x;
    y = L.forwardSolve(b);
    //cout << "y:\n" << y << endl;
    x = U.backSolve(y);
    cout << "x:\n" << x << endl;*//*

    //cout << "Magnitude of the residual ||b - Ax||: " << A.calculateResidual(x, b) << endl;
    //x = A.jacobi(b, 1000);

*/
/*    double lambda = A.powerIteration(b);
    cout << "Lamba: " << lambda << endl;*//*


*/
/*    x = A.gaussSeidel(b, 10);
    cout << "x: \n" << x;

    cout << A * x << endl;

    cout << endl << endl;

    x = A.sor(b, 10, 1.2);
    cout << "x: \n" << x;*//*

*/
/*
    A.augmentMatrix(b);
    A.rowEchelonForm();
    x = A.backSolve();
    cout << endl << endl << "************\n" << "X: " << x << endl;*//*


    //x = A.jacobi(b, 1000);
    //x = A.gaussSeidel(b, 1000);
    //x = A.sor(b, 1000, 1.10);
    //cout << x << endl;
    Matrix A(4, 4);
    readFile(A);
    Matrix b(4, 1);
    Matrix x = b;

    for (int i = 0; i < b.getRows(); i++) {
        b.setValue(i, 0, 1000);
    }

    cout << "A: \n" << A << "\n\nb: \n" << b << endl;
    double eigenvalue = A.powerIteration(b, 1000);
    cout << b << endl;
    //cout << "A * b = " << (A * b) << endl;
    //cout << "lamda * b = " << (b * eigenvalue) << endl;

    writeSolutionToFile(b);
    return 0;
}*/
