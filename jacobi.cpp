//
// Created by Adam on 4/20/2017.
//
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
    input.open("data.txt");

    if (input.is_open()) {
        //cout << "input opened" << endl;
        while (!input.eof()) {
            int num;
            if (isFirstNum) {
                input >> n;
                a.setRows(n);
                a.setCols(n);
                isFirstNum = false;
                //cout << "N: " << n << endl;
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
int main() {
    Matrix A, b, x;
    readFile(A, b);
    x = A.jacobi(b, 1000);
    cout << x << endl;
    cout << "Magnitude of the residual ||b - Ax||: " << A.calculateResidual(x, b) << endl;
    return 0;
}
