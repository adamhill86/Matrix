//
// Created by Adam on 4/20/2017.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Matrix.h"

using namespace std;

void buildMatrix(Matrix& A, int nx) {
    A.identityMatrix();
    //cout << A << endl;

    for (int k = nx + 1; k < 379; k++) {
        if ((k % 20 == 0) || (k == 39) || (k == 59) || (k == 79) || (k == 99) ||
            (k == 119) || (k == 139) || (k == 159) || (k == 179) || (k == 199) ||
            (k == 219) || (k == 239) || (k == 259) || (k == 279) || (k == 299) ||
            (k == 319) || (k == 339) || (k == 359)) {
            continue;
        } else {
            // values from stencil
            A.setValue(k, k - nx, 1.0);
            A.setValue(k, k - 1, 1.0);
            A.setValue(k, k, -4.0);
            A.setValue(k, k + 1, 1.0);
            A.setValue(k, k + nx, 1.0);
        }
    }

    ofstream stream;
    stream.open("matrixA.txt");
    int counter = 0;

    for (int i = 0; i < A.getRows(); i++) {
        for (int j = 0; j < A.getCols(); j++) {
            stream << A(i, j) << " ";
            counter++;

            if (counter == 400) {
                counter = 0;
                stream << endl;
            }
        }
    }

    stream.close();
}

void buildAMatrices(Matrix& Ai, Matrix& Aj, int nx) {
    Ai.identityMatrix();
    Aj = Ai;

    // 219 is the right boundary of row 10
    // This gives Ai rows 0 - 10 (11 total rows)
    for (int k = nx + 1; k < 219; k++) {
        if ((k % 20 == 0) || (k == 39) || (k == 59) || (k == 79) || (k == 99) ||
            (k == 119) || (k == 139) || (k == 159) || (k == 179) || (k == 199)) {
            continue;
        } else {
            // values from stencil
            Ai.setValue(k, k - nx, 1.0);
            Ai.setValue(k, k - 1, 1.0);
            Ai.setValue(k, k, -4.0);
            Ai.setValue(k, k + 1, 1.0);
            Ai.setValue(k, k + nx, 1.0);
        }
    }

/*    // Artificial boundary for Ai and Aj
    // set domain row 9 (matrix row 180) to 0 (initial random guess essentially)
    for (int j = 0; j < Ai.getCols(); j++) {
        cout << "Ai row 180:\n" << "j: " << j << " value: " << Ai(180, j) << " ";
        //Ai.setValue(181, j, 0.0);
        Aj.setValue(200 - 180, j, 0.0); // starts at 200 on grid but -180 to adjust for matrix indexing
    }*/

    for (int k = 181 + nx; k < 379; k++) {
        int row = k - 180;
        if ((k % 20 == 0) || (k == 199) ||
            (k == 219) || (k == 239) || (k == 259) || (k == 279) || (k == 299) ||
            (k == 319) || (k == 339) || (k == 359)) {
            continue;
        } else {
            Aj.setValue(row, row - nx, 1.0);
            Aj.setValue(row, row - 1, 1.0);
            Aj.setValue(row, row, -4.0);
            Aj.setValue(row, row + 1, 1.0);
            Aj.setValue(row, row + nx, 1.0);
        }
    }
}

void buildInitialBoundaryConditions(Matrix& b, int nx, double rhs) {
    b.setValue(19, 0, -100.0);

    for (int i = nx + 1; i < 379; i++) {
        if ((i % 20 == 0) || (i == 39) || (i == 59) || (i == 79) || (i == 99) ||
                (i == 119) || (i == 139) || (i == 159) || (i == 179) || (i == 199) ||
                (i == 219) || (i == 239) || (i == 259) || (i == 279) || (i == 299) ||
                (i == 319) || (i == 339) || (i == 359)) {
            continue;
        } else {
            b.setValue(i, 0, rhs);
        }
    }

    cout << b << endl;

    ofstream stream;
    stream.open("initalBC.txt");
    int counter = 0;

    for (int i = 0; i < b.getRows(); i++) {
        stream << b(i, 0) << " ";
        counter++;

        if (counter == 20) {
            counter = 0;
            stream << endl;
        }
    }

    stream.close();
}

void buildInitialSplitBoundaryConditions(Matrix& bi, Matrix& bj, int nx, double rhs) {
    bi.setValue(19, 0, -100.0);

    for (int i = nx + 1; i < 219; i++) {
        if ((i % 20 == 0) || (i == 39) || (i == 59) || (i == 79) || (i == 99) ||
            (i == 119) || (i == 139) || (i == 159) || (i == 179) || (i == 199)) {
            continue;
        } else {
            bi.setValue(i, 0, rhs);
        }
    }

    for (int i = 181 + nx; i < 379; i++) {
        int row = i - 180;
        if ((i % 20 == 0) || (i == 199) ||
            (i == 219) || (i == 239) || (i == 259) || (i == 279) || (i == 299) ||
            (i == 319) || (i == 339) || (i == 359)) {
            continue;
        } else {
            bj.setValue(row, 0, rhs);
        }
    }

    // artificial boundary conditions
    for (int i = 181; i < 199; i++) {
        bi.setValue(i, 0, 0.0);
    }

    for (int i = 201; i < 219; i++) {
        bj.setValue(i - 180, 0, 0.0);
    }
}

void writeSolutionToFile(const Matrix& soln, string filename) {
    ofstream stream;
    stream.open(filename);

    if (stream.is_open()) {
        stream << soln;
        stream.close();
    }
}


int main() {
    Matrix Ai(220, 220), Aj(220, 220);
    Matrix bi(220, 1), bj(220, 1);
    Matrix Li(220, 220), Lj(220, 220), Ui(220, 220), Uj(220, 220);
    Matrix xi(220, 1), xj(220, 1), yi(220, 1), yj(220, 1);
    int nx = 20, ny = nx;
    double h = 1.0 / (nx - 1);
    double rhs = -4 * (h * h);
    buildAMatrices(Ai, Aj, nx);
    buildInitialSplitBoundaryConditions(bi, bj, nx, rhs);


    Ai.decompose(Li, Ui);
    yi = Li.forwardSolve(bi);
    xi = Ui.backSolve(yi);

    Aj.decompose(Lj, Uj);
    yj = Lj.forwardSolve(bj);
    xj = Uj.backSolve(yj);

    cout << "Residual for Ai: " << Ai.calculateResidual(xi, bi) << endl;
    cout << "Redisual for Aj: " << Aj.calculateResidual(xj, bj) << endl;

    int iterations = 50;
    while (iterations > 0) {
        // get row 9 values from bottom of Aj's solution vector x
        for (int i = 181; i < 199; i++) {
            bi.setValue(i, 0, xj(i - 181, 0));
        }

        // get row 10 values from top of Ai's solution vector x
        for (int i = 201; i < 219; i++) {
            bj.setValue(i - 180, 0, xi(i, 0));
        }

        yi = Li.forwardSolve(bi);
        xi = Ui.backSolve(yi);
        yj = Lj.forwardSolve(bj);
        xj = Uj.backSolve(yj);

        cout << "Residual for Ai: " << Ai.calculateResidual(xi, bi) << endl;
        cout << "Redisual for Aj: " << Aj.calculateResidual(xj, bj) << endl;

        iterations--;
    }

    writeSolutionToFile(xi, "xi_soln.txt");
    writeSolutionToFile(xj, "xj_soln.txt");

    return 0;
}
