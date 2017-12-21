//
// Created by Adam on 1/13/2017.
//

#include <random>
#include <ctime>
#include <iomanip>
#include <fstream>
#include "Matrix.h"

Matrix::Matrix() {
    rows = 0;
    cols = 0;
    M = nullptr;
}

Matrix::Matrix(int r, int c) {
    rows = r;
    cols = c;
    M = new double*[r];

    for (int i = 0; i < r; i++) {
        M[i] = new double[c];
    }

    // fill new matrix with 0s
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            setValue(i, j, 0.0);
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < rows; i++) {
        delete M[i];
    }

    delete[] M;
}

void Matrix::setValue(int r, int c, double value) {
    M[r][c] = value;
}

int Matrix::getCols() const {
    return cols;
}

int Matrix::getRows() const {
    return rows;
}

void Matrix::setCols(int c) {
    if (c >= 0) {
        cols = c;
    }
    // TODO: actually increase the size in the M
    for (int i = 0; i < rows; i++) {
        M[i] = new double[c];
    }
}

void Matrix::setRows(int r) {
    if (r >= 0) {
        rows = r;
    }
    // TODO: actually increase the size in the M
    M = new double*[r];
}

ostream& operator<<(ostream& os, const Matrix& matrix) {
    for (int i = 0; i < matrix.getRows(); i++) {
        for (int j = 0; j < matrix.getCols(); j++) {
            //os << setprecision(5) << setw(10) << matrix(i,j);
            os << matrix(i, j);
        }
        os << endl;
    }

    return os;
}

/**
 * Access the matrix element from row r, column c.
 */
double Matrix::operator()(int r, int c) const {
    return M[r][c];
}

Matrix::Matrix(const Matrix &rhs) {
    rows = rhs.getRows();
    cols = rhs.getCols();

    M = new double*[rows];
    for (int i = 0; i < rows; i++) {
        M[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            M[i][j] = rhs(i, j);
        }
    }
}

void Matrix::operator=(const Matrix &rhs) {
    rows = rhs.getRows();
    cols = rhs.getCols();

    M = new double*[rows];
    for (int i = 0; i < rows; i++) {
        M[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            M[i][j] = rhs(i, j);
        }
    }
}

/**
 * Add two matrices together. Only works if both matrices have the same dimensions.
 * @return The resultant matrix; null if incompatible sizes.
 */
Matrix Matrix::operator+(const Matrix &rhs) {
    Matrix result;

    if (rhs.getCols() == cols && rhs.getRows() == rows) {
        result.setRows(rows);
        result.setCols(cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.setValue(i, j, M[i][j] + rhs(i, j));
            }
        }
    }

    return result;
}

/**
 * Scalar multiplication
 * @param scalar
 * @return
 */
Matrix Matrix::operator*(double scalar) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            M[i][j] *= scalar;
        }
    }
    return *this;
}

/**
 * Matrix multiplication. Matrix A with dims m x n, matrix B with dims n x p. n must = n.
 * Resulting matrix will be m x p.
 * @param rhs
 * @return The left-hand side matrix if the calculation cannot be done.
 */
Matrix Matrix::operator*(const Matrix &rhs) {
    if (cols == rhs.getRows()) {
        Matrix result(rows, rhs.getCols());
        // TODO: the actual calculations

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rhs.getCols(); j++) {
                for (int k = 0; k < cols; k++) {
                    result.setValue(i, j, result(i, j) + (M[i][k] * rhs(k, j)));
                }
            }
        }

        return result;
    }
    return *this;
}

Matrix Matrix::operator/(double scalar) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            M[i][j] /= scalar;
        }
    }
    return *this;
}

void Matrix::identityMatrix() {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                M[i][j] = 1.0;
            } else {
                M[i][j] = 0.0;
            }
        }
    }
}

/**
 * Fills the matrix with random values between -max and max.
 * The values on the main diagonal will be the largest value in each row
 * @param max The value used to determine the max and min random values to be generated.
 */
void Matrix::randomMatrix(double max) {
    double sum = 0.0;
    double randomNumber;

    default_random_engine gen(time(NULL));
    uniform_real_distribution<double> rdist(-max, max);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            randomNumber = rdist(gen);
            M[i][j] = randomNumber;
            sum += abs(randomNumber);
        }
        M[i][i] = sum;
        sum = 0.0;
    }
}

void Matrix::randomIntMatrix(int max) {
    int sum = 0;
    int randomNumber;

    default_random_engine gen(time(NULL));
    uniform_int_distribution<int> rdist(-max, max);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            randomNumber = rdist(gen);
            M[i][j] = randomNumber;
            sum += 10 * abs(randomNumber); // make the diagonal VERY dominant
        }
        M[i][i] = sum;
        sum = 0;
    }
}

void Matrix::randomIntVector(int max) {
    int randomNumber;

    default_random_engine gen(time(NULL));
    uniform_int_distribution<int> rdist(-max, max);

    for (int i = 0; i < rows; i++) {
        randomNumber = rdist(gen);
        M[i][0] = randomNumber;
    }
}

Matrix Matrix::operator-(const Matrix &rhs) {
    Matrix result;

    if (rhs.getCols() == cols && rhs.getRows() == rows) {
        result.setRows(rows);
        result.setCols(cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.setValue(i, j, M[i][j] - rhs(i, j));
            }
        }
    }

    return result;
}

/**
 * Pre: M must have been an n x n square matrix before augmentation
 */
void Matrix::rowEchelonForm() {
    int n = cols - 1; // subtract 1 because the matrix is augmented

    for (int i = 0; i < n; i++) {
        // find largest value in current column
        double max = M[i][i];
        int maxRow = i;

        for (int row = i + 1; row < n; row++) {
            if (M[row][i] > max) {
                max = M[row][i];
                maxRow = row;
            }
        }

        // Swap rows if necessary
        for (int k = i; k < n + 1; k++) {
            double temp = M[maxRow][k];
            M[maxRow][k] = M[i][k];
            M[i][k] = temp;
        }

        // set all values to 0 below this row in current column
        for (int k = i + 1; k < rows; k++) {
            double c = -M[k][i] / M[i][i];

            for (int j = i; j < n + 1; j++) {
                if (i == j) {
                    M[k][j] = 0.0;
                } else {
                    M[k][j] += c * M[i][j];
                }
            }
        }
        // divide through by leading term
        double leadingTerm = M[i][i];

        for (int k = i; k < cols; k++) {
            M[i][k] = M[i][k] / leadingTerm;
        }
    }
}

/**
 * Pre: the rhs matrix must be a vector of n elements
 * @param rhs
 */
void Matrix::augmentMatrix(const Matrix &rhs) {
    Matrix temp = *this;

    // erase the current matrix
    for (int i = 0; i < rows; i++) {
        delete M[i];
    }

    cols++;

    for (int i = 0; i < rows; i++) {
        M[i] = new double[cols];
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (j == cols - 1) {
                M[i][j] = rhs(i, 0);
            } else {
                M[i][j] = temp(i, j);
            }
        }
    }
}

Matrix Matrix::backSolve() {
    Matrix result(rows, 1);
    int n = cols - 1;
    double answer = 0.0;

    for (int i = n - 1; i >= 0; i--) {
        answer = M[i][n] / M[i][i];
        result.setValue(i, 0, answer);
        for (int k = i - 1; k >= 0; k--) {
            double value = M[k][i] * result(i, 0);
            M[k][n] -= value;
        }
    }

    return result;
}

/**
 *
 * @param calculatedSolution Our calculated matrix x from Ax = b
 * @param rhs matrix b
 * @return
 */
double Matrix::calculateResidual(Matrix& calculatedSolution, Matrix& rhs) {
    Matrix errorVector(rhs.getRows(), rhs.getCols());
    Matrix ax(rhs.getRows(), rhs.getCols());

    ax = (*this) * calculatedSolution;
    errorVector = rhs - ax;


    double sum = 0.0;
    for (int i = 0; i < errorVector.getRows(); i++) {
        for (int j = 0; j <  errorVector.getCols(); j++) {
            sum += pow(errorVector(i, j), 2.0);
        }
    }

    return sqrt(sum);
}

void Matrix::decompose(Matrix& L, Matrix& U) {
    int n = rows;

    for (int j = 0; j < n; j++) {
        L.setValue(j, j, 1.0);

        for (int i = 0; i < j + 1; i++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += U(k, j) * L(i, k);
                //cout << "Sum: " << sum << endl;
            }
            double value = M[i][j] - sum;
            U.setValue(i, j, value);
        }
        for (int i = j; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += U(k, j) * L(i, k);
            }
            L.setValue(i, j, (M[i][j] - sum) / U(j, j));
        }
    }
    //cout << L * U << endl;
}

Matrix Matrix::forwardSolve(const Matrix& b) {
    Matrix result(rows, 1);
    int n = rows;
    double sum;

    for (int i = 0; i < n; i++) {
        sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += M[i][j] * result(j, 0);
        }
        if (M[i][i] == 0) {
            cout << "***********Forward solve division by 0*********\n";
        }
        double answer = (b(i, 0) - sum) / M[i][i];
        result.setValue(i, 0, answer);
    }
    return result;
}

Matrix Matrix::backSolve(const Matrix &b) {
    Matrix result(rows, 1);
    int n = rows;
    double sum;

    for (int i = n - 1; i >= 0; i--) {
        sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += M[i][j] * result(j, 0);
        }
        if (M[i][i] == 0) {
            cout << "***********Backsolve division by 0*********\n";
        }
        double answer = (b(i, 0) - sum) / M[i][i];
        result.setValue(i, 0, answer);
    }
    return result;
}

Matrix Matrix::jacobi(Matrix &b, int numIterations) {
    Matrix dInverse(rows, cols);
    Matrix C = *this;
    Matrix xOld(b.getRows(), b.getCols());
    Matrix xNew = xOld;
    xNew.randomIntVector(5);
    double previous = 0.0;
    int counter = 1;

    // form matrices D and C
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                dInverse.setValue(i, j, 1 / M[i][j]);
                C.setValue(i, j, 0.0);
            }
        }
    }

    //cout << "D inverse:\n" << dInverse << endl;
    //cout << "C:\n" << C << endl;
    //cout << "xNew:\n" << xNew << endl;

    for (int i = 0; i < numIterations; i++) {
        xNew = dInverse * (b - (C * xOld));
        double residual = calculateResidual(xNew, b);

        if (previous != residual) {
            cout << "Residual: " << residual << ", iteration " << counter << endl;
            previous = residual;
            counter++;
        } else {
            cout << "Residual: " << residual << ", iteration " << counter << endl;
            xOld = xNew;
            break;
        }
        xOld = xNew;
        //cout << "xNew iteration " << i << endl << xNew << endl;
    }
    return xNew;
}

// Calculates the magnitude of a 1D Matrix (vector)
double Matrix::magnitude() {
    double square = 0.0;

    for (int i = 0; i < cols; i++) {
        square += (M[i][i] * M[i][i]);
    }

    return sqrt(square);
}

double Matrix::powerIteration(Matrix &b, int iterations) {
    double norm = 0.0;
    double previous = norm;

    ofstream stream;
    stream.open("q8-lamba-values.txt");

    for (int i = 0; i < iterations; i++) {
        Matrix temp = b;
        temp = (*this) * b;
        //cout << temp << endl;

        norm = temp.magnitude();
        cout << "Norm: " << setprecision(64) << norm << endl;
        stream << "Lamba: " << setprecision(64) << norm << endl;
        b = temp / norm;

        if (previous == norm) {
            cout << "Converged\n";
            break;
        } else {
            previous = norm;
        }
    }
    stream.close();
    return norm;
}

Matrix Matrix::gaussSeidel(Matrix &b, int iterations) {
    int n = rows;
    int counter = 1;
    double sum = 0.0;
    double residual, previous = 0.0;
    Matrix x(b.getRows(), b.getCols());
    //cout <<"x: \n" << x << endl;

    while (iterations > 0) {
        for (int i = 0; i < n; i++) {
            sum = b(i, 0);
            for (int j = 0; j < cols; j++) {
                if (i != j) {
                    sum -= (M[i][j] * x(j, 0));
                }
            }
            x.setValue(i, 0, sum / M[i][i]);
            //cout << "x[" << i << "]: " << x(i, 0) << endl;
        }
        residual = calculateResidual(x, b);
        cout << "Residual: " << residual << endl;
        if (residual == previous) {
            cout << "Converged after " << counter << " iterations\n";
            break;
        } else {
            previous = residual;
            counter++;
        }
        iterations--;
    }

    return x;
}

Matrix Matrix::sor(Matrix &b, int iterations, double omega) {
    /*int n = rows;
    int counter = 1;
    double sum = 0.0;
    double residual, previous = 0.0;
    Matrix x(b.getRows(), b.getCols());
    //cout <<"x: \n" << x << endl;

    while (iterations > 0) {
        for (int i = 0; i < n; i++) {
            sum = b(i, 0);
            for (int j = 0; j < cols; j++) {
                if (i != j) {
                    sum -= (M[i][j] * x(j, 0));
                }
            }
            x.setValue(i, 0, ((1 - omega) * x(i, 0)) + ((omega * sum) / M[i][i]));
            //cout << "x[" << i << "]: " << x(i, 0) << endl;
        }
        residual = calculateResidual(x, b);
        cout << "Residual: " << residual << endl;
        if (residual == previous) {
            cout << "Converged after " << counter << " iterations\n";
            break;
        } else {
            previous = residual;
            counter++;
        }
        iterations--;
    }
*/
    int n = rows;
    int counter = 1;
    double sum = 0.0;
    double residual, previous = 0.0;
    Matrix x(b.getRows(), b.getCols());

    while (iterations > 0) {
        for (int i = 0; i < n; i++) {
            sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += M[i][j] * x(j, 0);
                }
            }
            x.setValue(i, 0, (1 - omega) * x(i, 0) + ((omega / M[i][i]) * (b(i, 0) - sum)));
            //cout << "x[" << i << "]: " << x(i, 0) << endl;
        }
        residual = calculateResidual(x, b);
        cout << "Residual: " << residual << endl;
        if (residual == previous) {
            cout << "Converged after " << counter << " iterations\n";
            break;
        } else {
            previous = residual;
            counter++;
        }
        iterations--;
    }
    return x;
}