//
// Created by Adam on 1/13/2017.
//

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
using namespace std;

class Matrix {
public:
    Matrix();
    Matrix(int r, int c);
    Matrix(const Matrix& rhs);
    ~Matrix();

    void setValue(int r, int c, double value);

    int getRows() const;
    int getCols() const;
    void setRows(int r);
    void setCols(int c);

    /**
     * Turns the matrix into an identity matrix
     */
    void identityMatrix();

    /**
     * Fills the matrix with random values between -max and max
     * @param max
     */
    void randomMatrix(double max);

    void randomIntMatrix(int max);
    void randomIntVector(int max);

    /**
     * Use Guassian elimination to reduce a matrix to row echelon form.
     */
    void rowEchelonForm();

    /**
     * Create an augmented matrix, adding Matrix rhs onto the right of an existing matrix.
     * @param rhs
     */
    void augmentMatrix(const Matrix& rhs);

    /**
     * Backsolve an augmented matrix
     * @return The solution matrix
     */
    Matrix backSolve();

    /**
     * Backsolve a non-augmented matrix given a rhs matrix b
     * @param b rhs matrix
     * @return The solution matrix
     */
    Matrix backSolve(const Matrix& b);
    double calculateResidual(Matrix& calculatedSolution, Matrix& rhs);

    double magnitude();

    // LU Decomposition - Doolittle method
    void decompose(Matrix& L, Matrix& U);
    Matrix forwardSolve(const Matrix& b);

    // Jacobi iteration
    Matrix jacobi(Matrix& b, int numIterations);

    // Gauss-Seidel
    Matrix gaussSeidel(Matrix &b, int iterations);

    // Successive Over-Relaxation
    Matrix sor(Matrix &b, int iterations, double omega);

    // Power iteration
    double powerIteration(Matrix& b, int iterations);

    friend ostream& operator<<(ostream& os, const Matrix& matrix);
    double operator()(int r, int c) const;
    void operator=(const Matrix& rhs);
    Matrix operator-(const Matrix& rhs);
    Matrix operator+(const Matrix& rhs);
    Matrix operator*(const Matrix& rhs);
    Matrix operator*(double scalar);
    Matrix operator/(double scalar);

private:
    int rows, cols;
    double** M;
};


#endif //MATRIX_H