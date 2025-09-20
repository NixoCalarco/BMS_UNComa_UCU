#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

Matrix transpose(Matrix A) {
    Matrix T = {A.col, A.fil, {}};

    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < A.col; ++j) {
            T.mat[j][i] = A.mat[i][j];
        }
    }
    return T;
}

Matrix matrix_mult(Matrix A, Matrix B) {
    if (A.col != B.fil) {
        puts("oops");
        exit(1);
    }

    Matrix C = {A.fil, B.col, {}};

    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < B.col; ++j) {
            C.mat[i][j] = 0;
            for (size_t k = 0; k < B.fil; ++k) {
                C.mat[i][j] += A.mat[i][k] * B.mat[k][j];
            }
        }
    }

    return C;
}

Matrix matrix_sum(Matrix A, Matrix B) {
    return (Matrix){
        2,
        2,
        {{A.mat[0][0] + B.mat[0][0], A.mat[0][1] + B.mat[0][1]},
         {A.mat[1][0] + B.mat[1][0], A.mat[1][1] + B.mat[1][1]}},
    };
}

Matrix scalar_sum(double s, Matrix A) {
    Matrix B = {A.fil, A.col, {}};
    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < A.col; ++j) {
            B.mat[i][j] = s + A.mat[i][j];
        }
    }

    return B;
}

Matrix scalar_mult(double s, Matrix A) {
    Matrix B = {A.fil, A.col, {}};
    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < A.col; ++j) {
            B.mat[i][j] = s * A.mat[i][j];
        }
    }

    return B;
}

Matrix scalar_div(double s, Matrix A) {
    Matrix B = {A.fil, A.col, {}};
    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < A.col; ++j) {
            B.mat[i][j] = A.mat[i][j] / s;
        }
    }

    return B;
}

Matrix inverse(Matrix A) {
    if (A.fil != 2 || A.col != 2) {
        puts("oops");
        exit(1);
    }

    // por suerte invertir una 2x2 es una pavada
    double determinante =
        (A.mat[0][0] * A.mat[1][1]) - (A.mat[0][1] * A.mat[1][0]);

    if (determinante == 0) {
        puts("oops");
        exit(1);
    }

    Matrix B = {
        2, 2, {{A.mat[1][1], -A.mat[0][1]}, {-A.mat[1][0], A.mat[0][0]}}};

    return scalar_div(determinante, B);
}
