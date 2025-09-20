#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "bateria.h"
#include "emf.h"
#include "matrix.h"

typedef struct {
    size_t len, capacity;
    double *data;
} DArray;

DArray DArray_new() {
    DArray v = {0, 10, malloc(sizeof(double) * 10)};

    if (!v.data) {
        puts("oops");
        exit(1);
    }

    return v;
}

void DArray_push(DArray *this, double x) {
    if (this->len == this->capacity) {
        this->capacity *= 2;
        this->data = malloc(sizeof(double) * this->capacity);
        if (!this->data) {
            puts("oops");
            exit(1);
        }
    }

    this->data[this->len] = x;
    this->len += 1;
}

void dump_to_file() {
    FILE *f = fopen("data.txt", "w");
    for (;;) {
        fprintf(f, "", "");
    }
    fclose(f);
}

void test() {
    Matrix A = {2, 2, {{1, 2}, {3, 4}}};

    Matrix T = transpose(A);

    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < A.col; ++j)
            printf("%f ", A.mat[i][j]);
        puts("");
    }

    puts("");

    for (size_t i = 0; i < T.fil; ++i) {
        for (size_t j = 0; j < T.col; ++j)
            printf("%f ", T.mat[i][j]);
        puts("");
    }

    puts("");

    Matrix a = {2, 1, {{1, 2}}};

    Matrix t = transpose(a);

    for (size_t i = 0; i < a.fil; ++i) {
        for (size_t j = 0; j < a.col; ++j)
            printf("%f ", a.mat[i][j]);
        puts("");
    }

    puts("");

    for (size_t i = 0; i < t.fil; ++i) {
        for (size_t j = 0; j < t.col; ++j)
            printf("%f ", t.mat[i][j]);
        puts("");
    }

    puts("");

    Matrix B = {2, 2, {{2, 4}, {6, 8}}};
    Matrix C = matrix_mult(A, B);

    for (size_t i = 0; i < C.fil; ++i) {
        for (size_t j = 0; j < C.col; ++j)
            printf("%f ", C.mat[i][j]);
        puts("");
    }

    Matrix D = {2, 1, {{5}, {8}}};
    Matrix E = matrix_mult(A, D);

    for (size_t i = 0; i < E.fil; ++i) {
        for (size_t j = 0; j < E.col; ++j)
            printf("%f ", E.mat[i][j]);
        puts("");
    }
}

void print_matrix(Matrix A) {
    puts("");
    for (size_t i = 0; i < A.fil; ++i) {
        for (size_t j = 0; j < A.col; ++j)
            printf("%f ", A.mat[i][j]);
        puts("");
    }
}

int main() {
    constexpr int cant_muestras = BAT_N;
    double muestreo_deseado = 1. / 3600.;

    double bat_ci = .1;
    double capacidad = 7.7; // A/h

    double tension_min = 3.05;

    double a = 9.25e-02;
    double b = 3.8125e-02;
    double r = 5.125e-02;
    Matrix A = {
        2,
        2,
        {{1, 0}, {1 - exp(-muestreo_deseado / b), exp(-muestreo_deseado / b)}},
    };
    Matrix B = {
        1,
        2,
        {{-1, (b - a) * (1 - exp(-muestreo_deseado / b)) - 1}},
    };
    B = scalar_mult(muestreo_deseado / capacidad, B);

    for (int i = 0; i < cant_muestras; ++i) {
        corriente[i] *= -1;
        voltaje[i] /= 1'000;
    }

    // observador
    double R2 = 1.;
    Matrix Pb = {2, 2, {{1, 0}, {0, 1}}};
    Matrix Xb = {1, 2, {{bat_ci, bat_ci}}};

    double jk = 1.;
    double Vhat[cant_muestras] = {};
    double Vh[cant_muestras] = {};
    Vhat[0] = inverse_emf(bat_ci);

    for (int EKF = 0; EKF < 2; ++EKF) {
        double xxmin[cant_muestras] = {};
        double TT[cant_muestras] = {};
        double TE[cant_muestras] = {};
        double Se[cant_muestras] = {};
        double XXe[cant_muestras] = {};
        double Ps[cant_muestras] = {};
        double Px[cant_muestras] = {};
        double J[2] = {};
        DArray ne = DArray_new();
        DArray e = DArray_new();

        for (int i = 0; i < cant_muestras; ++i) {
            double Xm = emf(voltaje[i] + corriente[i] * r);
            double df = derivative_emf(Xb.mat[0][1]);

            if (!EKF) {
                Matrix R1 = {2, 2, {{1, 0}, {0, 10}}};
                Matrix C = {1, 2, {{0, 1}}};

                Matrix K = scalar_mult(
                    (1. / (R2 + (matrix_mult(matrix_mult(C, Pb), transpose(C)))
                                    .mat[0][0])),
                    matrix_mult(Pb, transpose(C)));

                Matrix Xe = matrix_sum(Xb, scalar_mult(Xm - Xb.mat[0][1], K));

                Matrix P = inverse(matrix_sum(
                    inverse(Pb), scalar_div(R2, matrix_mult(transpose(C), C))));

                Xb = matrix_sum(matrix_mult(A, Xe),
                                scalar_mult(corriente[i], B));
                Pb = matrix_sum(matrix_mult(matrix_mult(A, P), transpose(A)),
                                R1);

                Vh[i] = inverse_emf(Xb.mat[0][1]) - corriente[i] * r;
            } else {
                Matrix R1 = {2, 2, {{.1, 0}, {0, 10}}};
                Matrix C = {1, 2, {{0, df}}};
            }
        }
    }

    return 0;
}
