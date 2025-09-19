#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "bateria.h"
#include "emf.h"
#include "matrix.h"

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

int main() {
    return 0;
}
