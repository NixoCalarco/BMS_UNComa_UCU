#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bateria.h"
#include "emf.h"
#include "lamb.h"
#include "matrix.h"

/**
 * algunas utilidades
 */

constexpr int cant_muestras = BAT_N

    typedef double
        d_buff[cant_muestras]; // tipo array para facilitar el manejo de memoria

static void *alloc(size_t size) {
    void *p = malloc(size);
    if (!p) {
        puts("falló malloc!");
        exit(1);
    }
    memset(p, 0, size);
    return p;
}

static void clear(d_buff p) { memset(p, 0, sizeof(d_buff)); }

static void cleanup_free(void *p) { free(*(void **)p); }

#define defer_cleanup __attribute__((cleanup(cleanup_free)))

// tipo array con un índice para poder insertar valores la final
// typedef struct {
//     size_t len;
//     d_buff data;
// } d_vec;
//
// static void push_value(double val, d_vec *v) {
//     if (v->len >= sizeof(d_buff) / sizeof(double)) {
//         puts("te quedaste sin espacio");
//         exit(1);
//     }
//
//     v->data[v->len] = val;
//
//     ++(v->len);
// }
//
// static double norm(size_t n, double v[static n]) {
//     double x = 0.;
//     for (size_t i = 0; i < n; ++i)
//         x += v[i] * v[i];
//
//     return sqrt(x);
// }

// static int comp(const void *a, const void *b) {
//     double *_a = (double *)a;
//     double *_b = (double *)b;
//
//     return (*_a < *_b) ? -1 : (*_a > *_b) ? 1 : 0;
// }
//
// static double median(size_t n, double v[static n]) {
//     qsort((void *)v, n, sizeof(double), comp);
//
//     double m = 0;
//     if (n % 2 == 1) {
//         m = v[n / 2];
//     } else {
//         m = (v[n / 2 - 1] + v[n / 2]) / 2.;
//     }
//
//     return m;
// }
//
// static void print_matrix(Matrix A) {
//     puts("");
//     for (size_t i = 0; i < A.fil; ++i) {
//         for (size_t j = 0; j < A.col; ++j)
//             printf("%f ", A.mat[i][j]);
//         puts("");
//     }
// }

static void pyprint(size_t n, double v[static n], char *name) {
    printf("%s = [", name);
    for (size_t i = 0; i < n; ++i)
        printf("%f,", v[i]);
    puts("]");
}

struct update_params {
    const double R2;
    const Matrix C;
    const Matrix Xb;
    const Matrix Pb;
    const double diff; // Xm - Xb en normal, V - Vhat en extendido
};

static void kalman_update(const struct update_params *p, Matrix *Xe,
                          Matrix *P) {
    Matrix K = scalar_mult(
        (1. /
         (p->R2 +
          (matrix_mult(matrix_mult(p->C, p->Pb), transpose(p->C))).mat[0][0])),
        matrix_mult(p->Pb, transpose(p->C)));
    *Xe = matrix_sum(p->Xb, scalar_mult(p->diff, K));

    *P = inverse(matrix_sum(
        inverse(p->Pb), scalar_div(p->R2, matrix_mult(transpose(p->C), p->C))));
}

struct predict_params {
    const double I;
    const Matrix A;
    const Matrix B;
    const Matrix R1;
    const Matrix P;
    const Matrix Xe;
};

static void kalman_predict(const struct predict_params *p, Matrix *Xb,
                           Matrix *Pb) {

    *Xb = matrix_sum(matrix_mult(p->A, p->Xe), scalar_mult(p->I, p->B));
    *Pb = matrix_sum(matrix_mult(matrix_mult(p->A, p->P), transpose(p->A)),
                     p->R1);
}

int main() {
    constexpr double muestreo_deseado = 1. / 3600.;

    constexpr double bat_ci = .1;
    constexpr double capacidad = 7.7; // A/h

    constexpr double tension_min = 3.05;

    constexpr double a = 9.25e-02;
    constexpr double b = 3.8125e-02;
    constexpr double r = 5.125e-02;
    const Matrix A = {
        2,
        2,
        {{1, 0}, {1 - exp(-muestreo_deseado / b), exp(-muestreo_deseado / b)}},
    };
    const Matrix B = scalar_mult(
        muestreo_deseado / capacidad,
        (Matrix){
            2,
            1,
            {{-1}, {(b - a) * (1 - exp(-muestreo_deseado / b)) - 1}},
        });

    for (int i = 0; i < cant_muestras; ++i) {
        corriente[i] *= -1;
        voltaje[i] /= 1'000;
    }

    // observador
    constexpr double R2 = 1.;
    // Matrix Pb = {2, 2, {{1, 0}, {0, 1}}};
    // Matrix Xb = {1, 2, {{bat_ci, bat_ci}}};

    // double jk = 1.;
    defer_cleanup double *Vhat = alloc(sizeof(d_buff));
    defer_cleanup double *Vh = alloc(sizeof(d_buff));

    // double J[2] = {};

    defer_cleanup double *Se = alloc(sizeof(d_buff));
    defer_cleanup double *XXe = alloc(sizeof(d_buff));

    defer_cleanup double *TT = alloc(sizeof(d_buff));
    defer_cleanup double *TE = alloc(sizeof(d_buff));

    for (int EKF = 0; EKF < 2; ++EKF) {
        Matrix Pb = {2, 2, {{1, 0}, {0, 1}}};
        Matrix Xb = {2, 1, {{bat_ci}, {bat_ci}}};
        // double xxmin[cant_muestras] = {};

        clear(Se);
        clear(XXe);
        clear(TT);
        clear(TE);

        // double Ps[cant_muestras] = {};
        // double Px[cant_muestras] = {};
        // DArray ne = DArray_new();
        // DArray e = DArray_new();

        // defer_cleanup d_vec *ne = alloc(sizeof(d_vec));
        // defer_cleanup d_vec *e = alloc(sizeof(d_vec));

        for (int i = 0; i < cant_muestras; ++i) {
            const double Xm = emf(voltaje[i] + corriente[i] * r);
            const double df = derivative_emf(Xb.mat[0][1]);

            if (!EKF) {
                Vh[i] = inverse_emf(Xb.mat[1][0]) - corriente[i] * r;

                constexpr Matrix R1 = {2, 2, {{1, 0}, {0, 10}}};
                constexpr Matrix C = {1, 2, {{0, 1}}};

                // corregimos estado actual
                Matrix Xe, P;
                struct update_params up = {R2, C, Xb, Pb, Xm - Xb.mat[1][0]};
                kalman_update(&up, &Xe, &P);

                // predecimos estado futuro
                struct predict_params pp = {corriente[i], A, B, R1, P, Xe};
                kalman_predict(&pp, &Xb, &Pb);
            } else {
                Vhat[i] = inverse_emf(Xb.mat[1][0]) - corriente[i] * r;

                constexpr Matrix R1 = {2, 2, {{0.1, 0}, {0, 10}}};
                const Matrix C = {1, 2, {{0, df}}};

                // corregimos estado actual
                Matrix Xe, P;
                struct update_params up = {R2, C, Xb, Pb, voltaje[i] - Vhat[i]};
                kalman_update(&up, &Xe, &P);

                // Ps[i] = sqrt(P.mat[0][0]);
                // Px[i] = sqrt(P.mat[1][1]);

                // predecimos estado futuro
                struct predict_params pp = {corriente[i], A, B, R1, P, Xe};
                kalman_predict(&pp, &Xb, &Pb);
            }

            Xb.mat[0][0] = Xb.mat[0][0] < 0   ? 0
                           : Xb.mat[0][0] > 1 ? 1
                                              : Xb.mat[0][0];
            Xb.mat[1][0] = Xb.mat[1][0] < 0   ? 0
                           : Xb.mat[1][0] > 1 ? 1
                                              : Xb.mat[1][0];

            Se[i] = Xb.mat[0][0];
            XXe[i] = Xb.mat[1][0];

            int Te = 0, j = 1;
            while (corriente[i] > 0 && i + j < cant_muestras &&
                   voltaje[i + j] > tension_min) {
                ++Te;
                ++j;
            }

            if (corriente[i] > 0 && voltaje[i] > tension_min) {
                TE[i] = muestreo_deseado * Te;
                double xmin = emf(tension_min + corriente[i] * r);
                // xxmin[i] = xmin;

                double t1 = (xmin - Se[i]) * capacidad / corriente[i] - (b - a);
                double t2 =
                    (Se[i] - XXe[i]) * capacidad / corriente[i] + (b - a);

                double arg = -(t2 / b) * exp(t1 / b);

                int m2 = 0;
                for (int i = 0; i < LAMB_N; ++i) {
                    double x1 = fabs(lambert[m2].x - arg);
                    double x2 = fabs(lambert[i].x - arg);
                    if (x2 < x1)
                        m2 = i;
                }

                double w = lambert[m2].y;

                TT[i] = w * b - t1;
                // push_value(TE[i] - TT[i], e);
                // ++jk;
            }

            // if (i > 0 && corriente[i] < 0 && corriente[i - 1] > 0) {
            //     push_value(100 * norm(e->len, e->data) / sqrt(e->len), ne);
            //     clear(e->data);
            //     e->len = 0;
            //     // jk = 1;
            // }
        }

        // DArray ne_copy = DArray_new();
        // for (size_t i = 0; i < ne.len; ++i) {
        //     double elem = ne.data[i];
        //     if (elem < 60)
        //         DArray_push(&ne_copy, elem);
        // }
        // // J[EKF] = median(ne_copy.len, ne_copy.data);
        //
        // free(ne.data);
        // free(ne_copy.data);
        // free(e.data);

        if (EKF) {
            pyprint(cant_muestras, TE, "TE1");
            pyprint(cant_muestras, TT, "TT1");
            pyprint(cant_muestras, Se, "Se1");
            pyprint(cant_muestras, XXe, "XXe1");
        } else {
            pyprint(cant_muestras, TE, "TE0");
            pyprint(cant_muestras, TT, "TT0");
            pyprint(cant_muestras, Se, "Se0");
            pyprint(cant_muestras, XXe, "XXe0");
        }
    }

    // double minJ = J[0] < J[1] ? J[0] : J[1];
    // int idx = J[0] < J[1] ? 0 : 1;

    pyprint(cant_muestras, corriente, "Corriente");

    pyprint(cant_muestras, voltaje, "tension");
    pyprint(cant_muestras, Vhat, "Vhat");
    pyprint(cant_muestras, Vh, "Vh");

    return 0;
}
