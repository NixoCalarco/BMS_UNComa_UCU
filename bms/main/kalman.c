#include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>

// #include "bateria.h"
#include "bat.h"
// #include "bat_stub.h"
#include "emf.h"
#include "kalman.h"
#include "lamb.h"
#include "matrix.h"

/**
 * algunas utilidades
 */

constexpr int cant_muestras = BAT_N;

// typedef double
//     d_buff[cant_muestras]; // tipo array para facilitar el manejo de memoria

// static void *alloc(size_t size) {
//     void *p = malloc(size);
//     if (!p) {
//         puts("falló malloc!");
//         exit(1);
//     }
//     memset(p, 0, size);
//     return p;
// }

// static void clear(d_buff p) { memset(p, 0, sizeof(d_buff)); }

// static void cleanup_free(void *p) { free(*(void **)p); }

// #define defer_cleanup __attribute__((cleanup(cleanup_free)))

// static void pyprint(size_t n, double v[static n], char *name) {
//     printf("%s = [", name);
//     for (size_t i = 0; i < n; ++i)
//         printf("%f,", v[i]);
//     puts("]");
// }

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

void kalman(void (*send_values)(KalmanVals vals)) {
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

    // for (int i = 0; i < cant_muestras; ++i) {
    //     corriente[i] *= -1;
    //     voltaje[i] /= 1'000;
    // }

    // observador
    constexpr double R2 = 1.;

    // defer_cleanup double *Vhat = alloc(sizeof(d_buff));
    // defer_cleanup double *Vh = alloc(sizeof(d_buff));
    //
    // defer_cleanup double *Se = alloc(sizeof(d_buff));
    // defer_cleanup double *XXe = alloc(sizeof(d_buff));
    //
    // defer_cleanup double *TT = alloc(sizeof(d_buff));
    // defer_cleanup double *TE = alloc(sizeof(d_buff));

    // for (int EKF = 0; EKF < 2; ++EKF) {
    for (int EKF = 1; EKF >= 0; --EKF) {
        Matrix Pb = {2, 2, {{1, 0}, {0, 1}}};
        Matrix Xb = {2, 1, {{bat_ci}, {bat_ci}}};

        // clear(Se);
        // clear(XXe);
        // clear(TT);
        // clear(TE);

        for (int i = 0; i < cant_muestras; ++i) {
            const double Xm = emf(voltaje[i] + corriente[i] * r);
            const double df = derivative_emf(Xb.mat[0][1]);

            double V = 0, Se = 0, XXe = 0, TT = 0, TE = 0;

            if (!EKF) {
                V = inverse_emf(Xb.mat[1][0]) - corriente[i] * r;

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
                V = inverse_emf(Xb.mat[1][0]) - corriente[i] * r;

                constexpr Matrix R1 = {2, 2, {{0.1, 0}, {0, 10}}};
                const Matrix C = {1, 2, {{0, df}}};

                // corregimos estado actual
                Matrix Xe, P;
                struct update_params up = {R2, C, Xb, Pb, voltaje[i] - V};
                kalman_update(&up, &Xe, &P);

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

            Se = Xb.mat[0][0];
            XXe = Xb.mat[1][0];

            int Te = 0, j = 1;
            while (corriente[i] > 0 && i + j < cant_muestras &&
                   voltaje[i + j] > tension_min) {
                ++Te;
                ++j;
            }

            if (corriente[i] > 0 && voltaje[i] > tension_min) {
                TE = muestreo_deseado * Te;
                double xmin = emf(tension_min + corriente[i] * r);

                double t1 = (xmin - Se) * capacidad / corriente[i] - (b - a);
                double t2 = (Se - XXe) * capacidad / corriente[i] + (b - a);

                double arg = -(t2 / b) * exp(t1 / b);

                int m2 = 0;
                for (int i = 0; i < LAMB_N; ++i) {
                    double x1 = fabs(lambert[m2].x - arg);
                    double x2 = fabs(lambert[i].x - arg);
                    if (x2 < x1)
                        m2 = i;
                }

                double w = lambert[m2].y;

                TT = w * b - t1;
            }

            // pateamos todo de vuelta
            send_values((KalmanVals){
                .algorithm = EKF ? EXTENDED : NORMAL,
                .V = V,
                .Se = Se,
                .XXe = XXe,
                .TT = TT,
                .TE = TE,
            });
        }

        // if (EKF) {
        //     pyprint(cant_muestras, TE, "TE1");
        //     pyprint(cant_muestras, TT, "TT1");
        //     pyprint(cant_muestras, Se, "Se1");
        //     pyprint(cant_muestras, XXe, "XXe1");
        // } else {
        //     pyprint(cant_muestras, TE, "TE0");
        //     pyprint(cant_muestras, TT, "TT0");
        //     pyprint(cant_muestras, Se, "Se0");
        //     pyprint(cant_muestras, XXe, "XXe0");
        // }
    }

    // pyprint(cant_muestras, corriente, "Corriente");
    //
    // pyprint(cant_muestras, voltaje, "tension");
    // pyprint(cant_muestras, Vhat, "Vhat");
    // pyprint(cant_muestras, Vh, "Vh");
}
