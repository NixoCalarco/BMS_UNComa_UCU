#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

/*
 * las matrices van a ser de 2*2, 2*1 o 1*2
 * desperdiciamos 16 bytes cada tanto (2 doubles), pero los tamaños son fijos
 * durante la compilación y no tenemos que jorobar con memoria dinámica
 *
 * además, definir vectores fila / columna es trivial, porque las celdas no
 * especificadas quedan en 0:
 *  --> Matrix v1 = {
 *          1, 2,
 *          {
 *              {x, y}
 *          }
 *      };
 *  --> Matrix v2 = {
 *          2, 1,
 *          {
 *              {x},
 *              {y}
 *          }
 *      };
 *
 */
typedef struct {
    size_t fil, col;
    double mat[2][2];
} Matrix;

Matrix transpose(Matrix A);

Matrix matrix_mult(Matrix A, Matrix B);

Matrix scalar_sum(double s, Matrix A);

Matrix scalar_mult(double s, Matrix A);

Matrix scalar_div(double s, Matrix A);

#endif
