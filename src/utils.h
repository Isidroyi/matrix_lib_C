#ifndef UTILS_H
#define UTILS_H

#include "s21_matrix.h"

typedef struct matrix_struct matrix_t;
void s21_create_minor_matrix(int I, int J, matrix_t *A, matrix_t *minor);
int s21_is_valid(matrix_t const *matrix);
void initialization_matrix(double n, matrix_t *A);

#endif