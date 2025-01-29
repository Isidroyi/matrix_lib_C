
#include "utils.h"

#include "s21_matrix.h"

void s21_create_minor_matrix(int I, int J, matrix_t *A, matrix_t *minor) {
  int i_minor = 0;
  int j_minor = 0;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (i != I && j != J) {
        minor->matrix[i_minor][j_minor] = A->matrix[i][j];
        j_minor++;
        if (j_minor == A->columns - 1) {
          j_minor = 0;
          i_minor++;
        }
      }
    }
  }
}

int s21_is_valid(matrix_t const *matrix) {
  int error = SUCCESS;
  if (matrix == NULL || matrix->matrix == NULL || matrix->rows <= 0 ||
      matrix->columns <= 0) {
    error = FAILURE;
  }
  return error;
}

void initialization_matrix(double n, matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        A->matrix[i][j] = n;
      }
    }
  }
}