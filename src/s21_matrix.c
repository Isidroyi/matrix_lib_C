#include "s21_matrix.h"

#include "utils.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = OK;
  if (result == NULL || rows <= 0 || columns <= 0) {
    error = INCORRECT_MATRIX;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix != NULL) {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
        if (result->matrix[i] == NULL) {
          error = INCORRECT_MATRIX;
          break;
        }
      }
    } else {
      error = INCORRECT_MATRIX;
    }
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
    A->columns = 0;
    A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int result_eq = SUCCESS;
  if (s21_is_valid(A) && s21_is_valid(B) && A->rows == B->rows &&
      A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1.0E-7) {
          result_eq = FAILURE;
          break;
        }
      }
    }
  } else {
    result_eq = FAILURE;
  }
  return result_eq;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int result_sum = OK;
  if (s21_is_valid(A) && s21_is_valid(B) && result != NULL) {
    if (A->rows == B->rows && A->columns == B->columns) {
      if (s21_create_matrix(A->rows, A->columns, result) == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
          }
        }
      } else {
        result_sum = INCORRECT_MATRIX;
      }
    } else {
      result_sum = CALCULATION_ERROR;
    }
  } else {
    result_sum = INCORRECT_MATRIX;
  }
  return result_sum;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int result_sub = OK;
  if (s21_is_valid(A) && s21_is_valid(B) && result != NULL) {
    if (A->rows == B->rows && A->columns == B->columns) {
      if (s21_create_matrix(A->rows, A->columns, result) == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
          }
        }
      } else {
        result_sub = INCORRECT_MATRIX;
      }
    } else {
      result_sub = CALCULATION_ERROR;
    }
  } else {
    result_sub = INCORRECT_MATRIX;
  }
  return result_sub;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int result_mult_num = OK;
  if (s21_is_valid(A) && result != NULL) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else {
      result_mult_num = INCORRECT_MATRIX;
    }
  } else {
    result_mult_num = INCORRECT_MATRIX;
  }
  return result_mult_num;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int result_mult = OK;
  if (s21_is_valid(A) && s21_is_valid(B) && result != NULL) {
    if (A->columns == B->rows) {
      if (s21_create_matrix(A->rows, B->columns, result) == OK) {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < B->columns; j++) {
            for (int k = 0; k < A->columns; k++) {
              result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
            }
          }
        }
      } else {
        result_mult = INCORRECT_MATRIX;
      }
    } else {
      result_mult = CALCULATION_ERROR;
    }
  } else {
    result_mult = INCORRECT_MATRIX;
  }
  return result_mult;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int result_transp = OK;
  if (s21_is_valid(A) && result != NULL) {
    if (s21_create_matrix(A->columns, A->rows, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[j][i] = A->matrix[i][j];
        }
      }
    } else {
      result_transp = INCORRECT_MATRIX;
    }
  } else {
    result_transp = INCORRECT_MATRIX;
  }
  return result_transp;
}

int s21_determinant(matrix_t *A, double *result) {
  int result_det = OK;
  if (s21_is_valid(A)) {
    if (A->rows == A->columns) {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[1][0] * A->matrix[0][1];
      } else {
        int sign = 1;
        matrix_t minor = {0};
        double minor_det;
        *result = 0;
        for (int i = 0; i < A->columns; i++) {
          s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
          s21_create_minor_matrix(0, i, A, &minor);
          s21_determinant(&minor, &minor_det);
          *result += sign * A->matrix[0][i] * minor_det;
          sign *= -1;
          s21_remove_matrix(&minor);
        }
      }
    } else {
      result_det = CALCULATION_ERROR;
    }
  } else {
    result_det = INCORRECT_MATRIX;
  }
  return result_det;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int result_compl = OK;
  if (s21_is_valid(A) && result != NULL) {
    if (A->rows == A->columns) {
      if (s21_create_matrix(A->rows, A->columns, result) == OK) {
        matrix_t minor = {0};
        double det;
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
            s21_create_minor_matrix(i, j, A, &minor);
            s21_determinant(&minor, &det);
            result->matrix[i][j] = pow(-1, (i + j)) * det;
            s21_remove_matrix(&minor);
          }
        }
      } else {
        result_compl = INCORRECT_MATRIX;
      }
    } else {
      result_compl = CALCULATION_ERROR;
    }
  } else {
    result_compl = INCORRECT_MATRIX;
  }
  return result_compl;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int result_inv = OK;
  if (s21_is_valid(A) && result != NULL) {
    if (A->rows == A->columns) {
      double det;
      s21_determinant(A, &det);
      if (det != 0) {
        if (s21_create_matrix(A->rows, A->columns, result) == OK) {
          matrix_t tmp_calc = {0};
          s21_calc_complements(A, &tmp_calc);
          matrix_t transp_tmp_calc = {0};
          s21_transpose(&tmp_calc, &transp_tmp_calc);
          s21_mult_number(&transp_tmp_calc, 1 / det, result);
          s21_remove_matrix(&transp_tmp_calc);
          s21_remove_matrix(&tmp_calc);
        } else {
          result_inv = INCORRECT_MATRIX;
        }
      } else {
        result_inv = CALCULATION_ERROR;
      }
    } else {
      result_inv = CALCULATION_ERROR;
    }
  } else {
    result_inv = INCORRECT_MATRIX;
  }
  return result_inv;
}
