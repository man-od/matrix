#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = OK;
  int i = 0;
  if (rows < 1 || columns < 1) {
    status = ERROR_MATRIX;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    for (i = 0; i < rows; i++) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
    }
  }
  return status;
}

void s21_remove_matrix(matrix_t *A) {
  int i = 0;
  if (A->matrix) {
    for (i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = SUCCESS;
  int i = 0;
  int j = 0;
  if ((A->matrix != NULL || A->rows > 0 || A->columns > 0) ||
      (B->matrix != NULL || B->rows > 0 || B->columns > 0)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      for (i = 0; i < A->rows; i++) {
        for (j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
            status = FAILURE;
          }
        }
      }
    } else {
      status = FAILURE;
    }
  } else {
    status = FAILURE;
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  result->matrix = NULL;
  int i = 0;
  int j = 0;
  if ((A->matrix != NULL || A->rows > 0 || A->columns > 0) ||
      (B->matrix != NULL || B->rows > 0 || B->columns > 0)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (i = 0; i < result->rows; i++) {
        for (j = 0; j < result->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      status = ERROR_CALC;
    }
  } else {
    status = ERROR_MATRIX;
  }
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  result->matrix = NULL;
  int i = 0;
  int j = 0;
  if ((A->matrix != NULL || A->rows > 0 || A->columns > 0) ||
      (B->matrix != NULL || B->rows > 0 || B->columns > 0)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      for (i = 0; i < result->rows; i++) {
        for (j = 0; j < result->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      status = ERROR_CALC;
    }
  } else {
    status = ERROR_MATRIX;
  }
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = OK;
  result->matrix = NULL;
  int i = 0;
  int j = 0;
  if (A->matrix != NULL || A->rows > 0 || A->columns > 0) {
    s21_create_matrix(A->rows, A->columns, result);
    for (i = 0; i < result->rows; i++) {
      for (j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  } else {
    status = ERROR_MATRIX;
  }
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  result->matrix = NULL;
  int i = 0;
  int j = 0;
  int k = 0;
  if ((A->matrix != NULL || A->rows > 0 || A->columns > 0) ||
      (B->matrix != NULL || B->rows > 0 || B->columns > 0)) {
    if (A->columns == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (i = 0; i < result->rows; i++) {
        for (j = 0; j < result->columns; j++) {
          for (k = 0; k < A->columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      status = ERROR_CALC;
    }
  } else {
    status = ERROR_MATRIX;
  }
  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = OK;
  result->matrix = NULL;
  int i = 0;
  int j = 0;
  if (A->matrix != NULL || (A->rows > 0 && A->columns > 0)) {
    s21_create_matrix(A->columns, A->rows, result);
    for (i = 0; i < A->columns; i++) {
      for (j = 0; j < A->rows; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  } else {
    status = ERROR_MATRIX;
  }
  return status;
}

void s21_getMinorFunc(int ii, int jj, matrix_t *A, matrix_t *minormatrix_t) {
  s21_create_matrix(A->rows - 1, A->columns - 1, minormatrix_t);
  int i_m = 0;
  for (int i = 0; i < A->rows; i++) {
    int j_m = 0;
    if (i == ii) continue;
    for (int j = 0; j < A->columns; j++) {
      if (j != jj) {
        minormatrix_t->matrix[i_m][j_m] = A->matrix[i][j];
        j_m++;
      }
    }
    i_m++;
  }
}

double s21_getdet(int n, matrix_t *A) {
  double res = 0;
  if (n == 1) {
    res = A->matrix[0][0];
  } else if (n == 2) {
    res = A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
  } else {
    for (int c = 0; c < n; c++) {
      matrix_t submat;
      s21_create_matrix(n - 1, n - 1, &submat);
      int subi = 0;
      for (int i = 1; i < n; i++) {
        int subj = 0;
        for (int j = 0; j < n; j++) {
          if (j != c) {
            submat.matrix[subi][subj] = A->matrix[i][j];
            subj++;
          }
        }
        subi++;
      }
      double sign = (c % 2 == 0) ? 1.0 : -1.0;
      res += sign * A->matrix[0][c] * s21_getdet(n - 1, &submat);
      s21_remove_matrix(&submat);
    }
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = OK;
  if (A == NULL || A->rows < 1 || A->columns < 1) {
    status = ERROR_MATRIX;
  } else if (A->columns != A->rows) {
    status = ERROR_MATRIX;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else {
      *result = s21_getdet(A->rows, A);
    }
  }
  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;
  int i = 0, j = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1) {
    status = ERROR_MATRIX;
  } else {
    s21_create_matrix(A->columns, A->rows, result);
    for (i = 0; i < result->rows; i++) {
      for (j = 0; j < result->columns; j++) {
        matrix_t minor_mat;
        s21_getMinorFunc(i, j, A, &minor_mat);
        double res;
        s21_determinant(&minor_mat, &res);
        result->matrix[i][j] = pow((-1), i + j) * res;
        s21_remove_matrix(&minor_mat);
      }
    }
  }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  result->matrix = NULL;
  int status = OK;
  double res;
  if (A == NULL || (A->rows < 1 || A->columns < 1)) {
    status = ERROR_MATRIX;
  } else if (A->columns != A->rows) {
    status = ERROR_CALC;
  } else if (!s21_determinant(A, &res) && !res) {
    status = ERROR_CALC;
  } else {
    matrix_t inverse, transpose;
    s21_calc_complements(A, &inverse);
    s21_transpose(&inverse, &transpose);
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = 1 / res * transpose.matrix[i][j];
      }
    }
    s21_remove_matrix(&inverse);
    s21_remove_matrix(&transpose);
  }
  return status;
}
