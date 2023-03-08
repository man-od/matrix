#include <check.h>

#include "s21_matrix.h"

START_TEST(test_s21_create_matrix_3) {
  matrix_t A;
  ck_assert_int_eq(s21_create_matrix(3, 3, &A), 0);
  ck_assert_int_eq(A.matrix[0][0], 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_create_matrix_4) {
  matrix_t A;
  ck_assert_int_eq(s21_create_matrix(5, 2, &A), 0);
  ck_assert_int_eq(A.matrix[4][0], 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_create_matrix_5) {
  matrix_t A;
  ck_assert_int_eq(s21_create_matrix(10, 10, &A), 0);
  ck_assert_int_eq(A.matrix[9][9], 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_eq_matrix_2) {
  matrix_t A;
  matrix_t B;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(1, 1, &B);
  int fail3 = s21_eq_matrix(&A, &B);  // 1 - sucess

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_eq_matrix_3) {
  matrix_t A;
  matrix_t B;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int fail3 = s21_eq_matrix(&A, &B);  // 1 - sucess

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_eq_matrix_4) {
  matrix_t A;
  matrix_t B;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][1] = 1;

  int fail3 = s21_eq_matrix(&A, &B);  // 1 - sucess

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_eq_matrix_5) {
  matrix_t A;
  matrix_t B;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 3;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_eq_matrix(&A, &B);  // 1 - sucess

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_eq_matrix_6) {
  matrix_t A;
  matrix_t B;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 5;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_eq_matrix_7) {
  matrix_t A;
  matrix_t B;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = -5;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_eq_matrix(&A, &B);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_sum_matrix_2) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(1, 1, &B);
  int fail3 = s21_sum_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sum_matrix_3) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int fail3 = s21_sum_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sum_matrix_4) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][1] = 1;
  int fail3 = s21_sum_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 1);
  ck_assert_double_eq(result.matrix[0][1], 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sum_matrix_5) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 3;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_sum_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 2);
  ck_assert_double_eq(result.matrix[1][0], 4);
  ck_assert_double_eq(result.matrix[2][0], 6);
  ck_assert_double_eq(result.matrix[3][0], 4);
  ck_assert_double_eq(result.matrix[4][0], 2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sum_matrix_6) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 5;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_sum_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 2);
  ck_assert_double_eq(result.matrix[1][0], 4);
  ck_assert_double_eq(result.matrix[2][0], 8);
  ck_assert_double_eq(result.matrix[3][0], 4);
  ck_assert_double_eq(result.matrix[4][0], 2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sum_matrix_8) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 24700;
  A.matrix[1][0] = -0.999;
  A.matrix[2][0] = 123.091;
  A.matrix[3][0] = 0.12355;
  A.matrix[4][0] = 0.11;
  B.matrix[0][0] = -24673;
  B.matrix[1][0] = -0.001;
  B.matrix[2][0] = -0.091;
  B.matrix[3][0] = 0.00145;
  B.matrix[4][0] = -0.11;
  int fail3 = s21_sum_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 27);
  ck_assert_double_eq(result.matrix[1][0], -1);
  ck_assert_double_eq(result.matrix[2][0], 123);
  ck_assert_double_eq(result.matrix[3][0], 0.125);
  ck_assert_double_eq(result.matrix[4][0], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix_2) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(1, 1, &B);
  int fail3 = s21_sub_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix_3) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int fail3 = s21_sub_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix_4) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][1] = 1;
  int fail3 = s21_sub_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 1);
  ck_assert_double_eq(result.matrix[0][1], -1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix_5) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 3;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_sub_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 0);
  ck_assert_double_eq(result.matrix[1][0], 0);
  ck_assert_double_eq(result.matrix[2][0], 0);
  ck_assert_double_eq(result.matrix[3][0], 0);
  ck_assert_double_eq(result.matrix[4][0], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix_6) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 5;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_sub_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 0);
  ck_assert_double_eq(result.matrix[1][0], 0);
  ck_assert_double_eq(result.matrix[2][0], -2);
  ck_assert_double_eq(result.matrix[3][0], 0);
  ck_assert_double_eq(result.matrix[4][0], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_sub_matrix_8) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 24700;
  A.matrix[1][0] = -0.999;
  A.matrix[2][0] = 123.091;
  A.matrix[3][0] = 0.12355;
  A.matrix[4][0] = 0.11;
  B.matrix[0][0] = -24673;
  B.matrix[1][0] = -0.001;
  B.matrix[2][0] = -0.091;
  B.matrix[3][0] = 0.00145;
  B.matrix[4][0] = -0.11;
  int fail3 = s21_sub_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 49373, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -0.998, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 123.182, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][0], 0.12210, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][0], 0.22, 1E-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_2) {
  matrix_t A;
  double temp1 = 0;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[1][1] = 1;
  A.matrix[1][0] = 1;
  matrix_t result;

  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 0);
  ck_assert_double_eq(result.matrix[0][1], 0);
  ck_assert_double_eq(result.matrix[1][1], 0);
  ck_assert_double_eq(result.matrix[1][0], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_3) {
  matrix_t A;
  double temp1 = 123.0123;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[1][1] = 1;
  A.matrix[1][0] = 1;
  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], temp1);
  ck_assert_double_eq(result.matrix[0][1], temp1);
  ck_assert_double_eq(result.matrix[1][1], temp1);
  ck_assert_double_eq(result.matrix[1][0], temp1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_4) {
  matrix_t A;
  double temp1 = 1234;
  matrix_t result;
  int fail1 = s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[0][2] = 1;
  A.matrix[1][0] = 1;
  A.matrix[1][1] = 1;
  A.matrix[1][2] = 1;
  A.matrix[2][0] = 1;
  A.matrix[2][1] = 1;
  A.matrix[2][2] = 1;
  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], temp1);
  ck_assert_double_eq(result.matrix[0][1], temp1);
  ck_assert_double_eq(result.matrix[0][2], temp1);
  ck_assert_double_eq(result.matrix[1][0], temp1);
  ck_assert_double_eq(result.matrix[1][1], temp1);
  ck_assert_double_eq(result.matrix[1][2], temp1);
  ck_assert_double_eq(result.matrix[2][0], temp1);
  ck_assert_double_eq(result.matrix[2][1], temp1);
  ck_assert_double_eq(result.matrix[2][2], temp1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_5) {
  matrix_t A;
  double temp1 = 0.1;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 0.1, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.2, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 0.3, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][0], 0.2, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][0], 0.1, 1E-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_6) {
  matrix_t A;
  double temp1 = 0.1;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  A.matrix[0][0] = 40;
  A.matrix[1][0] = 3.0;
  A.matrix[2][0] = 0.3;
  A.matrix[3][0] = 0.2;
  A.matrix[4][0] = 0.1;
  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 4, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 0.3, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 0.03, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][0], 0.02, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][0], 0.01, 1E-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_7) {
  matrix_t A;
  double temp1 = 0.00001;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  A.matrix[0][0] = 100000;
  A.matrix[1][0] = 200000;
  A.matrix[2][0] = 300000;
  A.matrix[3][0] = -200000;
  A.matrix[4][0] = -100000;
  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(result.matrix[0][0], 1);
  ck_assert_int_eq(result.matrix[1][0], 2);
  ck_assert_int_eq(result.matrix[2][0], 3);
  ck_assert_int_eq(result.matrix[3][0], -2);
  ck_assert_int_eq(result.matrix[4][0], -1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_number_8) {
  matrix_t A;
  double temp1 = 98.1;
  matrix_t result;
  int fail1 = s21_create_matrix(10, 1, &A);
  A.matrix[0][0] = 24700;
  A.matrix[1][0] = -0.999;
  A.matrix[2][0] = 123.091;
  A.matrix[3][0] = 0.12355;
  A.matrix[4][0] = 0.11;
  A.matrix[5][0] = -24673;
  A.matrix[6][0] = -0.001;
  A.matrix[7][0] = -0.091;
  A.matrix[8][0] = 0.00145;
  A.matrix[9][0] = -0.11;
  int fail2 = s21_mult_number(&A, temp1, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], 2423070, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][0], -98.0019, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][0], 12075.2271, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][0], 12.120255, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][0], 10.791, 1E-6);
  ck_assert_double_eq_tol(result.matrix[5][0], -2420421.3, 1E-6);
  ck_assert_double_eq_tol(result.matrix[6][0], -0.0981, 1E-6);
  ck_assert_double_eq_tol(result.matrix[7][0], -8.9271, 1E-6);
  ck_assert_double_eq_tol(result.matrix[8][0], 0.142245, 1E-6);
  ck_assert_double_eq_tol(result.matrix[9][0], -10.791, 1E-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix_2) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(1, 1, &B);
  int fail3 = s21_mult_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 2);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix_3) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int fail3 = s21_mult_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq(result.matrix[0][0], 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix_4) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;

  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;

  int fail3 = s21_mult_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_int_eq(result.matrix[0][0], 7);
  ck_assert_int_eq(result.matrix[0][1], 10);
  ck_assert_int_eq(result.matrix[1][0], 15);
  ck_assert_int_eq(result.matrix[1][1], 22);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix_5) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(1, 5, &B);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[0][2] = 3;
  B.matrix[0][3] = 2;
  B.matrix[0][4] = 1;
  int fail3 = s21_mult_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_int_eq(result.matrix[0][0], 1);
  ck_assert_int_eq(result.matrix[0][1], 2);
  ck_assert_int_eq(result.matrix[0][2], 3);
  ck_assert_int_eq(result.matrix[0][3], 2);
  ck_assert_int_eq(result.matrix[0][4], 1);
  ck_assert_int_eq(result.matrix[1][0], 2);
  ck_assert_int_eq(result.matrix[1][1], 4);
  ck_assert_int_eq(result.matrix[1][2], 6);
  ck_assert_int_eq(result.matrix[1][3], 4);
  ck_assert_int_eq(result.matrix[1][4], 2);
  ck_assert_int_eq(result.matrix[2][0], 3);
  ck_assert_int_eq(result.matrix[2][1], 6);
  ck_assert_int_eq(result.matrix[2][2], 9);
  ck_assert_int_eq(result.matrix[2][3], 6);
  ck_assert_int_eq(result.matrix[2][4], 3);
  ck_assert_int_eq(result.matrix[3][0], 2);
  ck_assert_int_eq(result.matrix[3][1], 4);
  ck_assert_int_eq(result.matrix[3][2], 6);
  ck_assert_int_eq(result.matrix[3][3], 4);
  ck_assert_int_eq(result.matrix[3][4], 2);
  ck_assert_int_eq(result.matrix[4][0], 1);
  ck_assert_int_eq(result.matrix[4][1], 2);
  ck_assert_int_eq(result.matrix[4][2], 3);
  ck_assert_int_eq(result.matrix[4][3], 2);
  ck_assert_int_eq(result.matrix[4][4], 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix_6) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(1, 5, &A);
  int fail2 = s21_create_matrix(5, 1, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = 2;
  A.matrix[0][4] = 1;
  B.matrix[0][0] = 1;
  B.matrix[1][0] = 2;
  B.matrix[2][0] = 5;
  B.matrix[3][0] = 2;
  B.matrix[4][0] = 1;
  int fail3 = s21_mult_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_int_eq(result.rows, 1);
  ck_assert_int_eq(result.columns, 1);
  ck_assert_int_eq(result.matrix[0][0], 25);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_mult_matrix_8) {
  matrix_t A;
  matrix_t B;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  int fail2 = s21_create_matrix(1, 5, &B);
  A.matrix[0][0] = 247;
  A.matrix[1][0] = -0.999;
  A.matrix[2][0] = 123.091;
  A.matrix[3][0] = 0.12355;
  A.matrix[4][0] = 0.11;
  B.matrix[0][0] = -246;
  B.matrix[0][1] = -0.001;
  B.matrix[0][2] = -0.091;
  B.matrix[0][3] = 0.00145;
  B.matrix[0][4] = -0.11;
  int fail3 = s21_mult_matrix(&A, &B, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(fail3, 0);
  ck_assert_double_eq_tol(result.matrix[0][0], -60762, 1E-6);
  ck_assert_double_eq_tol(result.matrix[0][1], -0.247, 1E-6);
  ck_assert_double_eq_tol(result.matrix[0][2], -22.477, 1E-6);
  ck_assert_double_eq_tol(result.matrix[0][3], 0.35815, 1E-6);
  ck_assert_double_eq_tol(result.matrix[0][4], -27.17, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][0], 245.754, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][1], 0.000999, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][2], 0.090909, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][3], -0.001448, 1E-6);
  ck_assert_double_eq_tol(result.matrix[1][4], 0.10989, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][0], -30280.386, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][1], -0.123091, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][2], -11.201281, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][3], 0.178481, 1E-6);
  ck_assert_double_eq_tol(result.matrix[2][4], -13.54001, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][0], -30.3933, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][1], -0.000123, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][2], -0.011243, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][3], 0.000179, 1E-6);
  ck_assert_double_eq_tol(result.matrix[3][4], -0.013590, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][0], -27.06, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][1], -0.00011, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][2], -0.01001, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][3], 0.000159, 1E-6);
  ck_assert_double_eq_tol(result.matrix[4][4], -0.0121, 1E-6);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose_2) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_transpose(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose_3) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 1;
  int fail2 = s21_transpose(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 1);
  ck_assert_double_eq(result.matrix[0][1], 0);
  ck_assert_double_eq(result.matrix[1][0], 1);
  ck_assert_double_eq(result.matrix[1][1], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose_4) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][1] = -1.1;
  int fail2 = s21_transpose(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][1], 0);
  ck_assert_double_eq(result.matrix[1][0], -1.1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose_5) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  int fail2 = s21_transpose(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 1);
  ck_assert_double_eq(result.matrix[0][1], 2);
  ck_assert_double_eq(result.matrix[0][2], 3);
  ck_assert_double_eq(result.matrix[0][3], 2);
  ck_assert_double_eq(result.matrix[0][4], 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose_6) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[1][0] = 2;
  A.matrix[2][0] = 3;
  A.matrix[3][0] = 2;
  A.matrix[4][0] = 1;
  A.matrix[0][1] = 1;
  A.matrix[1][1] = 2;
  A.matrix[2][1] = 3;
  A.matrix[3][1] = 2;
  A.matrix[4][1] = 1;
  int fail2 = s21_transpose(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 1);
  ck_assert_double_eq(result.matrix[0][1], 2);
  ck_assert_double_eq(result.matrix[0][2], 3);
  ck_assert_double_eq(result.matrix[0][3], 2);
  ck_assert_double_eq(result.matrix[0][4], 1);
  ck_assert_double_eq(result.matrix[1][0], 1);
  ck_assert_double_eq(result.matrix[1][1], 2);
  ck_assert_double_eq(result.matrix[1][2], 3);
  ck_assert_double_eq(result.matrix[1][3], 2);
  ck_assert_double_eq(result.matrix[1][4], 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_transpose_8) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(5, 1, &A);
  A.matrix[0][0] = 24700;
  A.matrix[1][0] = -0.999;
  A.matrix[2][0] = 123.091;
  A.matrix[3][0] = 0.12355;
  A.matrix[4][0] = 0.11;
  int fail2 = s21_transpose(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);

  ck_assert_double_eq(result.matrix[0][0], 24700);
  ck_assert_double_eq(result.matrix[0][1], -0.999);
  ck_assert_double_eq(result.matrix[0][2], 123.091);
  ck_assert_double_eq(result.matrix[0][3], 0.12355);
  ck_assert_double_eq(result.matrix[0][4], 0.11);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

// START_TEST(test_s21_determinant_1) {
//   matrix_t A;
//   double result = 0;
//   int fail1 = s21_create_matrix(5, 1, &A);
//   int fail2 = s21_determinant(&A, &result);
//   ck_assert_int_eq(fail1, 0);
//   ck_assert_int_eq(fail2, 2);
//   ck_assert_double_eq(result, 0);
//   s21_remove_matrix(&A);
// }
// END_TEST

START_TEST(test_s21_determinant_2) {
  matrix_t A;
  double result = 0;
  int fail1 = s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 2;
  int fail2 = s21_determinant(&A, &result);
  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result, 2);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_determinant_3) {
  matrix_t A;
  double result = 0;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 8;
  A.matrix[1][0] = 2.4;
  A.matrix[1][1] = -0.98765;
  int fail2 = s21_determinant(&A, &result);
  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result, -21.1753);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_determinant_4) {
  matrix_t A;
  double result = 0;
  int fail1 = s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 81.07652;
  A.matrix[0][1] = 0.931209;
  A.matrix[0][2] = 701;
  A.matrix[1][0] = 2.4;
  A.matrix[1][1] = -0.98765;
  A.matrix[1][2] = 32.09133;
  A.matrix[2][0] = 2.333;
  A.matrix[2][1] = -91.01;
  A.matrix[2][2] = 73.0001;
  int fail2 = s21_determinant(&A, &result);
  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_int_eq(result, 79355);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_determinant_5) {
  matrix_t A;
  double result = 0;
  int fail1 = s21_create_matrix(4, 4, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = -9.1;
  A.matrix[1][0] = 13.3;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 2;
  A.matrix[1][3] = 10;
  A.matrix[2][0] = 2;
  A.matrix[2][1] = -1;
  A.matrix[2][2] = -1;
  A.matrix[2][3] = -1;
  A.matrix[3][0] = 2;
  A.matrix[3][1] = -1.7;
  A.matrix[3][2] = -1.6;
  A.matrix[3][3] = -1.3;
  int fail2 = s21_determinant(&A, &result);
  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq_tol(result, -28.903, 1E-6);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(test_s21_determinant_6) {
  matrix_t matrix;
  s21_create_matrix(1, 1, &matrix);
  matrix.matrix[0][0] = 1.232;

  double det = 0;

  int res = s21_determinant(&matrix, &det);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(det, 1.232, 1E-6);
  s21_remove_matrix(&matrix);
}
END_TEST

START_TEST(test_s21_determinant_7) {
  matrix_t matrix;
  s21_create_matrix(2, 2, &matrix);
  matrix.matrix[0][0] = 2342.234;
  matrix.matrix[0][1] = 2.312;
  matrix.matrix[1][0] = 24.234;
  matrix.matrix[1][1] = 424.2;

  double det = 0;

  int res = s21_determinant(&matrix, &det);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(det, 993519.633792, 1E-6);
  s21_remove_matrix(&matrix);
}
END_TEST

START_TEST(test_s21_determinant_8) {
  matrix_t matrix;
  s21_create_matrix(3, 3, &matrix);
  matrix.matrix[0][0] = 2.34;
  matrix.matrix[0][1] = 27.345;
  matrix.matrix[0][2] = 25.23;
  matrix.matrix[1][0] = 2.55;
  matrix.matrix[1][1] = 8.56;
  matrix.matrix[1][2] = 45.75;
  matrix.matrix[2][0] = 55.3;
  matrix.matrix[2][1] = 4.34;
  matrix.matrix[2][2] = 24.33;

  double det = 0;

  int res = s21_determinant(&matrix, &det);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(det, 55844.5082595, 1E-6);
  s21_remove_matrix(&matrix);
}
END_TEST

START_TEST(test_s21_determinant_9) {
  matrix_t matrix;
  s21_create_matrix(4, 4, &matrix);
  matrix.matrix[0][0] = 1.34;
  matrix.matrix[0][1] = 7.423;
  matrix.matrix[0][2] = 2.45;
  matrix.matrix[0][3] = 2.42;
  matrix.matrix[1][0] = 34;
  matrix.matrix[1][1] = 2;
  matrix.matrix[1][2] = 5;
  matrix.matrix[1][3] = 24.52;
  matrix.matrix[2][0] = 3;
  matrix.matrix[2][1] = 4.5;
  matrix.matrix[2][2] = 2.4;
  matrix.matrix[2][3] = 2.4;
  matrix.matrix[3][0] = 8.43;
  matrix.matrix[3][1] = 3.5;
  matrix.matrix[3][2] = 53.6;
  matrix.matrix[3][3] = 24.5;

  double det = 0;

  int res = s21_determinant(&matrix, &det);
  ck_assert_int_eq(res, 0);
  ck_assert_double_eq_tol(det, -10769.66942328, 1E-6);
  s21_remove_matrix(&matrix);
}
END_TEST

START_TEST(test_s21_calc_complements_2) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  int fail2 = s21_calc_complements(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 0);
  ck_assert_double_eq(result.matrix[0][1], 0);
  ck_assert_double_eq(result.matrix[1][0], 0);
  ck_assert_double_eq(result.matrix[1][1], 0);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}

END_TEST

START_TEST(test_s21_calc_complements_3) {
  matrix_t test, test_2, result;
  s21_create_matrix(3, 3, &test);
  test.matrix[0][0] = 1;
  test.matrix[0][1] = 2;
  test.matrix[0][2] = 3;
  test.matrix[1][0] = 0;
  test.matrix[1][1] = 4;
  test.matrix[1][2] = 2;
  test.matrix[2][0] = 5;
  test.matrix[2][1] = 2;
  test.matrix[2][2] = 1;
  s21_calc_complements(&test, &result);
  s21_create_matrix(3, 3, &test_2);
  test_2.matrix[0][0] = 0;
  test_2.matrix[0][1] = 10;
  test_2.matrix[0][2] = -20;
  test_2.matrix[1][0] = 4;
  test_2.matrix[1][1] = -14;
  test_2.matrix[1][2] = 8;
  test_2.matrix[2][0] = -8;
  test_2.matrix[2][1] = -2;
  test_2.matrix[2][2] = 4;
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      ck_assert_double_eq(result.matrix[i][j], test_2.matrix[i][j]);
    }
  }
  s21_remove_matrix(&test);
  s21_remove_matrix(&test_2);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_calc_complements_4) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;

  int fail2 = s21_calc_complements(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 4);
  ck_assert_double_eq(result.matrix[0][1], -3);
  ck_assert_double_eq(result.matrix[1][0], -2);
  ck_assert_double_eq(result.matrix[1][1], 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_calc_complements_5) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(2, 2, &A);
  A.matrix[0][0] = 0.12;
  A.matrix[0][1] = -0.93;
  A.matrix[1][0] = -1.9;
  A.matrix[1][1] = 0.33;

  int fail2 = s21_calc_complements(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], 0.33);
  ck_assert_double_eq(result.matrix[0][1], 1.9);
  ck_assert_double_eq(result.matrix[1][0], 0.93);
  ck_assert_double_eq(result.matrix[1][1], 0.12);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_calc_complements_6) {
  matrix_t A;
  matrix_t result;
  int fail1 = s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = -1;
  A.matrix[1][2] = 0;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 5;

  int fail2 = s21_calc_complements(&A, &result);

  ck_assert_int_eq(fail1, 0);
  ck_assert_int_eq(fail2, 0);
  ck_assert_double_eq(result.matrix[0][0], -5);
  ck_assert_double_eq(result.matrix[0][1], -10);
  ck_assert_double_eq(result.matrix[0][2], 4);
  ck_assert_double_eq(result.matrix[1][0], -4);
  ck_assert_double_eq(result.matrix[1][1], 5);
  ck_assert_double_eq(result.matrix[1][2], -2);
  ck_assert_double_eq(result.matrix[2][0], 3);
  ck_assert_double_eq(result.matrix[2][1], 6);
  ck_assert_double_eq(result.matrix[2][2], -5);

  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(test_s21_inverse_matrix_1) {
  matrix_t matrix1 = {0};
  matrix_t result = {0};
  s21_create_matrix(1, 1, &matrix1);

  matrix1.matrix[0][0] = 0;
  ck_assert_int_eq(s21_inverse_matrix(&matrix1, &result), 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&result);
}

END_TEST

START_TEST(test_s21_inverse_matrix_2) {
  matrix_t A;
  matrix_t result;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &result);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 0);
  ck_assert_double_eq(result.matrix[0][0], 1);
  ck_assert_double_eq(result.matrix[0][1], -1);
  ck_assert_double_eq(result.matrix[0][2], 1);
  ck_assert_double_eq(result.matrix[1][0], -38);
  ck_assert_double_eq(result.matrix[1][1], 41);
  ck_assert_double_eq(result.matrix[1][2], -34);
  ck_assert_double_eq(result.matrix[2][0], 27);
  ck_assert_double_eq(result.matrix[2][1], -29);
  ck_assert_double_eq(result.matrix[2][2], 24);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}

END_TEST

START_TEST(test_s21_inverse_matrix_3) {
  matrix_t matrix1 = {0};
  matrix_t result = {0};
  s21_create_matrix(3, 3, &matrix1);

  matrix1.matrix[0][0] = 1;
  matrix1.matrix[0][1] = 2;
  matrix1.matrix[0][2] = 3;
  matrix1.matrix[1][0] = 4;
  matrix1.matrix[1][1] = 5;
  matrix1.matrix[1][2] = 6;
  matrix1.matrix[2][0] = 7;
  matrix1.matrix[2][1] = 8;
  matrix1.matrix[2][2] = 9;

  ck_assert_int_eq(s21_inverse_matrix(&matrix1, &result), 2);

  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&result);
}

END_TEST

START_TEST(test_s21_inverse_matrix_4) {
  matrix_t matrix1 = {0};
  matrix_t result = {0};
  s21_create_matrix(2, 3, &matrix1);

  matrix1.matrix[0][0] = 1;
  matrix1.matrix[0][1] = 2;
  matrix1.matrix[0][2] = 3;
  matrix1.matrix[1][0] = 0;
  matrix1.matrix[1][1] = 4;
  matrix1.matrix[1][2] = 2;
  ck_assert_int_eq(s21_inverse_matrix(&matrix1, &result), 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&result);
}

END_TEST

int main() {
  Suite *s1 = suite_create("s21_matrix");
  TCase *tc1_1 = tcase_create("s21_matrix");
  SRunner *sr = srunner_create(s1);
  int nf;

  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, test_s21_create_matrix_3);
  tcase_add_test(tc1_1, test_s21_create_matrix_4);
  tcase_add_test(tc1_1, test_s21_create_matrix_5);

  tcase_add_test(tc1_1, test_s21_eq_matrix_2);
  tcase_add_test(tc1_1, test_s21_eq_matrix_3);
  tcase_add_test(tc1_1, test_s21_eq_matrix_4);
  tcase_add_test(tc1_1, test_s21_eq_matrix_5);
  tcase_add_test(tc1_1, test_s21_eq_matrix_6);
  tcase_add_test(tc1_1, test_s21_eq_matrix_7);

  tcase_add_test(tc1_1, test_s21_sum_matrix_2);
  tcase_add_test(tc1_1, test_s21_sum_matrix_3);
  tcase_add_test(tc1_1, test_s21_sum_matrix_4);
  tcase_add_test(tc1_1, test_s21_sum_matrix_5);
  tcase_add_test(tc1_1, test_s21_sum_matrix_6);
  tcase_add_test(tc1_1, test_s21_sum_matrix_8);

  tcase_add_test(tc1_1, test_s21_sub_matrix_2);
  tcase_add_test(tc1_1, test_s21_sub_matrix_3);
  tcase_add_test(tc1_1, test_s21_sub_matrix_4);
  tcase_add_test(tc1_1, test_s21_sub_matrix_5);
  tcase_add_test(tc1_1, test_s21_sub_matrix_6);
  tcase_add_test(tc1_1, test_s21_sub_matrix_8);

  tcase_add_test(tc1_1, test_s21_mult_number_2);
  tcase_add_test(tc1_1, test_s21_mult_number_3);
  tcase_add_test(tc1_1, test_s21_mult_number_4);
  tcase_add_test(tc1_1, test_s21_mult_number_5);
  tcase_add_test(tc1_1, test_s21_mult_number_6);
  tcase_add_test(tc1_1, test_s21_mult_number_7);
  tcase_add_test(tc1_1, test_s21_mult_number_8);

  tcase_add_test(tc1_1, test_s21_mult_matrix_2);
  tcase_add_test(tc1_1, test_s21_mult_matrix_3);
  tcase_add_test(tc1_1, test_s21_mult_matrix_4);
  tcase_add_test(tc1_1, test_s21_mult_matrix_5);
  tcase_add_test(tc1_1, test_s21_mult_matrix_6);
  tcase_add_test(tc1_1, test_s21_mult_matrix_8);

  tcase_add_test(tc1_1, test_s21_transpose_2);
  tcase_add_test(tc1_1, test_s21_transpose_3);
  tcase_add_test(tc1_1, test_s21_transpose_4);
  tcase_add_test(tc1_1, test_s21_transpose_5);
  tcase_add_test(tc1_1, test_s21_transpose_6);
  tcase_add_test(tc1_1, test_s21_transpose_8);

  tcase_add_test(tc1_1, test_s21_determinant_2);
  tcase_add_test(tc1_1, test_s21_determinant_3);
  tcase_add_test(tc1_1, test_s21_determinant_4);
  tcase_add_test(tc1_1, test_s21_determinant_5);
  tcase_add_test(tc1_1, test_s21_determinant_6);
  tcase_add_test(tc1_1, test_s21_determinant_7);
  tcase_add_test(tc1_1, test_s21_determinant_8);
  tcase_add_test(tc1_1, test_s21_determinant_9);

  tcase_add_test(tc1_1, test_s21_calc_complements_2);
  tcase_add_test(tc1_1, test_s21_calc_complements_3);
  tcase_add_test(tc1_1, test_s21_calc_complements_4);
  tcase_add_test(tc1_1, test_s21_calc_complements_5);
  tcase_add_test(tc1_1, test_s21_calc_complements_6);

  tcase_add_test(tc1_1, test_s21_inverse_matrix_1);
  tcase_add_test(tc1_1, test_s21_inverse_matrix_2);
  tcase_add_test(tc1_1, test_s21_inverse_matrix_3);
  tcase_add_test(tc1_1, test_s21_inverse_matrix_4);

  srunner_run_all(sr, CK_VERBOSE);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}
