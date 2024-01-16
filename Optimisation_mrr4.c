/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Optimisation_mrr4.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 17-Nov-2023 16:44:14
 */

/* Include Files */
#include "Optimisation_mrr4.h"
#include "diff.h"
#include "find.h"
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static double rt_hypotd_snf(double u0, double u1);

static double rt_roundd_snf(double u);

/* Function Definitions */
/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

/*
 * Arguments    : double u
 * Return Type  : double
 */
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }
  return y;
}

/*
 * Arguments    : double neffK1
 *                double neffK2
 *                double neffK3
 *                double neffK4
 *                double t1
 *                double t2
 *                double t3
 *                double t4
 *                double t5
 * Return Type  : double
 */
double Optimisation_mrr4(double neffK1, double neffK2, double neffK3,
                         double neffK4, double t1, double t2, double t3,
                         double t4, double t5)
{
  creal_T b_z[4];
  creal_T z[4];
  double Td[400];
  double Td_data[400];
  double lamda_0[400];
  double tmp_data[399];
  double k[5];
  double t[5];
  double Goal_ftotal;
  double b_z_im;
  double b_z_re;
  double c_z_im;
  double c_z_re;
  double d;
  double d1;
  double d10;
  double d11;
  double d12;
  double d13;
  double d14;
  double d15;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d_z_im;
  double d_z_re;
  double delta1;
  double delta2;
  double e_z_im;
  double e_z_re;
  double m;
  double z_im;
  double z_re;
  int ind0_data[400];
  int ind30_data[400];
  int ind3_data[400];
  int Td_size[2];
  int ind30_size[2];
  int Goal_C;
  int b_k;
  int i;
  int i1;
  int x_size_idx_1;
  int y;
  boolean_T b_Td[400];
  boolean_T x_data[399];
  /*  四环级联优化：目标函数计算 */
  /* 利用for循环且考虑速度预分配节省算力 */
  /* 微环的环长 */
  /*  neffB=4.221;%等效折射率截距 */
  t[0] = t1;
  t[1] = t2;
  t[2] = t3;
  t[3] = t4;
  t[4] = t5;
  /* t=[round(t1,2), round(t2,2),round(t3,2), round(t4,2),round(t5,2)]; */
  /* t=[round(t1,3), round(t2,3),round(t3,3), round(t4,3),round(t5,3)]; */
  for (b_k = 0; b_k < 5; b_k++) {
    d = t[b_k];
    k[b_k] = sqrt(1.0 - d * d);
  }
  /* t为自耦合系数，k为交叉耦合系数 */
  /*  设定lamada扫描步长 */
  m = rt_roundd_snf(
      ((((-neffK1 * 1550.0 + 3.79312) * 1.0053096491487338E+6 / 1550.0 +
         (-neffK2 * 1550.0 + 3.79312) * 1.0053096491487338E+6 / 1550.0) +
        (-neffK3 * 1550.0 + 3.79312) * 1.0053096491487338E+6 / 1550.0) +
       (-neffK4 * 1550.0 + 3.79312) * 1.0053096491487338E+6 / 1550.0) /
      4.0);
  /* 谐振等级 */
  m = (((3.8132601363790454E+6 / (m + neffK1 * 1.0053096491487338E+6) +
         3.8132601363790454E+6 / (m + neffK2 * 1.0053096491487338E+6)) +
        3.8132601363790454E+6 / (m + neffK3 * 1.0053096491487338E+6)) +
       3.8132601363790454E+6 / (m + neffK4 * 1.0053096491487338E+6)) /
      4.0;
  /* 谐振峰中心位置 */
  /* 自由光谱程 */
  lamda_0[399] = m + 0.31501915868259384;
  lamda_0[0] = m - 0.31501915868259384;
  if (m - 0.31501915868259384 == -(m + 0.31501915868259384)) {
    m = (m + 0.31501915868259384) / 399.0;
    for (b_k = 0; b_k < 398; b_k++) {
      lamda_0[b_k + 1] = (2.0 * ((double)b_k + 2.0) - 401.0) * m;
    }
  } else if (((m - 0.31501915868259384 < 0.0) !=
              (m + 0.31501915868259384 < 0.0)) &&
             ((fabs(m - 0.31501915868259384) > 8.9884656743115785E+307) ||
              (fabs(m + 0.31501915868259384) > 8.9884656743115785E+307))) {
    delta1 = (m - 0.31501915868259384) / 399.0;
    delta2 = (m + 0.31501915868259384) / 399.0;
    for (b_k = 0; b_k < 398; b_k++) {
      lamda_0[b_k + 1] =
          ((m - 0.31501915868259384) + delta2 * ((double)b_k + 1.0)) -
          delta1 * ((double)b_k + 1.0);
    }
  } else {
    delta1 = ((m + 0.31501915868259384) - (m - 0.31501915868259384)) / 399.0;
    for (b_k = 0; b_k < 398; b_k++) {
      lamda_0[b_k + 1] =
          (m - 0.31501915868259384) + ((double)b_k + 1.0) * delta1;
    }
  }
  /* 设定lamada扫描步长 */
  /* 预分配数组 */
  /*  传输函数计算表达式 */
  m = 0.0 * k[0];
  if (k[0] == 0.0) {
    z_re = 1.0 / m;
    z_im = 0.0;
  } else if (m == 0.0) {
    z_re = 0.0;
    z_im = -(1.0 / k[0]);
  } else {
    z_re = rtNaN;
    z_im = rtNaN;
  }
  m = 0.0 * k[1];
  if (k[1] == 0.0) {
    b_z_re = 1.0 / m;
    b_z_im = 0.0;
  } else if (m == 0.0) {
    b_z_re = 0.0;
    b_z_im = -(1.0 / k[1]);
  } else {
    b_z_re = rtNaN;
    b_z_im = rtNaN;
  }
  m = 0.0 * k[2];
  if (k[2] == 0.0) {
    c_z_re = 1.0 / m;
    c_z_im = 0.0;
  } else if (m == 0.0) {
    c_z_re = 0.0;
    c_z_im = -(1.0 / k[2]);
  } else {
    c_z_re = rtNaN;
    c_z_im = rtNaN;
  }
  m = 0.0 * k[3];
  if (k[3] == 0.0) {
    d_z_re = 1.0 / m;
    d_z_im = 0.0;
  } else if (m == 0.0) {
    d_z_re = 0.0;
    d_z_im = -(1.0 / k[3]);
  } else {
    d_z_re = rtNaN;
    d_z_im = rtNaN;
  }
  m = 0.0 * k[4];
  if (k[4] == 0.0) {
    e_z_re = 1.0 / m;
    e_z_im = 0.0;
  } else if (m == 0.0) {
    e_z_re = 0.0;
    e_z_im = -(1.0 / k[4]);
  } else {
    e_z_re = rtNaN;
    e_z_im = rtNaN;
  }
  d = t4 * d_z_re;
  d1 = t4 * d_z_im;
  d2 = -t4 * d_z_re;
  d3 = -t4 * d_z_im;
  d4 = t3 * c_z_re;
  d5 = t3 * c_z_im;
  d6 = -t3 * c_z_re;
  d7 = -t3 * c_z_im;
  d8 = t2 * b_z_re;
  d9 = t2 * b_z_im;
  d10 = -t2 * b_z_re;
  d11 = -t2 * b_z_im;
  d12 = t1 * z_re;
  d13 = t1 * z_im;
  d14 = -t1 * z_re;
  d15 = -t1 * z_im;
  for (Goal_C = 0; Goal_C < 400; Goal_C++) {
    creal_T M[4];
    double M_im;
    double M_re;
    double b_x_im;
    double b_x_re;
    double c_x_im;
    double c_x_re;
    double d16;
    double d17;
    double d18;
    double d_x_im;
    double d_x_re;
    double e_x_im;
    double e_x_re;
    double f_x_im;
    double f_x_re;
    double g_x_im;
    double g_x_re;
    double neffK_idx_0;
    double x_im;
    double x_re;
    m = lamda_0[Goal_C];
    neffK_idx_0 = 6.31654681669719E+6 * (-neffK1 * m + 3.79312) / m;
    delta2 = 6.31654681669719E+6 * (-neffK2 * m + 3.79312) / m;
    delta1 = 6.31654681669719E+6 * (-neffK3 * m + 3.79312) / m;
    m = 6.31654681669719E+6 * (-neffK4 * m + 3.79312) / m;
    /* 环程相移 */
    /* lamda在1550nm附近，谐振波长也需要在这附近  */
    /* 耦合矩阵 */
    /* 环内传输矩阵 */
    M_re = 0.0 * m;
    if (-m == 0.0) {
      x_re = M_re / 2.0;
      x_im = 0.0;
    } else if (M_re == 0.0) {
      x_re = 0.0;
      x_im = -m / 2.0;
    } else {
      x_re = rtNaN;
      x_im = -m / 2.0;
    }
    if (x_re == 0.0) {
      x_re = cos(x_im);
      x_im = sin(x_im);
    } else if (x_im == 0.0) {
      x_re = rtNaN;
      x_im = 0.0;
    } else {
      x_re = rtNaN;
      x_im = rtNaN;
    }
    if (m == 0.0) {
      M_im = M_re / 2.0;
      m = 0.0;
    } else if (M_re == 0.0) {
      M_im = 0.0;
      m /= 2.0;
    } else {
      M_im = rtNaN;
      m /= 2.0;
    }
    if (M_im == 0.0) {
      M_im = cos(m);
      m = sin(m);
    } else if (m == 0.0) {
      M_im = rtNaN;
      m = 0.0;
    } else {
      M_im = rtNaN;
      m = rtNaN;
    }
    M_re = 0.0 * delta1;
    if (-delta1 == 0.0) {
      b_x_re = M_re / 2.0;
      b_x_im = 0.0;
    } else if (M_re == 0.0) {
      b_x_re = 0.0;
      b_x_im = -delta1 / 2.0;
    } else {
      b_x_re = rtNaN;
      b_x_im = -delta1 / 2.0;
    }
    if (b_x_re == 0.0) {
      b_x_re = cos(b_x_im);
      b_x_im = sin(b_x_im);
    } else if (b_x_im == 0.0) {
      b_x_re = rtNaN;
      b_x_im = 0.0;
    } else {
      b_x_re = rtNaN;
      b_x_im = rtNaN;
    }
    if (delta1 == 0.0) {
      e_x_re = M_re / 2.0;
      e_x_im = 0.0;
    } else if (M_re == 0.0) {
      e_x_re = 0.0;
      e_x_im = delta1 / 2.0;
    } else {
      e_x_re = rtNaN;
      e_x_im = delta1 / 2.0;
    }
    if (e_x_re == 0.0) {
      e_x_re = cos(e_x_im);
      e_x_im = sin(e_x_im);
    } else if (e_x_im == 0.0) {
      e_x_re = rtNaN;
      e_x_im = 0.0;
    } else {
      e_x_re = rtNaN;
      e_x_im = rtNaN;
    }
    M_re = 0.0 * delta2;
    if (-delta2 == 0.0) {
      c_x_re = M_re / 2.0;
      c_x_im = 0.0;
    } else if (M_re == 0.0) {
      c_x_re = 0.0;
      c_x_im = -delta2 / 2.0;
    } else {
      c_x_re = rtNaN;
      c_x_im = -delta2 / 2.0;
    }
    if (c_x_re == 0.0) {
      c_x_re = cos(c_x_im);
      c_x_im = sin(c_x_im);
    } else if (c_x_im == 0.0) {
      c_x_re = rtNaN;
      c_x_im = 0.0;
    } else {
      c_x_re = rtNaN;
      c_x_im = rtNaN;
    }
    if (delta2 == 0.0) {
      f_x_re = M_re / 2.0;
      f_x_im = 0.0;
    } else if (M_re == 0.0) {
      f_x_re = 0.0;
      f_x_im = delta2 / 2.0;
    } else {
      f_x_re = rtNaN;
      f_x_im = delta2 / 2.0;
    }
    if (f_x_re == 0.0) {
      f_x_re = cos(f_x_im);
      f_x_im = sin(f_x_im);
    } else if (f_x_im == 0.0) {
      f_x_re = rtNaN;
      f_x_im = 0.0;
    } else {
      f_x_re = rtNaN;
      f_x_im = rtNaN;
    }
    M_re = 0.0 * neffK_idx_0;
    if (-neffK_idx_0 == 0.0) {
      d_x_re = M_re / 2.0;
      d_x_im = 0.0;
    } else if (M_re == 0.0) {
      d_x_re = 0.0;
      d_x_im = -neffK_idx_0 / 2.0;
    } else {
      d_x_re = rtNaN;
      d_x_im = -neffK_idx_0 / 2.0;
    }
    if (d_x_re == 0.0) {
      d_x_re = cos(d_x_im);
      d_x_im = sin(d_x_im);
    } else if (d_x_im == 0.0) {
      d_x_re = rtNaN;
      d_x_im = 0.0;
    } else {
      d_x_re = rtNaN;
      d_x_im = rtNaN;
    }
    if (neffK_idx_0 == 0.0) {
      g_x_re = M_re / 2.0;
      g_x_im = 0.0;
    } else if (M_re == 0.0) {
      g_x_re = 0.0;
      g_x_im = neffK_idx_0 / 2.0;
    } else {
      g_x_re = rtNaN;
      g_x_im = neffK_idx_0 / 2.0;
    }
    if (g_x_re == 0.0) {
      g_x_re = cos(g_x_im);
      g_x_im = sin(g_x_im);
    } else if (g_x_im == 0.0) {
      g_x_re = rtNaN;
      g_x_im = 0.0;
    } else {
      g_x_re = rtNaN;
      g_x_im = rtNaN;
    }
    z[0].re = t5 * e_z_re;
    z[0].im = t5 * e_z_im;
    z[2].re = -e_z_re;
    z[2].im = -e_z_im;
    z[1].re = e_z_re;
    z[1].im = e_z_im;
    z[3].re = -t5 * e_z_re;
    z[3].im = -t5 * e_z_im;
    d16 = 1.0015033834597085 * M_im;
    d17 = 1.0015033834597085 * m;
    d18 = 0.99849887330932929 * x_re;
    delta1 = 0.99849887330932929 * x_im;
    for (i = 0; i < 2; i++) {
      delta2 = z[i].re;
      M_re = delta2 * 0.0;
      M_im = z[i].im;
      neffK_idx_0 = M_im * 0.0;
      x_re = z[i + 2].re;
      x_im = z[i + 2].im;
      b_z[i].re = (M_re - neffK_idx_0) + (x_re * d16 - x_im * d17);
      b_z[i].im = (M_re + neffK_idx_0) + (x_re * d17 + x_im * d16);
      M_re = x_re * 0.0;
      neffK_idx_0 = x_im * 0.0;
      b_z[i + 2].re = (delta2 * d18 - M_im * delta1) + (M_re - neffK_idx_0);
      b_z[i + 2].im = (delta2 * delta1 + M_im * d18) + (M_re + neffK_idx_0);
    }
    for (i = 0; i < 2; i++) {
      M_re = b_z[i].re;
      neffK_idx_0 = b_z[i].im;
      x_re = b_z[i + 2].re;
      x_im = b_z[i + 2].im;
      M[i].re = (M_re * d - neffK_idx_0 * d1) + (x_re * d_z_re - x_im * d_z_im);
      M[i].im = (M_re * d1 + neffK_idx_0 * d) + (x_re * d_z_im + x_im * d_z_re);
      M[i + 2].re =
          (M_re * -d_z_re - neffK_idx_0 * -d_z_im) + (x_re * d2 - x_im * d3);
      M[i + 2].im =
          (M_re * -d_z_im + neffK_idx_0 * -d_z_re) + (x_re * d3 + x_im * d2);
    }
    d16 = 1.0015033834597085 * e_x_re;
    d17 = 1.0015033834597085 * e_x_im;
    d18 = 0.99849887330932929 * b_x_re;
    delta1 = 0.99849887330932929 * b_x_im;
    for (i = 0; i < 2; i++) {
      delta2 = M[i].re;
      M_re = delta2 * 0.0;
      M_im = M[i].im;
      neffK_idx_0 = M_im * 0.0;
      x_re = M[i + 2].re;
      x_im = M[i + 2].im;
      m = (M_re - neffK_idx_0) + (x_re * d16 - x_im * d17);
      M_re = (M_re + neffK_idx_0) + (x_re * d17 + x_im * d16);
      neffK_idx_0 = x_re * 0.0;
      x_re = x_im * 0.0;
      x_im = (delta2 * d18 - M_im * delta1) + (neffK_idx_0 - x_re);
      neffK_idx_0 = (delta2 * delta1 + M_im * d18) + (neffK_idx_0 + x_re);
      M[i].re = (m * d4 - M_re * d5) + (x_im * c_z_re - neffK_idx_0 * c_z_im);
      M[i].im = (m * d5 + M_re * d4) + (x_im * c_z_im + neffK_idx_0 * c_z_re);
      M[i + 2].re =
          (m * -c_z_re - M_re * -c_z_im) + (x_im * d6 - neffK_idx_0 * d7);
      M[i + 2].im =
          (m * -c_z_im + M_re * -c_z_re) + (x_im * d7 + neffK_idx_0 * d6);
    }
    d16 = 1.0015033834597085 * f_x_re;
    d17 = 1.0015033834597085 * f_x_im;
    d18 = 0.99849887330932929 * c_x_re;
    delta1 = 0.99849887330932929 * c_x_im;
    for (i = 0; i < 2; i++) {
      delta2 = M[i].re;
      M_re = delta2 * 0.0;
      M_im = M[i].im;
      neffK_idx_0 = M_im * 0.0;
      x_re = M[i + 2].re;
      x_im = M[i + 2].im;
      m = (M_re - neffK_idx_0) + (x_re * d16 - x_im * d17);
      M_re = (M_re + neffK_idx_0) + (x_re * d17 + x_im * d16);
      neffK_idx_0 = x_re * 0.0;
      x_re = x_im * 0.0;
      x_im = (delta2 * d18 - M_im * delta1) + (neffK_idx_0 - x_re);
      neffK_idx_0 = (delta2 * delta1 + M_im * d18) + (neffK_idx_0 + x_re);
      M[i].re = (m * d8 - M_re * d9) + (x_im * b_z_re - neffK_idx_0 * b_z_im);
      M[i].im = (m * d9 + M_re * d8) + (x_im * b_z_im + neffK_idx_0 * b_z_re);
      M[i + 2].re =
          (m * -b_z_re - M_re * -b_z_im) + (x_im * d10 - neffK_idx_0 * d11);
      M[i + 2].im =
          (m * -b_z_im + M_re * -b_z_re) + (x_im * d11 + neffK_idx_0 * d10);
    }
    d16 = 1.0015033834597085 * g_x_re;
    d17 = 1.0015033834597085 * g_x_im;
    d18 = 0.99849887330932929 * d_x_re;
    delta1 = 0.99849887330932929 * d_x_im;
    for (i = 0; i < 2; i++) {
      delta2 = M[i].re;
      M_re = delta2 * 0.0;
      M_im = M[i].im;
      neffK_idx_0 = M_im * 0.0;
      x_re = M[i + 2].re;
      x_im = M[i + 2].im;
      z[i].re = (M_re - neffK_idx_0) + (x_re * d16 - x_im * d17);
      z[i].im = (M_re + neffK_idx_0) + (x_re * d17 + x_im * d16);
      M_re = x_re * 0.0;
      neffK_idx_0 = x_im * 0.0;
      z[i + 2].re = (delta2 * d18 - M_im * delta1) + (M_re - neffK_idx_0);
      z[i + 2].im = (delta2 * delta1 + M_im * d18) + (M_re + neffK_idx_0);
    }
    for (i = 0; i < 2; i++) {
      M_re = z[i].re;
      neffK_idx_0 = z[i].im;
      x_re = z[i + 2].re;
      x_im = z[i + 2].im;
      M[i].re = (M_re * d12 - neffK_idx_0 * d13) + (x_re * z_re - x_im * z_im);
      M[i].im = (M_re * d13 + neffK_idx_0 * d12) + (x_re * z_im + x_im * z_re);
      M[i + 2].re =
          (M_re * -z_re - neffK_idx_0 * -z_im) + (x_re * d14 - x_im * d15);
      M[i + 2].im =
          (M_re * -z_im + neffK_idx_0 * -z_re) + (x_re * d15 + x_im * d14);
    }
    /* 传输矩阵 */
    neffK_idx_0 = M[0].re * M[3].re - M[0].im * M[3].im;
    M_im = M[0].re * M[3].im + M[0].im * M[3].re;
    if (M[1].im == 0.0) {
      if (M_im == 0.0) {
        M_re = neffK_idx_0 / M[1].re;
        M_im = 0.0;
      } else if (neffK_idx_0 == 0.0) {
        M_re = 0.0;
        M_im /= M[1].re;
      } else {
        M_re = neffK_idx_0 / M[1].re;
        M_im /= M[1].re;
      }
    } else if (M[1].re == 0.0) {
      if (neffK_idx_0 == 0.0) {
        M_re = M_im / M[1].im;
        M_im = 0.0;
      } else if (M_im == 0.0) {
        M_re = 0.0;
        M_im = -(neffK_idx_0 / M[1].im);
      } else {
        M_re = M_im / M[1].im;
        M_im = -(neffK_idx_0 / M[1].im);
      }
    } else {
      delta2 = fabs(M[1].re);
      m = fabs(M[1].im);
      if (delta2 > m) {
        m = M[1].im / M[1].re;
        delta1 = M[1].re + m * M[1].im;
        M_re = (neffK_idx_0 + m * M_im) / delta1;
        M_im = (M_im - m * neffK_idx_0) / delta1;
      } else if (m == delta2) {
        if (M[1].re > 0.0) {
          m = 0.5;
        } else {
          m = -0.5;
        }
        if (M[1].im > 0.0) {
          delta1 = 0.5;
        } else {
          delta1 = -0.5;
        }
        M_re = (neffK_idx_0 * m + M_im * delta1) / delta2;
        M_im = (M_im * m - neffK_idx_0 * delta1) / delta2;
      } else {
        m = M[1].re / M[1].im;
        delta1 = M[1].im + m * M[1].re;
        M_re = (m * neffK_idx_0 + M_im) / delta1;
        M_im = (m * M_im - neffK_idx_0) / delta1;
      }
    }
    m = rt_hypotd_snf(M[2].re - M_re, M[2].im - M_im);
    Td[Goal_C] = 10.0 * log10(m * m);
    /* 下载端DROP传递函数 */
  }
  /* figure(100);plot(lamda_0,Td_0,'LineWidth',1.5,'Color',[1 0 0]); */
  /* text(wave01+fsr/20,min(Td_0)+5,['t_{1}=',num2str(t(1)),'
   * t_{2}=',num2str(t(2)),' t_{3}=',num2str(t(3)),' t_{4}=',num2str(t(4)),'
   * t_{5}=',num2str(t(5))]); */
  /*  构建评价函数及定义传输特性参数  */
  if (!rtIsNaN(Td[0])) {
    Goal_C = 1;
  } else {
    boolean_T exitg1;
    Goal_C = 0;
    b_k = 2;
    exitg1 = false;
    while ((!exitg1) && (b_k < 401)) {
      if (!rtIsNaN(Td[b_k - 1])) {
        Goal_C = b_k;
        exitg1 = true;
      } else {
        b_k++;
      }
    }
  }
  if (Goal_C == 0) {
    delta2 = Td[0];
  } else {
    delta2 = Td[Goal_C - 1];
    i = Goal_C + 1;
    for (b_k = i; b_k < 401; b_k++) {
      d = Td[b_k - 1];
      if (delta2 < d) {
        delta2 = d;
      }
    }
  }
  /* 计算通带带内波动 */
  for (i = 0; i < 400; i++) {
    b_Td[i] = (Td[i] >= delta2 - 0.5);
  }
  eml_find(b_Td, ind30_data, ind30_size);
  b_k = ind30_size[1];
  if (b_k - 1 >= 0) {
    memcpy(&ind0_data[0], &ind30_data[0], (unsigned int)b_k * sizeof(int));
  }
  /*  找到第一个功率大于等于0.5dB降低的频率 */
  /* 找到最后一个功率大于等于0.5dB降低的频率 */
  i = ind0_data[ind30_size[1] - 1];
  if (ind0_data[0] > i) {
    i1 = 0;
    Goal_C = 0;
  } else {
    i1 = ind0_data[0] - 1;
    Goal_C = i;
  }
  Td_size[0] = 1;
  b_k = Goal_C - i1;
  Td_size[1] = b_k;
  for (Goal_C = 0; Goal_C < b_k; Goal_C++) {
    Td_data[Goal_C] = Td[i1 + Goal_C];
  }
  m = maximum(Td_data, Td_size, &Goal_C);
  if (ind0_data[0] > i) {
    i1 = 0;
    Goal_C = 0;
  } else {
    i1 = ind0_data[0] - 1;
    Goal_C = i;
  }
  Td_size[0] = 1;
  b_k = Goal_C - i1;
  Td_size[1] = b_k;
  for (Goal_C = 0; Goal_C < b_k; Goal_C++) {
    Td_data[Goal_C] = Td[i1 + Goal_C];
  }
  delta1 = minimum(Td_data, Td_size, &Goal_C);
  /*  计算3dB通带 */
  for (i1 = 0; i1 < 400; i1++) {
    b_Td[i1] = (Td[i1] >= delta2 - 3.0);
  }
  eml_find(b_Td, ind30_data, ind30_size);
  Goal_C = ind30_size[1];
  b_k = ind30_size[1];
  if (b_k - 1 >= 0) {
    memcpy(&ind3_data[0], &ind30_data[0], (unsigned int)b_k * sizeof(int));
  }
  /*  找到第一个截止频率 */
  /*  找到第二个截止频率 */
  /*  计算10dB通带 */
  /*  找到第一个截止频率 */
  /*  计算30dB通带 */
  for (i1 = 0; i1 < 400; i1++) {
    b_Td[i1] = (Td[i1] >= delta2 - 30.0);
  }
  eml_find(b_Td, ind30_data, ind30_size);
  /*  找到第一个截止频率 */
  /*  找到第二个截止频率 */
  /* TndB= Td_min+10; */
  /*  找到第一个截止频率 */
  /*  找到第二个截止频率 */
  /*  构建目标函数 */
  /* 插损 */
  /* 构建目标函数_带内波动 */
  m = fabs(m - delta1);
  /* 构建目标函数_抑制比 */
  /* 构建目标函数_3dB带宽 */
  delta1 = fabs(lamda_0[ind3_data[Goal_C - 1] - 1] - lamda_0[ind3_data[0] - 1]);
  /* 目标函数：10dB带宽 */
  /* 目标函数：30dB带宽 */
  delta2 = fabs(lamda_0[ind30_data[ind30_size[1] - 1] - 1] -
                lamda_0[ind30_data[0] - 1]);
  /* 目标函数：阻带带宽 */
  /* 目标函数：形状因子 */
  /* 目标函数：Q参数 */
  /* 目标函数：传输谱平滑程度 */
  b_k = ind0_data[0];
  Td_size[0] = 1;
  Td_size[1] = ind0_data[0];
  if (b_k - 1 >= 0) {
    memcpy(&Td_data[0], &Td[0], (unsigned int)b_k * sizeof(double));
  }
  diff(Td_data, Td_size, tmp_data, ind30_size);
  x_size_idx_1 = ind30_size[1];
  b_k = ind30_size[1];
  for (i1 = 0; i1 < b_k; i1++) {
    x_data[i1] = (tmp_data[i1] < 0.0);
  }
  if (ind30_size[1] == 0) {
    y = 0;
  } else {
    y = x_data[0];
    for (b_k = 2; b_k <= x_size_idx_1; b_k++) {
      y += x_data[b_k - 1];
    }
  }
  Td_size[0] = 1;
  Td_size[1] = 401 - i;
  b_k = 401 - i;
  for (i1 = 0; i1 < b_k; i1++) {
    Td_data[i1] = Td[(i + i1) - 1];
  }
  diff(Td_data, Td_size, tmp_data, ind30_size);
  x_size_idx_1 = ind30_size[1];
  b_k = ind30_size[1];
  for (i = 0; i < b_k; i++) {
    x_data[i] = (tmp_data[i] > 0.0);
  }
  if (ind30_size[1] == 0) {
    Goal_C = 0;
  } else {
    Goal_C = x_data[0];
    for (b_k = 2; b_k <= x_size_idx_1; b_k++) {
      Goal_C += x_data[b_k - 1];
    }
  }
  Goal_C += y;
  /* 计算总的目标函数 */
  if (Goal_C > 2) {
    Goal_ftotal = ((fabs(delta1 - 0.0079) + fabs(delta2 - 0.049)) + m * 10.0) +
                  (double)Goal_C;
  } else if (m > 0.6) {
    Goal_ftotal = (fabs(delta1 - 0.0079) + fabs(delta2 - 0.049)) + m * 10.0;
  } else {
    Goal_ftotal = fabs(delta1 - 0.0079) + fabs(delta2 - 0.049);
  }
  /*  Goal_ftotal=2*Goal_fw_f*100+5*Goal_fbd*10+10*Goal_10dB/10+Goal_fc*8+Goal_fy/60;
   */
  /*  Goal_ftotal=Goal_3dBw^2+Goal_30dBw^2+Goal_10dBw^2; */
  return Goal_ftotal;
}

/*
 * File trailer for Optimisation_mrr4.c
 *
 * [EOF]
 */
