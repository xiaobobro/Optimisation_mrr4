/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: minOrMax.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 17-Nov-2023 16:44:14
 */

/* Include Files */
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const double x_data[]
 *                const int x_size[2]
 *                int *idx
 * Return Type  : double
 */
double maximum(const double x_data[], const int x_size[2], int *idx)
{
  double ex;
  int k;
  int last;
  last = x_size[1];
  if (x_size[1] <= 2) {
    if (x_size[1] == 1) {
      ex = x_data[0];
      *idx = 1;
    } else {
      ex = x_data[x_size[1] - 1];
      if ((x_data[0] < ex) || (rtIsNaN(x_data[0]) && (!rtIsNaN(ex)))) {
        *idx = x_size[1];
      } else {
        ex = x_data[0];
        *idx = 1;
      }
    }
  } else {
    if (!rtIsNaN(x_data[0])) {
      *idx = 1;
    } else {
      boolean_T exitg1;
      *idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!rtIsNaN(x_data[k - 1])) {
          *idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (*idx == 0) {
      ex = x_data[0];
      *idx = 1;
    } else {
      int i;
      ex = x_data[*idx - 1];
      i = *idx + 1;
      for (k = i; k <= last; k++) {
        double d;
        d = x_data[k - 1];
        if (ex < d) {
          ex = d;
          *idx = k;
        }
      }
    }
  }
  return ex;
}

/*
 * Arguments    : const double x_data[]
 *                const int x_size[2]
 *                int *idx
 * Return Type  : double
 */
double minimum(const double x_data[], const int x_size[2], int *idx)
{
  double ex;
  int k;
  int last;
  last = x_size[1];
  if (x_size[1] <= 2) {
    if (x_size[1] == 1) {
      ex = x_data[0];
      *idx = 1;
    } else {
      ex = x_data[x_size[1] - 1];
      if ((x_data[0] > ex) || (rtIsNaN(x_data[0]) && (!rtIsNaN(ex)))) {
        *idx = x_size[1];
      } else {
        ex = x_data[0];
        *idx = 1;
      }
    }
  } else {
    if (!rtIsNaN(x_data[0])) {
      *idx = 1;
    } else {
      boolean_T exitg1;
      *idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!rtIsNaN(x_data[k - 1])) {
          *idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (*idx == 0) {
      ex = x_data[0];
      *idx = 1;
    } else {
      int i;
      ex = x_data[*idx - 1];
      i = *idx + 1;
      for (k = i; k <= last; k++) {
        double d;
        d = x_data[k - 1];
        if (ex > d) {
          ex = d;
          *idx = k;
        }
      }
    }
  }
  return ex;
}

/*
 * File trailer for minOrMax.c
 *
 * [EOF]
 */
