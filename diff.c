/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: diff.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 17-Nov-2023 16:44:14
 */

/* Include Files */
#include "diff.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const double x_data[]
 *                const int x_size[2]
 *                double y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
void diff(const double x_data[], const int x_size[2], double y_data[],
          int y_size[2])
{
  int dimSize;
  int u0;
  dimSize = x_size[1];
  u0 = x_size[1] - 1;
  if (u0 > 1) {
    u0 = 1;
  }
  if (u0 < 1) {
    y_size[0] = 1;
    y_size[1] = 0;
  } else {
    y_size[0] = 1;
    y_size[1] = x_size[1] - 1;
    if (x_size[1] - 1 != 0) {
      double work_data;
      work_data = x_data[0];
      for (u0 = 2; u0 <= dimSize; u0++) {
        double d;
        double tmp1;
        tmp1 = x_data[u0 - 1];
        d = tmp1;
        tmp1 -= work_data;
        work_data = d;
        y_data[u0 - 2] = tmp1;
      }
    }
  }
}

/*
 * File trailer for diff.c
 *
 * [EOF]
 */
