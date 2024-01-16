/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: find.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 17-Nov-2023 16:44:14
 */

/* Include Files */
#include "find.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : const boolean_T x[400]
 *                int i_data[]
 *                int i_size[2]
 * Return Type  : void
 */
void eml_find(const boolean_T x[400], int i_data[], int i_size[2])
{
  int idx;
  int ii;
  boolean_T exitg1;
  idx = 0;
  i_size[0] = 1;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii < 400)) {
    if (x[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= 400) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }
  if (idx < 1) {
    i_size[1] = 0;
  } else {
    i_size[1] = idx;
  }
}

/*
 * File trailer for find.c
 *
 * [EOF]
 */
