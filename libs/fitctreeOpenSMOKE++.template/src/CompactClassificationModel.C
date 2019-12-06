/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CompactClassificationModel.c
 *
 * Code generation for function 'CompactClassificationModel'
 *
 */

/* Include files */
#include "CompactClassificationModel.h"
#include "classifyPoint.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void c_CompactClassificationModel_se(c_classreg_learning_coder_class *obj)
{
  signed char b_I[2916];
  int k;
  memset(&b_I[0], 0, 2916U * sizeof(signed char));
  for (k = 0; k < 54; k++) {
    b_I[k + 54 * k] = 1;
  }

  for (k = 0; k < 2916; k++) {
    obj->Cost[k] = 1.0 - (double)b_I[k];
  }
}

/* End of code generation (CompactClassificationModel.c) */
