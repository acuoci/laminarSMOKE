/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CompactClassificationTree.c
 *
 * Code generation for function 'CompactClassificationTree'
 *
 */

/* Include files */
#include "CompactClassificationTree.h"
#include "classifyPoint.h"
#include "rt_nonfinite.h"

/* Function Definitions */
double c_CompactClassificationTree_pre(const double obj_CutPredictorIndex[1273],
  const double obj_Children[2546], const double obj_CutPoint[1273], const double
  obj_PruneList_data[], const boolean_T obj_NanCutPoints[1273], const double
  obj_ClassNames[54], const double obj_Cost[2916], const double
  obj_ClassProbability[68742], const double X[21])
{
  int m;
  int i;
  double unusedU4[54];
  double d;
  int k;
  boolean_T exitg1;
  double ex;
  m = 0;
  while (!((obj_PruneList_data[m] <= 0.0) || rtIsNaN(X[(int)
           obj_CutPredictorIndex[m] - 1]) || obj_NanCutPoints[m])) {
    if (X[(int)obj_CutPredictorIndex[m] - 1] < obj_CutPoint[m]) {
      m = (int)obj_Children[m << 1] - 1;
    } else {
      m = (int)obj_Children[(m << 1) + 1] - 1;
    }
  }

  for (i = 0; i < 54; i++) {
    d = 0.0;
    for (k = 0; k < 54; k++) {
      d += obj_ClassProbability[m + 1273 * k] * obj_Cost[k + 54 * i];
    }

    unusedU4[i] = d;
  }

  if (!rtIsNaN(unusedU4[0])) {
    m = 1;
  } else {
    m = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 55)) {
      if (!rtIsNaN(unusedU4[k - 1])) {
        m = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (m == 0) {
    m = 1;
  } else {
    ex = unusedU4[m - 1];
    i = m + 1;
    for (k = i; k < 55; k++) {
      d = unusedU4[k - 1];
      if (ex > d) {
        ex = d;
        m = k;
      }
    }
  }

  return obj_ClassNames[m - 1];
}

/* End of code generation (CompactClassificationTree.c) */
