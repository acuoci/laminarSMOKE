/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * classifyPoint_initialize.c
 *
 * Code generation for function 'classifyPoint_initialize'
 *
 */

/* Include files */
#include "classifyPoint_initialize.h"
#include "classifyPoint.h"
#include "classifyPoint_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void classifyPoint_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_classifyPoint = true;
}

/* End of code generation (classifyPoint_initialize.c) */
