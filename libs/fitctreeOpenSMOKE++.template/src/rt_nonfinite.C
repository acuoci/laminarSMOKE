/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rt_nonfinite.c
 *
 * Code generation for function 'classifyPoint'
 *
 */

/*
 * Abstract:
 *      MATLAB for code generation function to initialize non-finites,
 *      (Inf, NaN and -Inf).
 */
#include "rt_nonfinite.h"
#include <math.h>

real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;

/* Function: rt_InitInfAndNaN ==================================================
 * Abstract:
 * Initialize the rtInf, rtMinusInf, and rtNaN needed by the
 * generated code. NaN is initialized as non-signaling. Assumes IEEE.
 */

/* Suppress Visual Studio 2013 INFINITY macro expansion compiler warning. */
#if defined(_MSC_VER) && _MSC_VER == 1800

#pragma warning(disable: 4756 56)

#endif

void rt_InitInfAndNaN()
{
  rtNaN = nan("");
  rtNaNF = nanf("");
  rtInf = (real_T)INFINITY;
  rtInfF = (real32_T)INFINITY;
  rtMinusInf = -(real_T)INFINITY;
  rtMinusInfF = -(real32_T)INFINITY;

#if defined(_MSC_VER) && _MSC_VER == 1800

#pragma warning(default: 4756 56)

#endif

}

/* Function: rtIsInf ==================================================
 * Abstract:
 * Test if value is infinite
 */
boolean_T rtIsInf(real_T value)
{
  return (isinf(value) ? 1U : 0U);
}

/* Function: rtIsInfF =================================================
 * Abstract:
 * Test if single-precision value is infinite
 */
boolean_T rtIsInfF(real32_T value)
{
  return (isinf((real_T)value) ? 1U : 0U);
}

/* Function: rtIsNaN ==================================================
 * Abstract:
 * Test if value is not a number
 */
boolean_T rtIsNaN(real_T value)
{
  return (isnan(value) ? 1U : 0U);
}

/* Function: rtIsNaNF =================================================
 * Abstract:
 * Test if single-precision value is not a number
 */
boolean_T rtIsNaNF(real32_T value)
{
  return (isnan((real_T)value) ? 1U : 0U);
}

/* End of code generation (rt_nonfinite.c) */
