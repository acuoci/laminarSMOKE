/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * classifyPoint_types.h
 *
 * Code generation for function 'classifyPoint_types'
 *
 */

#ifndef CLASSIFYPOINT_TYPES_H
#define CLASSIFYPOINT_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef enum_c_classreg_learning_coderutils_
#define enum_c_classreg_learning_coderutils_

enum c_classreg_learning_coderutils_
{
  Logit = 0,                           /* Default value */
  Doublelogit,
  Invlogit,
  Ismax,
  Sign,
  Symmetric,
  Symmetricismax,
  Symmetriclogit,
  Identity
};

#endif                                 /*enum_c_classreg_learning_coderutils_*/

#ifndef typedef_c_classreg_learning_coderutils_
#define typedef_c_classreg_learning_coderutils_

typedef enum c_classreg_learning_coderutils_ c_classreg_learning_coderutils_;

#endif                                 /*typedef_c_classreg_learning_coderutils_*/

#ifndef struct_emxArray_real_T_1273
#define struct_emxArray_real_T_1273

struct emxArray_real_T_1273
{
  double data[1273];
  int size[1];
};

#endif                                 /*struct_emxArray_real_T_1273*/

#ifndef typedef_emxArray_real_T_1273
#define typedef_emxArray_real_T_1273

typedef struct emxArray_real_T_1273 emxArray_real_T_1273;

#endif                                 /*typedef_emxArray_real_T_1273*/

#ifndef struct_sV1mTyf4SwJVgoO0e5HZf8G_tag
#define struct_sV1mTyf4SwJVgoO0e5HZf8G_tag

struct sV1mTyf4SwJVgoO0e5HZf8G_tag
{
  double CutPredictorIndex[1273];
  double Children[2546];
  double CutPoint[1273];
  emxArray_real_T_1273 PruneList;
  boolean_T NanCutPoints[1273];
  boolean_T InfCutPoints[1273];
  double ClassNames[54];
  int ClassNamesLength[54];
  c_classreg_learning_coderutils_ ScoreTransform;
  double Prior[54];
  boolean_T ClassLogicalIndices[54];
  double Cost[2916];
  double ClassProbability[68742];
};

#endif                                 /*struct_sV1mTyf4SwJVgoO0e5HZf8G_tag*/

#ifndef typedef_c_classreg_learning_coder_class
#define typedef_c_classreg_learning_coder_class

typedef struct sV1mTyf4SwJVgoO0e5HZf8G_tag c_classreg_learning_coder_class;

#endif                                 /*typedef_c_classreg_learning_coder_class*/
#endif

/* End of code generation (classifyPoint_types.h) */
