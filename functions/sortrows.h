/*
 * sortrows.h
 *
 * Code generation for function 'sortrows'
 *
 */

#ifndef __SORTROWS_H__
#define __SORTROWS_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "FindDrainedBlock_types.h"

/* Function Declarations */
extern void sortrows(real_T y_data[], int32_T y_size[1], real_T ndx_data[],
                     int32_T ndx_size[1]);

#endif

/* End of code generation (sortrows.h) */
