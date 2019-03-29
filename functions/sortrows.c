/*
 * sortrows.c
 *
 * Code generation for function 'sortrows'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "FindDrainedBlock.h"
#include "sortrows.h"

/* Function Declarations */
static boolean_T eml_sort_le(const real_T v_data[], int32_T irow1, int32_T irow2);

/* Function Definitions */
static boolean_T eml_sort_le(const real_T v_data[], int32_T irow1, int32_T irow2)
{
  boolean_T p;
  boolean_T b4;
  p = true;
  if ((v_data[irow1 - 1] == v_data[irow2 - 1]) || (muDoubleScalarIsNaN
       (v_data[irow1 - 1]) && muDoubleScalarIsNaN(v_data[irow2 - 1]))) {
    b4 = true;
  } else {
    b4 = false;
  }

  if (!b4) {
    if ((v_data[irow1 - 1] <= v_data[irow2 - 1]) || muDoubleScalarIsNaN
        (v_data[irow2 - 1])) {
      p = true;
    } else {
      p = false;
    }
  }

  return p;
}

void sortrows(real_T y_data[], int32_T y_size[1], real_T ndx_data[], int32_T
              ndx_size[1])
{
  int32_T m;
  int32_T idx_data[4];
  int32_T k;
  int32_T b_m;
  int32_T idx0_data[4];
  int32_T i;
  int32_T j;
  int32_T pEnd;
  int32_T p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;
  real_T ycol_data[4];
  m = y_size[0];
  ndx_size[0] = y_size[0];
  for (k = 1; k <= y_size[0]; k++) {
    idx_data[k - 1] = k;
  }

  for (k = 1; k <= y_size[0] - 1; k += 2) {
    if (eml_sort_le(y_data, k, k + 1)) {
    } else {
      idx_data[k - 1] = k + 1;
      idx_data[k] = k;
    }
  }

  b_m = y_size[0];
  for (k = 0; k < b_m; k++) {
    idx0_data[k] = 1;
  }

  i = 2;
  while (i < y_size[0]) {
    b_m = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < y_size[0] + 1; pEnd = qEnd + i) {
      p = j;
      q = pEnd;
      qEnd = j + b_m;
      if (qEnd > y_size[0] + 1) {
        qEnd = y_size[0] + 1;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if (eml_sort_le(y_data, idx_data[p - 1], idx_data[q - 1])) {
          idx0_data[k] = idx_data[p - 1];
          p++;
          if (p == pEnd) {
            while (q < qEnd) {
              k++;
              idx0_data[k] = idx_data[q - 1];
              q++;
            }
          }
        } else {
          idx0_data[k] = idx_data[q - 1];
          q++;
          if (q == qEnd) {
            while (p < pEnd) {
              k++;
              idx0_data[k] = idx_data[p - 1];
              p++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx_data[(j + k) - 1] = idx0_data[k];
      }

      j = qEnd;
    }

    i = b_m;
  }

  b_m = y_size[0];
  for (i = 0; i + 1 <= y_size[0]; i++) {
    ycol_data[i] = y_data[idx_data[i] - 1];
  }

  for (i = 0; i + 1 <= b_m; i++) {
    y_data[i] = ycol_data[i];
  }

  for (k = 0; k + 1 <= m; k++) {
    ndx_data[k] = idx_data[k];
  }
}

/* End of code generation (sortrows.c) */
