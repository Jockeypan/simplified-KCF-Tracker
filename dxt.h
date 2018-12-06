#ifndef _DXT_H_
#define _DXT_H_
#include "basetype.h"

namespace dxt{

enum {
	DFT_INVERSE = 1, DFT_SCALE = 2, DFT_ROWS = 4, DFT_COMPLEX_OUTPUT = 16, DFT_REAL_OUTPUT = 32,
	DCT_INVERSE = DFT_INVERSE, DCT_ROWS = DFT_ROWS
};
#define MAT_DXT_FORWARD  0
#define MAT_DXT_INVERSE  1
#define MAT_DXT_SCALE    2 /* divide result by size of array */
#define MAT_DXT_INV_SCALE (MAT_DXT_INVERSE + MAT_DXT_SCALE)
#define MAT_DXT_INVERSE_SCALE MAT_DXT_INV_SCALE
#define MAT_DXT_ROWS     4 /* transform each row individually */
#define MAT_DXT_MUL_CONJ 8 /* conjugate the second argument of cvMulSpectrums */

void dft(CMat src, CMat& dst, int flags = 0, int nonzeroRows = 0);
//void idft(CMat src, CMat dst, int flags = 0, int nonzeroRows = 0);
//void mulSpectrums(CMat a, CMat b, CMat c,
//	int flags, bool conjB = false);
}
#endif