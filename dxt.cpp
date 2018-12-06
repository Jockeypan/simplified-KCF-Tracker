/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "basetype.h"
#include "complex.hpp"
#include "dxt.h"

namespace dxt
{
#undef USE_IPP_DFT



/****************************************************************************************\
                               Discrete Fourier Transform
\****************************************************************************************/

//#define CV_MAX_LOCAL_DFT_SIZE  (1 << 15)
#define DXT_PI   3.1415926535897932384626433832795
#define DXT_SWAP(a,b,t) ((t) = (a), (a) = (b), (b) = (t))

void* dxtAlignPtr(const void* ptr, int align = 32)
{
	assert((align & (align - 1)) == 0);
	return (void*)(((size_t)ptr + align - 1) & ~(size_t)(align - 1));
}

static unsigned char bitrevTab[] =
{
  0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0,
  0x08,0x88,0x48,0xc8,0x28,0xa8,0x68,0xe8,0x18,0x98,0x58,0xd8,0x38,0xb8,0x78,0xf8,
  0x04,0x84,0x44,0xc4,0x24,0xa4,0x64,0xe4,0x14,0x94,0x54,0xd4,0x34,0xb4,0x74,0xf4,
  0x0c,0x8c,0x4c,0xcc,0x2c,0xac,0x6c,0xec,0x1c,0x9c,0x5c,0xdc,0x3c,0xbc,0x7c,0xfc,
  0x02,0x82,0x42,0xc2,0x22,0xa2,0x62,0xe2,0x12,0x92,0x52,0xd2,0x32,0xb2,0x72,0xf2,
  0x0a,0x8a,0x4a,0xca,0x2a,0xaa,0x6a,0xea,0x1a,0x9a,0x5a,0xda,0x3a,0xba,0x7a,0xfa,
  0x06,0x86,0x46,0xc6,0x26,0xa6,0x66,0xe6,0x16,0x96,0x56,0xd6,0x36,0xb6,0x76,0xf6,
  0x0e,0x8e,0x4e,0xce,0x2e,0xae,0x6e,0xee,0x1e,0x9e,0x5e,0xde,0x3e,0xbe,0x7e,0xfe,
  0x01,0x81,0x41,0xc1,0x21,0xa1,0x61,0xe1,0x11,0x91,0x51,0xd1,0x31,0xb1,0x71,0xf1,
  0x09,0x89,0x49,0xc9,0x29,0xa9,0x69,0xe9,0x19,0x99,0x59,0xd9,0x39,0xb9,0x79,0xf9,
  0x05,0x85,0x45,0xc5,0x25,0xa5,0x65,0xe5,0x15,0x95,0x55,0xd5,0x35,0xb5,0x75,0xf5,
  0x0d,0x8d,0x4d,0xcd,0x2d,0xad,0x6d,0xed,0x1d,0x9d,0x5d,0xdd,0x3d,0xbd,0x7d,0xfd,
  0x03,0x83,0x43,0xc3,0x23,0xa3,0x63,0xe3,0x13,0x93,0x53,0xd3,0x33,0xb3,0x73,0xf3,
  0x0b,0x8b,0x4b,0xcb,0x2b,0xab,0x6b,0xeb,0x1b,0x9b,0x5b,0xdb,0x3b,0xbb,0x7b,0xfb,
  0x07,0x87,0x47,0xc7,0x27,0xa7,0x67,0xe7,0x17,0x97,0x57,0xd7,0x37,0xb7,0x77,0xf7,
  0x0f,0x8f,0x4f,0xcf,0x2f,0xaf,0x6f,0xef,0x1f,0x9f,0x5f,0xdf,0x3f,0xbf,0x7f,0xff
};

static const double DFTTab[][2] =
{
{ 1.00000000000000000, 0.00000000000000000 },
{-1.00000000000000000, 0.00000000000000000 },
{ 0.00000000000000000, 1.00000000000000000 },
{ 0.70710678118654757, 0.70710678118654746 },
{ 0.92387953251128674, 0.38268343236508978 },
{ 0.98078528040323043, 0.19509032201612825 },
{ 0.99518472667219693, 0.09801714032956060 },
{ 0.99879545620517241, 0.04906767432741802 },
{ 0.99969881869620425, 0.02454122852291229 },
{ 0.99992470183914450, 0.01227153828571993 },
{ 0.99998117528260111, 0.00613588464915448 },
{ 0.99999529380957619, 0.00306795676296598 },
{ 0.99999882345170188, 0.00153398018628477 },
{ 0.99999970586288223, 0.00076699031874270 },
{ 0.99999992646571789, 0.00038349518757140 },
{ 0.99999998161642933, 0.00019174759731070 },
{ 0.99999999540410733, 0.00009587379909598 },
{ 0.99999999885102686, 0.00004793689960307 },
{ 0.99999999971275666, 0.00002396844980842 },
{ 0.99999999992818922, 0.00001198422490507 },
{ 0.99999999998204725, 0.00000599211245264 },
{ 0.99999999999551181, 0.00000299605622633 },
{ 0.99999999999887801, 0.00000149802811317 },
{ 0.99999999999971945, 0.00000074901405658 },
{ 0.99999999999992983, 0.00000037450702829 },
{ 0.99999999999998246, 0.00000018725351415 },
{ 0.99999999999999567, 0.00000009362675707 },
{ 0.99999999999999889, 0.00000004681337854 },
{ 0.99999999999999978, 0.00000002340668927 },
{ 0.99999999999999989, 0.00000001170334463 },
{ 1.00000000000000000, 0.00000000585167232 },
{ 1.00000000000000000, 0.00000000292583616 }
};

#define BitRev(i,shift) \
   ((int)((((unsigned)bitrevTab[(i)&255] << 24)+ \
           ((unsigned)bitrevTab[((i)>> 8)&255] << 16)+ \
           ((unsigned)bitrevTab[((i)>>16)&255] <<  8)+ \
           ((unsigned)bitrevTab[((i)>>24)])) >> (shift)))

static int
DFTFactorize( int n, int* factors )
{
    int nf = 0, f, i, j;

    if( n <= 5 )
    {
        factors[0] = n;
        return 1;
    }

    f = (((n - 1)^n)+1) >> 1;
    if( f > 1 )
    {
        factors[nf++] = f;
        n = f == n ? 1 : n/f;
    }

    for( f = 3; n > 1; )
    {
        int d = n/f;
        if( d*f == n )
        {
            factors[nf++] = f;
            n = d;
        }
        else
        {
            f += 2;
            if( f*f > n )
                break;
        }
    }

    if( n > 1 )
        factors[nf++] = n;

    f = (factors[0] & 1) == 0;
    for( i = f; i < (nf+f)/2; i++ )
        DXT_SWAP( factors[i], factors[nf-i-1+f], j );

    return nf;
}

static void
DFTInit( int n0, int nf, int* factors, int* itab, int elem_size, void* _wave, int inv_itab )
{
    int digits[34], radix[34];
    int n = factors[0], m = 0;
    int* itab0 = itab;
    int i, j, k;
    Complex<double> w, w1;
    double t;

    if( n0 <= 5 )
    {
        itab[0] = 0;
        itab[n0-1] = n0-1;

        if( n0 != 4 )
        {
            for( i = 1; i < n0-1; i++ )
                itab[i] = i;
        }
        else
        {
            itab[1] = 2;
            itab[2] = 1;
        }
        if( n0 == 5 )
        {
            if( elem_size == sizeof(Complex<double>) )
                ((Complex<double>*)_wave)[0] = Complex<double>(1.,0.);
            else
                ((Complex<float>*)_wave)[0] = Complex<float>(1.f,0.f);
        }
        if( n0 != 4 )
            return;
        m = 2;
    }
    else
    {
        // radix[] is initialized from index 'nf' down to zero
        assert (nf < 34);
        radix[nf] = 1;
        digits[nf] = 0;
        for( i = 0; i < nf; i++ )
        {
            digits[i] = 0;
            radix[nf-i-1] = radix[nf-i]*factors[nf-i-1];
        }

        if( inv_itab && factors[0] != factors[nf-1] )
            itab = (int*)_wave;

        if( (n & 1) == 0 )
        {
            int a = radix[1], na2 = n*a>>1, na4 = na2 >> 1;
            for( m = 0; (unsigned)(1 << m) < (unsigned)n; m++ )
                ;
            if( n <= 2  )
            {
                itab[0] = 0;
                itab[1] = na2;
            }
            else if( n <= 256 )
            {
                int shift = 10 - m;
                for( i = 0; i <= n - 4; i += 4 )
                {
                    j = (bitrevTab[i>>2]>>shift)*a;
                    itab[i] = j;
                    itab[i+1] = j + na2;
                    itab[i+2] = j + na4;
                    itab[i+3] = j + na2 + na4;
                }
            }
            else
            {
                int shift = 34 - m;
                for( i = 0; i < n; i += 4 )
                {
                    int i4 = i >> 2;
                    j = BitRev(i4,shift)*a;
                    itab[i] = j;
                    itab[i+1] = j + na2;
                    itab[i+2] = j + na4;
                    itab[i+3] = j + na2 + na4;
                }
            }

            digits[1]++;

            if( nf >= 2 )
            {
                for( i = n, j = radix[2]; i < n0; )
                {
                    for( k = 0; k < n; k++ )
                        itab[i+k] = itab[k] + j;
                    if( (i += n) >= n0 )
                        break;
                    j += radix[2];
                    for( k = 1; ++digits[k] >= factors[k]; k++ )
                    {
                        digits[k] = 0;
                        j += radix[k+2] - radix[k];
                    }
                }
            }
        }
        else
        {
            for( i = 0, j = 0;; )
            {
                itab[i] = j;
                if( ++i >= n0 )
                    break;
                j += radix[1];
                for( k = 0; ++digits[k] >= factors[k]; k++ )
                {
                    digits[k] = 0;
                    j += radix[k+2] - radix[k];
                }
            }
        }

        if( itab != itab0 )
        {
            itab0[0] = 0;
            for( i = n0 & 1; i < n0; i += 2 )
            {
                int k0 = itab[i];
                int k1 = itab[i+1];
                itab0[k0] = i;
                itab0[k1] = i+1;
            }
        }
    }

    if( (n0 & (n0-1)) == 0 )
    {
        w.re = w1.re = DFTTab[m][0];
        w.im = w1.im = -DFTTab[m][1];
    }
    else
    {
        t = -DXT_PI*2/n0;
        w.im = w1.im = sin(t);
        w.re = w1.re = std::sqrt(1. - w1.im*w1.im);	// remove std::
    }
    n = (n0+1)/2;

    if( elem_size == sizeof(Complex<double>) )
    {
        Complex<double>* wave = (Complex<double>*)_wave;

        wave[0].re = 1.;
        wave[0].im = 0.;

        if( (n0 & 1) == 0 )
        {
            wave[n].re = -1.;
            wave[n].im = 0;
        }

        for( i = 1; i < n; i++ )
        {
            wave[i] = w;
            wave[n0-i].re = w.re;
            wave[n0-i].im = -w.im;

            t = w.re*w1.re - w.im*w1.im;
            w.im = w.re*w1.im + w.im*w1.re;
            w.re = t;
        }
    }
    else
    {
        Complex<float>* wave = (Complex<float>*)_wave;
        assert( elem_size == sizeof(Complex<float>) );

        wave[0].re = 1.f;
        wave[0].im = 0.f;

        if( (n0 & 1) == 0 )
        {
            wave[n].re = -1.f;
            wave[n].im = 0.f;
        }

        for( i = 1; i < n; i++ )
        {
            wave[i].re = (float)w.re;
            wave[i].im = (float)w.im;
            wave[n0-i].re = (float)w.re;
            wave[n0-i].im = (float)-w.im;

            t = w.re*w1.re - w.im*w1.im;
            w.im = w.re*w1.im + w.im*w1.re;
            w.re = t;
        }
    }
}

template<typename T> struct DFT_VecR4
{
    int operator()(Complex<T>*, int, int, int&, const Complex<T>*) const { return 1; }
};

enum { DFT_NO_PERMUTE=256, DFT_COMPLEX_INPUT_OR_OUTPUT=512 };


// mixed-radix complex discrete Fourier transform: double-precision version
template<typename T> static void
DFT( const Complex<T>* src, Complex<T>* dst, int n,
     int nf, const int* factors, const int* itab,
     const Complex<T>* wave, int tab_size,
     const void*
#ifdef USE_IPP_DFT
     spec
#endif
     , Complex<T>* buf,
     int flags, double _scale )
{
    static const T sin_120 = (T)0.86602540378443864676372317075294;
    static const T fft5_2 = (T)0.559016994374947424102293417182819;
    static const T fft5_3 = (T)-0.951056516295153572116439333379382;
    static const T fft5_4 = (T)-1.538841768587626701285145288018455;
    static const T fft5_5 = (T)0.363271264002680442947733378740309;

    int n0 = n, f_idx, nx;
    int inv = flags & DFT_INVERSE;
    int dw0 = tab_size, dw;
    int i, j, k;
    Complex<T> t;
    T scale = (T)_scale;
    int tab_step;

    tab_step = tab_size == n ? 1 : tab_size == n*2 ? 2 : tab_size/n;

    // 0. shuffle data
    if( dst != src )
    {
        assert( (flags & DFT_NO_PERMUTE) == 0 );
        if( !inv )
        {
            for( i = 0; i <= n - 2; i += 2, itab += 2*tab_step )
            {
                int k0 = itab[0], k1 = itab[tab_step];
                assert( (unsigned)k0 < (unsigned)n && (unsigned)k1 < (unsigned)n );
                dst[i] = src[k0]; dst[i+1] = src[k1];
            }

            if( i < n )
                dst[n-1] = src[n-1];
        }
        else
        {
            for( i = 0; i <= n - 2; i += 2, itab += 2*tab_step )
            {
                int k0 = itab[0], k1 = itab[tab_step];
                assert( (unsigned)k0 < (unsigned)n && (unsigned)k1 < (unsigned)n );
                t.re = src[k0].re; t.im = -src[k0].im;
                dst[i] = t;
                t.re = src[k1].re; t.im = -src[k1].im;
                dst[i+1] = t;
            }

            if( i < n )
            {
                t.re = src[n-1].re; t.im = -src[n-1].im;
                dst[i] = t;
            }
        }
    }
    else
    {
        if( (flags & DFT_NO_PERMUTE) == 0 )
        {
            assert( factors[0] == factors[nf-1] );
            if( nf == 1 )
            {
                if( (n & 3) == 0 )
                {
                    int n2 = n/2;
                    Complex<T>* dsth = dst + n2;

                    for( i = 0; i < n2; i += 2, itab += tab_step*2 )
                    {
                        j = itab[0];
                        assert( (unsigned)j < (unsigned)n2 );

                        DXT_SWAP(dst[i+1], dsth[j], t);
                        if( j > i )
                        {
                            DXT_SWAP(dst[i], dst[j], t);
                            DXT_SWAP(dsth[i+1], dsth[j+1], t);
                        }
                    }
                }
                // else do nothing
            }
            else
            {
                for( i = 0; i < n; i++, itab += tab_step )
                {
                    j = itab[0];
                    assert( (unsigned)j < (unsigned)n );
                    if( j > i )
                        DXT_SWAP(dst[i], dst[j], t);
                }
            }
        }

        if( inv )
        {
            for( i = 0; i <= n - 2; i += 2 )
            {
                T t0 = -dst[i].im;
                T t1 = -dst[i+1].im;
                dst[i].im = t0; dst[i+1].im = t1;
            }

            if( i < n )
                dst[n-1].im = -dst[n-1].im;
        }
    }

    n = 1;
    // 1. power-2 transforms
    if( (factors[0] & 1) == 0 )
    {
        //if( 0 && factors[0] >= 4 && checkHardwareSupport(CV_CPU_SSE3) )//no support
        //{
        //    DFT_VecR4<T> vr4;
        //    n = vr4(dst, factors[0], n0, dw0, wave);
        //}

        // radix-4 transform
        for( ; n*4 <= factors[0]; )
        {
            nx = n;
            n *= 4;
            dw0 /= 4;

            for( i = 0; i < n0; i += n )
            {
                Complex<T> *v0, *v1;
                T r0, i0, r1, i1, r2, i2, r3, i3, r4, i4;

                v0 = dst + i;
                v1 = v0 + nx*2;

                r0 = v1[0].re; i0 = v1[0].im;
                r4 = v1[nx].re; i4 = v1[nx].im;

                r1 = r0 + r4; i1 = i0 + i4;
                r3 = i0 - i4; i3 = r4 - r0;

                r2 = v0[0].re; i2 = v0[0].im;
                r4 = v0[nx].re; i4 = v0[nx].im;

                r0 = r2 + r4; i0 = i2 + i4;
                r2 -= r4; i2 -= i4;

                v0[0].re = r0 + r1; v0[0].im = i0 + i1;
                v1[0].re = r0 - r1; v1[0].im = i0 - i1;
                v0[nx].re = r2 + r3; v0[nx].im = i2 + i3;
                v1[nx].re = r2 - r3; v1[nx].im = i2 - i3;

                for( j = 1, dw = dw0; j < nx; j++, dw += dw0 )
                {
                    v0 = dst + i + j;
                    v1 = v0 + nx*2;

                    r2 = v0[nx].re*wave[dw*2].re - v0[nx].im*wave[dw*2].im;
                    i2 = v0[nx].re*wave[dw*2].im + v0[nx].im*wave[dw*2].re;
                    r0 = v1[0].re*wave[dw].im + v1[0].im*wave[dw].re;
                    i0 = v1[0].re*wave[dw].re - v1[0].im*wave[dw].im;
                    r3 = v1[nx].re*wave[dw*3].im + v1[nx].im*wave[dw*3].re;
                    i3 = v1[nx].re*wave[dw*3].re - v1[nx].im*wave[dw*3].im;

                    r1 = i0 + i3; i1 = r0 + r3;
                    r3 = r0 - r3; i3 = i3 - i0;
                    r4 = v0[0].re; i4 = v0[0].im;

                    r0 = r4 + r2; i0 = i4 + i2;
                    r2 = r4 - r2; i2 = i4 - i2;

                    v0[0].re = r0 + r1; v0[0].im = i0 + i1;
                    v1[0].re = r0 - r1; v1[0].im = i0 - i1;
                    v0[nx].re = r2 + r3; v0[nx].im = i2 + i3;
                    v1[nx].re = r2 - r3; v1[nx].im = i2 - i3;
                }
            }
        }

        for( ; n < factors[0]; )
        {
            // do the remaining radix-2 transform
            nx = n;
            n *= 2;
            dw0 /= 2;

            for( i = 0; i < n0; i += n )
            {
                Complex<T>* v = dst + i;
                T r0 = v[0].re + v[nx].re;
                T i0 = v[0].im + v[nx].im;
                T r1 = v[0].re - v[nx].re;
                T i1 = v[0].im - v[nx].im;
                v[0].re = r0; v[0].im = i0;
                v[nx].re = r1; v[nx].im = i1;

                for( j = 1, dw = dw0; j < nx; j++, dw += dw0 )
                {
                    v = dst + i + j;
                    r1 = v[nx].re*wave[dw].re - v[nx].im*wave[dw].im;
                    i1 = v[nx].im*wave[dw].re + v[nx].re*wave[dw].im;
                    r0 = v[0].re; i0 = v[0].im;

                    v[0].re = r0 + r1; v[0].im = i0 + i1;
                    v[nx].re = r0 - r1; v[nx].im = i0 - i1;
                }
            }
        }
    }

    // 2. all the other transforms
    for( f_idx = (factors[0]&1) ? 0 : 1; f_idx < nf; f_idx++ )
    {
        int factor = factors[f_idx];
        nx = n;
        n *= factor;
        dw0 /= factor;

        if( factor == 3 )
        {
            // radix-3
            for( i = 0; i < n0; i += n )
            {
                Complex<T>* v = dst + i;

                T r1 = v[nx].re + v[nx*2].re;
                T i1 = v[nx].im + v[nx*2].im;
                T r0 = v[0].re;
                T i0 = v[0].im;
                T r2 = sin_120*(v[nx].im - v[nx*2].im);
                T i2 = sin_120*(v[nx*2].re - v[nx].re);
                v[0].re = r0 + r1; v[0].im = i0 + i1;
                r0 -= (T)0.5*r1; i0 -= (T)0.5*i1;
                v[nx].re = r0 + r2; v[nx].im = i0 + i2;
                v[nx*2].re = r0 - r2; v[nx*2].im = i0 - i2;

                for( j = 1, dw = dw0; j < nx; j++, dw += dw0 )
                {
                    v = dst + i + j;
                    r0 = v[nx].re*wave[dw].re - v[nx].im*wave[dw].im;
                    i0 = v[nx].re*wave[dw].im + v[nx].im*wave[dw].re;
                    i2 = v[nx*2].re*wave[dw*2].re - v[nx*2].im*wave[dw*2].im;
                    r2 = v[nx*2].re*wave[dw*2].im + v[nx*2].im*wave[dw*2].re;
                    r1 = r0 + i2; i1 = i0 + r2;

                    r2 = sin_120*(i0 - r2); i2 = sin_120*(i2 - r0);
                    r0 = v[0].re; i0 = v[0].im;
                    v[0].re = r0 + r1; v[0].im = i0 + i1;
                    r0 -= (T)0.5*r1; i0 -= (T)0.5*i1;
                    v[nx].re = r0 + r2; v[nx].im = i0 + i2;
                    v[nx*2].re = r0 - r2; v[nx*2].im = i0 - i2;
                }
            }
        }
        else if( factor == 5 )
        {
            // radix-5
            for( i = 0; i < n0; i += n )
            {
                for( j = 0, dw = 0; j < nx; j++, dw += dw0 )
                {
                    Complex<T>* v0 = dst + i + j;
                    Complex<T>* v1 = v0 + nx*2;
                    Complex<T>* v2 = v1 + nx*2;

                    T r0, i0, r1, i1, r2, i2, r3, i3, r4, i4, r5, i5;

                    r3 = v0[nx].re*wave[dw].re - v0[nx].im*wave[dw].im;
                    i3 = v0[nx].re*wave[dw].im + v0[nx].im*wave[dw].re;
                    r2 = v2[0].re*wave[dw*4].re - v2[0].im*wave[dw*4].im;
                    i2 = v2[0].re*wave[dw*4].im + v2[0].im*wave[dw*4].re;

                    r1 = r3 + r2; i1 = i3 + i2;
                    r3 -= r2; i3 -= i2;

                    r4 = v1[nx].re*wave[dw*3].re - v1[nx].im*wave[dw*3].im;
                    i4 = v1[nx].re*wave[dw*3].im + v1[nx].im*wave[dw*3].re;
                    r0 = v1[0].re*wave[dw*2].re - v1[0].im*wave[dw*2].im;
                    i0 = v1[0].re*wave[dw*2].im + v1[0].im*wave[dw*2].re;

                    r2 = r4 + r0; i2 = i4 + i0;
                    r4 -= r0; i4 -= i0;

                    r0 = v0[0].re; i0 = v0[0].im;
                    r5 = r1 + r2; i5 = i1 + i2;

                    v0[0].re = r0 + r5; v0[0].im = i0 + i5;

                    r0 -= (T)0.25*r5; i0 -= (T)0.25*i5;
                    r1 = fft5_2*(r1 - r2); i1 = fft5_2*(i1 - i2);
                    r2 = -fft5_3*(i3 + i4); i2 = fft5_3*(r3 + r4);

                    i3 *= -fft5_5; r3 *= fft5_5;
                    i4 *= -fft5_4; r4 *= fft5_4;

                    r5 = r2 + i3; i5 = i2 + r3;
                    r2 -= i4; i2 -= r4;

                    r3 = r0 + r1; i3 = i0 + i1;
                    r0 -= r1; i0 -= i1;

                    v0[nx].re = r3 + r2; v0[nx].im = i3 + i2;
                    v2[0].re = r3 - r2; v2[0].im = i3 - i2;

                    v1[0].re = r0 + r5; v1[0].im = i0 + i5;
                    v1[nx].re = r0 - r5; v1[nx].im = i0 - i5;
                }
            }
        }
        else
        {
            // radix-"factor" - an odd number
            int p, q, factor2 = (factor - 1)/2;
            int d, dd, dw_f = tab_size/factor;
            Complex<T>* a = buf;
            Complex<T>* b = buf + factor2;

            for( i = 0; i < n0; i += n )
            {
                for( j = 0, dw = 0; j < nx; j++, dw += dw0 )
                {
                    Complex<T>* v = dst + i + j;
                    Complex<T> v_0 = v[0];
                    Complex<T> vn_0 = v_0;

                    if( j == 0 )
                    {
                        for( p = 1, k = nx; p <= factor2; p++, k += nx )
                        {
                            T r0 = v[k].re + v[n-k].re;
                            T i0 = v[k].im - v[n-k].im;
                            T r1 = v[k].re - v[n-k].re;
                            T i1 = v[k].im + v[n-k].im;

                            vn_0.re += r0; vn_0.im += i1;
                            a[p-1].re = r0; a[p-1].im = i0;
                            b[p-1].re = r1; b[p-1].im = i1;
                        }
                    }
                    else
                    {
                        const Complex<T>* wave_ = wave + dw*factor;
                        d = dw;

                        for( p = 1, k = nx; p <= factor2; p++, k += nx, d += dw )
                        {
                            T r2 = v[k].re*wave[d].re - v[k].im*wave[d].im;
                            T i2 = v[k].re*wave[d].im + v[k].im*wave[d].re;

                            T r1 = v[n-k].re*wave_[-d].re - v[n-k].im*wave_[-d].im;
                            T i1 = v[n-k].re*wave_[-d].im + v[n-k].im*wave_[-d].re;

                            T r0 = r2 + r1;
                            T i0 = i2 - i1;
                            r1 = r2 - r1;
                            i1 = i2 + i1;

                            vn_0.re += r0; vn_0.im += i1;
                            a[p-1].re = r0; a[p-1].im = i0;
                            b[p-1].re = r1; b[p-1].im = i1;
                        }
                    }

                    v[0] = vn_0;

                    for( p = 1, k = nx; p <= factor2; p++, k += nx )
                    {
                        Complex<T> s0 = v_0, s1 = v_0;
                        d = dd = dw_f*p;

                        for( q = 0; q < factor2; q++ )
                        {
                            T r0 = wave[d].re * a[q].re;
                            T i0 = wave[d].im * a[q].im;
                            T r1 = wave[d].re * b[q].im;
                            T i1 = wave[d].im * b[q].re;

                            s1.re += r0 + i0; s0.re += r0 - i0;
                            s1.im += r1 - i1; s0.im += r1 + i1;

                            d += dd;
                            d -= -(d >= tab_size) & tab_size;
                        }

                        v[k] = s0;
                        v[n-k] = s1;
                    }
                }
            }
        }
    }

    if( scale != 1 )
    {
        T re_scale = scale, im_scale = scale;
        if( inv )
            im_scale = -im_scale;

        for( i = 0; i < n0; i++ )
        {
            T t0 = dst[i].re*re_scale;
            T t1 = dst[i].im*im_scale;
            dst[i].re = t0;
            dst[i].im = t1;
        }
    }
    else if( inv )
    {
        for( i = 0; i <= n0 - 2; i += 2 )
        {
            T t0 = -dst[i].im;
            T t1 = -dst[i+1].im;
            dst[i].im = t0;
            dst[i+1].im = t1;
        }

        if( i < n0 )
            dst[n0-1].im = -dst[n0-1].im;
    }
}


/* FFT of real vector
   output vector format:
     re(0), re(1), im(1), ... , re(n/2-1), im((n+1)/2-1) [, re((n+1)/2)] OR ...
     re(0), 0, re(1), im(1), ..., re(n/2-1), im((n+1)/2-1) [, re((n+1)/2), 0] */
template<typename T> static void
RealDFT( const T* src, T* dst, int n, int nf, int* factors, const int* itab,
         const Complex<T>* wave, int tab_size, const void*
#ifdef USE_IPP_DFT
         spec
#endif
         ,
         Complex<T>* buf, int flags, double _scale )
{
    int complex_output = (flags & DFT_COMPLEX_INPUT_OR_OUTPUT) != 0;
    T scale = (T)_scale;
    int j, n2 = n >> 1;
    dst += complex_output;

    assert( tab_size == n );

    if( n == 1 )
    {
        dst[0] = src[0]*scale;
    }
    else if( n == 2 )
    {
        T t = (src[0] + src[1])*scale;
        dst[1] = (src[0] - src[1])*scale;
        dst[0] = t;
    }
    else if( n & 1 )
    {
        dst -= complex_output;
        Complex<T>* _dst = (Complex<T>*)dst;
        _dst[0].re = src[0]*scale;
        _dst[0].im = 0;
        for( j = 1; j < n; j += 2 )
        {
            T t0 = src[itab[j]]*scale;
            T t1 = src[itab[j+1]]*scale;
            _dst[j].re = t0;
            _dst[j].im = 0;
            _dst[j+1].re = t1;
            _dst[j+1].im = 0;
        }
        DFT( _dst, _dst, n, nf, factors, itab, wave,
             tab_size, 0, buf, DFT_NO_PERMUTE, 1 );
        if( !complex_output )
            dst[1] = dst[0];
    }
    else
    {
        T t0, t;
        T h1_re, h1_im, h2_re, h2_im;
        T scale2 = scale*(T)0.5;
        factors[0] >>= 1;

        DFT( (Complex<T>*)src, (Complex<T>*)dst, n2, nf - (factors[0] == 1),
             factors + (factors[0] == 1),
             itab, wave, tab_size, 0, buf, 0, 1 );
        factors[0] <<= 1;

        t = dst[0] - dst[1];
        dst[0] = (dst[0] + dst[1])*scale;
        dst[1] = t*scale;

        t0 = dst[n2];
        t = dst[n-1];
        dst[n-1] = dst[1];

        for( j = 2, wave++; j < n2; j += 2, wave++ )
        {
            /* calc odd */
            h2_re = scale2*(dst[j+1] + t);
            h2_im = scale2*(dst[n-j] - dst[j]);

            /* calc even */
            h1_re = scale2*(dst[j] + dst[n-j]);
            h1_im = scale2*(dst[j+1] - t);

            /* rotate */
            t = h2_re*wave->re - h2_im*wave->im;
            h2_im = h2_re*wave->im + h2_im*wave->re;
            h2_re = t;
            t = dst[n-j-1];

            dst[j-1] = h1_re + h2_re;
            dst[n-j-1] = h1_re - h2_re;
            dst[j] = h1_im + h2_im;
            dst[n-j] = h2_im - h1_im;
        }

        if( j <= n2 )
        {
            dst[n2-1] = t0*scale;
            dst[n2] = -t*scale;
        }
    }

    if( complex_output && (n & 1) == 0 )
    {
        dst[-1] = dst[0];
        dst[0] = 0;
        dst[n] = 0;
    }
}

/* Inverse FFT of complex conjugate-symmetric vector
   input vector format:
      re[0], re[1], im[1], ... , re[n/2-1], im[n/2-1], re[n/2] OR
      re(0), 0, re(1), im(1), ..., re(n/2-1), im((n+1)/2-1) [, re((n+1)/2), 0] */
template<typename T> static void
CCSIDFT( const T* src, T* dst, int n, int nf, int* factors, const int* itab,
         const Complex<T>* wave, int tab_size,
         const void*
#ifdef USE_IPP_DFT
         spec
#endif
         , Complex<T>* buf,
         int flags, double _scale )
{
    int complex_input = (flags & DFT_COMPLEX_INPUT_OR_OUTPUT) != 0;
    int j, k, n2 = (n+1) >> 1;
    T scale = (T)_scale;
    T save_s1 = 0.;
    T t0, t1, t2, t3, t;

    assert( tab_size == n );

    if( complex_input )
    {
        assert( src != dst );
        save_s1 = src[1];
        ((T*)src)[1] = src[0];
        src++;
    }

    if( n == 1 )
    {
        dst[0] = (T)(src[0]*scale);
    }
    else if( n == 2 )
    {
        t = (src[0] + src[1])*scale;
        dst[1] = (src[0] - src[1])*scale;
        dst[0] = t;
    }
    else if( n & 1 )
    {
        Complex<T>* _src = (Complex<T>*)(src-1);
        Complex<T>* _dst = (Complex<T>*)dst;

        _dst[0].re = src[0];
        _dst[0].im = 0;
        for( j = 1; j < n2; j++ )
        {
            int k0 = itab[j], k1 = itab[n-j];
            t0 = _src[j].re; t1 = _src[j].im;
            _dst[k0].re = t0; _dst[k0].im = -t1;
            _dst[k1].re = t0; _dst[k1].im = t1;
        }

        DFT( _dst, _dst, n, nf, factors, itab, wave,
             tab_size, 0, buf, DFT_NO_PERMUTE, 1. );
        dst[0] *= scale;
        for( j = 1; j < n; j += 2 )
        {
            t0 = dst[j*2]*scale;
            t1 = dst[j*2+2]*scale;
            dst[j] = t0;
            dst[j+1] = t1;
        }
    }
    else
    {
        int inplace = src == dst;
        const Complex<T>* w = wave;

        t = src[1];
        t0 = (src[0] + src[n-1]);
        t1 = (src[n-1] - src[0]);
        dst[0] = t0;
        dst[1] = t1;

        for( j = 2, w++; j < n2; j += 2, w++ )
        {
            T h1_re, h1_im, h2_re, h2_im;

            h1_re = (t + src[n-j-1]);
            h1_im = (src[j] - src[n-j]);

            h2_re = (t - src[n-j-1]);
            h2_im = (src[j] + src[n-j]);

            t = h2_re*w->re + h2_im*w->im;
            h2_im = h2_im*w->re - h2_re*w->im;
            h2_re = t;

            t = src[j+1];
            t0 = h1_re - h2_im;
            t1 = -h1_im - h2_re;
            t2 = h1_re + h2_im;
            t3 = h1_im - h2_re;

            if( inplace )
            {
                dst[j] = t0;
                dst[j+1] = t1;
                dst[n-j] = t2;
                dst[n-j+1]= t3;
            }
            else
            {
                int j2 = j >> 1;
                k = itab[j2];
                dst[k] = t0;
                dst[k+1] = t1;
                k = itab[n2-j2];
                dst[k] = t2;
                dst[k+1]= t3;
            }
        }

        if( j <= n2 )
        {
            t0 = t*2;
            t1 = src[n2]*2;

            if( inplace )
            {
                dst[n2] = t0;
                dst[n2+1] = t1;
            }
            else
            {
                k = itab[n2];
                dst[k*2] = t0;
                dst[k*2+1] = t1;
            }
        }

        factors[0] >>= 1;
        DFT( (Complex<T>*)dst, (Complex<T>*)dst, n2,
             nf - (factors[0] == 1),
             factors + (factors[0] == 1), itab,
             wave, tab_size, 0, buf,
             inplace ? 0 : DFT_NO_PERMUTE, 1. );
        factors[0] <<= 1;

        for( j = 0; j < n; j += 2 )
        {
            t0 = dst[j]*scale;
            t1 = dst[j+1]*(-scale);
            dst[j] = t0;
            dst[j+1] = t1;
        }
    }

#ifdef USE_IPP_DFT
finalize:
#endif
    if( complex_input )
        ((T*)src)[0] = (T)save_s1;
}

static void
CopyColumn( const uchar* _src, size_t src_step,
            uchar* _dst, size_t dst_step,
            int len, size_t elem_size )
{
    int i, t0, t1;
    const int* src = (const int*)_src;
    int* dst = (int*)_dst;
    src_step /= sizeof(src[0]);
    dst_step /= sizeof(dst[0]);

    if( elem_size == sizeof(int) )
    {
        for( i = 0; i < len; i++, src += src_step, dst += dst_step )
            dst[0] = src[0];
    }
    else if( elem_size == sizeof(int)*2 )
    {
        for( i = 0; i < len; i++, src += src_step, dst += dst_step )
        {
            t0 = src[0]; t1 = src[1];
            dst[0] = t0; dst[1] = t1;
        }
    }
    else if( elem_size == sizeof(int)*4 )
    {
        for( i = 0; i < len; i++, src += src_step, dst += dst_step )
        {
            t0 = src[0]; t1 = src[1];
            dst[0] = t0; dst[1] = t1;
            t0 = src[2]; t1 = src[3];
            dst[2] = t0; dst[3] = t1;
        }
    }
}


static void
CopyFrom2Columns( const uchar* _src, size_t src_step,
                  uchar* _dst0, uchar* _dst1,
                  int len, size_t elem_size )
{
    int i, t0, t1;
    const int* src = (const int*)_src;
    int* dst0 = (int*)_dst0;
    int* dst1 = (int*)_dst1;
    src_step /= sizeof(src[0]);

    if( elem_size == sizeof(int) )
    {
        for( i = 0; i < len; i++, src += src_step )
        {
            t0 = src[0]; t1 = src[1];
            dst0[i] = t0; dst1[i] = t1;
        }
    }
    else if( elem_size == sizeof(int)*2 )
    {
        for( i = 0; i < len*2; i += 2, src += src_step )
        {
            t0 = src[0]; t1 = src[1];
            dst0[i] = t0; dst0[i+1] = t1;
            t0 = src[2]; t1 = src[3];
            dst1[i] = t0; dst1[i+1] = t1;
        }
    }
    else if( elem_size == sizeof(int)*4 )
    {
        for( i = 0; i < len*4; i += 4, src += src_step )
        {
            t0 = src[0]; t1 = src[1];
            dst0[i] = t0; dst0[i+1] = t1;
            t0 = src[2]; t1 = src[3];
            dst0[i+2] = t0; dst0[i+3] = t1;
            t0 = src[4]; t1 = src[5];
            dst1[i] = t0; dst1[i+1] = t1;
            t0 = src[6]; t1 = src[7];
            dst1[i+2] = t0; dst1[i+3] = t1;
        }
    }
}


static void
CopyTo2Columns( const uchar* _src0, const uchar* _src1,
                uchar* _dst, size_t dst_step,
                int len, size_t elem_size )
{
    int i, t0, t1;
    const int* src0 = (const int*)_src0;
    const int* src1 = (const int*)_src1;
    int* dst = (int*)_dst;
    dst_step /= sizeof(dst[0]);

    if( elem_size == sizeof(int) )
    {
        for( i = 0; i < len; i++, dst += dst_step )
        {
            t0 = src0[i]; t1 = src1[i];
            dst[0] = t0; dst[1] = t1;
        }
    }
    else if( elem_size == sizeof(int)*2 )
    {
        for( i = 0; i < len*2; i += 2, dst += dst_step )
        {
            t0 = src0[i]; t1 = src0[i+1];
            dst[0] = t0; dst[1] = t1;
            t0 = src1[i]; t1 = src1[i+1];
            dst[2] = t0; dst[3] = t1;
        }
    }
    else if( elem_size == sizeof(int)*4 )
    {
        for( i = 0; i < len*4; i += 4, dst += dst_step )
        {
            t0 = src0[i]; t1 = src0[i+1];
            dst[0] = t0; dst[1] = t1;
            t0 = src0[i+2]; t1 = src0[i+3];
            dst[2] = t0; dst[3] = t1;
            t0 = src1[i]; t1 = src1[i+1];
            dst[4] = t0; dst[5] = t1;
            t0 = src1[i+2]; t1 = src1[i+3];
            dst[6] = t0; dst[7] = t1;
        }
    }
}


static void
ExpandCCS( uchar* _ptr, int n, int elem_size )
{
    int i;
    if( elem_size == (int)sizeof(float) )
    {
        float* p = (float*)_ptr;
        for( i = 1; i < (n+1)/2; i++ )
        {
            p[(n-i)*2] = p[i*2-1];
            p[(n-i)*2+1] = -p[i*2];
        }
        if( (n & 1) == 0 )
        {
            p[n] = p[n-1];
            p[n+1] = 0.f;
            n--;
        }
        for( i = n-1; i > 0; i-- )
            p[i+1] = p[i];
        p[1] = 0.f;
    }
    else
    {
        double* p = (double*)_ptr;
        for( i = 1; i < (n+1)/2; i++ )
        {
            p[(n-i)*2] = p[i*2-1];
            p[(n-i)*2+1] = -p[i*2];
        }
        if( (n & 1) == 0 )
        {
            p[n] = p[n-1];
            p[n+1] = 0.f;
            n--;
        }
        for( i = n-1; i > 0; i-- )
            p[i+1] = p[i];
        p[1] = 0.f;
    }
}


typedef void (*DFTFunc)(
     const void* src, void* dst, int n, int nf, int* factors,
     const int* itab, const void* wave, int tab_size,
     const void* spec, void* buf, int inv, double scale );

static void DFT_32f( const Complexf* src, Complexf* dst, int n,
    int nf, const int* factors, const int* itab,
    const Complexf* wave, int tab_size,
    const void* spec, Complexf* buf,
    int flags, double scale )
{
    DFT(src, dst, n, nf, factors, itab, wave, tab_size, spec, buf, flags, scale);
}

static void DFT_64f( const Complexd* src, Complexd* dst, int n,
    int nf, const int* factors, const int* itab,
    const Complexd* wave, int tab_size,
    const void* spec, Complexd* buf,
    int flags, double scale )
{
    DFT(src, dst, n, nf, factors, itab, wave, tab_size, spec, buf, flags, scale);
}


static void RealDFT_32f( const float* src, float* dst, int n, int nf, int* factors,
        const int* itab,  const Complexf* wave, int tab_size, const void* spec,
        Complexf* buf, int flags, double scale )
{
    RealDFT( src, dst, n, nf, factors, itab, wave, tab_size, spec, buf, flags, scale);
}

static void RealDFT_64f( const double* src, double* dst, int n, int nf, int* factors,
        const int* itab,  const Complexd* wave, int tab_size, const void* spec,
        Complexd* buf, int flags, double scale )
{
    RealDFT( src, dst, n, nf, factors, itab, wave, tab_size, spec, buf, flags, scale);
}

static void CCSIDFT_32f( const float* src, float* dst, int n, int nf, int* factors,
                         const int* itab,  const Complexf* wave, int tab_size, const void* spec,
                         Complexf* buf, int flags, double scale )
{
    CCSIDFT( src, dst, n, nf, factors, itab, wave, tab_size, spec, buf, flags, scale);
}

static void CCSIDFT_64f( const double* src, double* dst, int n, int nf, int* factors,
                         const int* itab,  const Complexd* wave, int tab_size, const void* spec,
                         Complexd* buf, int flags, double scale )
{
    CCSIDFT( src, dst, n, nf, factors, itab, wave, tab_size, spec, buf, flags, scale);
}

}

//#ifdef USE_IPP_DFT
//typedef IppStatus (CV_STDCALL* IppDFTGetSizeFunc)(int, int, IppHintAlgorithm, int*, int*, int*);
//typedef IppStatus (CV_STDCALL* IppDFTInitFunc)(int, int, IppHintAlgorithm, void*, uchar*);
//#endif

void dxt::dft(CMat _src0, CMat& _dst, int flags, int nonzero_rows)
{
    static DFTFunc dft_tbl[6] =
    {
        (DFTFunc)DFT_32f,
        (DFTFunc)RealDFT_32f,
        (DFTFunc)CCSIDFT_32f,
        (DFTFunc)DFT_64f,
        (DFTFunc)RealDFT_64f,
        (DFTFunc)CCSIDFT_64f
    };

    //AutoBuffer<uchar> buf;
	uchar *buff = NULL;
    void *spec = 0;

    //Mat src0 = _src0.getMat(), src = src0;
	CMat src = _src0;
    int prev_len = 0, stage = 0;
    bool inv = (flags & DFT_INVERSE) != 0;
    int nf = 0, real_transform = src.channels() == 1 || (inv && (flags & DFT_REAL_OUTPUT)!=0);
    int type = src.type(), depth = src.depth();
    //int elem_size = (int)src.elemSize1(), complex_elem_size = elem_size*2;
	int elem_size = sizeof(float), complex_elem_size = elem_size*2;
    int factors[34];
    bool inplace_transform = false;
//#ifdef USE_IPP_DFT
//    AutoBuffer<uchar> ippbuf;
//    int ipp_norm_flag = !(flags & DFT_SCALE) ? 8 : inv ? 2 : 1;
//#endif

    //CV_Assert( type == CV_32FC1 || type == CV_32FC2 || type == CV_64FC1 || type == CV_64FC2 );
	assert(type == MAT_32FC2);

    //if( !inv && src.channels() == 1 && (flags & DFT_COMPLEX_OUTPUT) )
    //    _dst.create( src.size(), CV_MAKETYPE(depth, 2) );
    //else if( inv && src.channels() == 2 && (flags & DFT_REAL_OUTPUT) )
    //    _dst.create( src.size(), depth );
    //else
    //    _dst.create( src.size(), type );
	_dst.create(src.size(), MAT_32FC2);

    //Mat dst = _dst.getMat();
	CMat dst = _dst;

    if( !real_transform )
        elem_size = complex_elem_size;

	if (src.cols == 1 && nonzero_rows > 0)
		assert(0);
        //CV_Error( CV_StsNotImplemented,
        //"This mode (using nonzero_rows with a single-column matrix) breaks the function's logic, so it is prohibited.\n"
        //"For fast convolution/correlation use 2-column matrix or single-row matrix instead" );

    // determine, which transform to do first - row-wise
    // (stage 0) or column-wise (stage 1) transform
    /*if( !(flags & DFT_ROWS) && src.rows > 1 &&
        ((src.cols == 1 && (!src.isContinuous() || !dst.isContinuous())) ||
         (src.cols > 1 && inv && real_transform)) )
        stage = 1;*/

    for(;;)
    {
        double scale = 1;
        uchar* wave = 0;
        int* itab = 0;
        uchar* ptr;
        int i, len, count, sz = 0;
        int use_buf = 0, odd_real = 0;
        DFTFunc dft_func;

        if( stage == 0 ) // row-wise transform
        {
            len = !inv ? src.cols : dst.cols;
            count = src.rows;
            if( len == 1 && !(flags & DFT_ROWS) )
            {
                len = !inv ? src.rows : dst.rows;
                count = 1;
            }
            odd_real = real_transform && (len & 1);
        }
        else
        {
            len = dst.rows;
            count = !inv ? src.cols : dst.cols;
            sz = 2*len*complex_elem_size;
        }

        spec = 0;
        {
            if( len != prev_len )
                nf = DFTFactorize( len, factors );

            inplace_transform = factors[0] == factors[nf-1];
            sz += len*(complex_elem_size + sizeof(int));
            i = nf > 1 && (factors[0] & 1) == 0;
            if( (factors[i] & 1) != 0 && factors[i] > 5 )
                sz += (factors[i]+1)*complex_elem_size;

            if( (stage == 0 && ((src.data == dst.data && !inplace_transform) || odd_real)) ||
                (stage == 1 && !inplace_transform) )
            {
                use_buf = 1;
                sz += len*complex_elem_size;
            }
        }

        //ptr = (uchar*)buf;
        //buf.allocate( sz + 32 );
        //if( ptr != (uchar*)buf )
        //    prev_len = 0; // because we release the buffer,
        //                  // force recalculation of
        //                  // twiddle factors and permutation table
		buff = (uchar*)malloc(sz + 32);
        ptr = (uchar*)buff;
        if( !spec )
        {
            wave = ptr;
            ptr += len*complex_elem_size;
            itab = (int*)ptr;
            ptr = (uchar*)dxtAlignPtr( ptr + len*sizeof(int), 16 );

            if( len != prev_len || (!inplace_transform && inv && real_transform))
                DFTInit( len, nf, factors, itab, complex_elem_size,
                            wave, stage == 0 && inv && real_transform );
            // otherwise reuse the tables calculated on the previous stage
        }

        if( stage == 0 )
        {
            uchar* tmp_buf = 0;
            int dptr_offset = 0;
            int dst_full_len = len*elem_size;
            int _flags = (int)inv + (src.channels() != dst.channels() ?
                         DFT_COMPLEX_INPUT_OR_OUTPUT : 0);
            if( use_buf )
            {
                tmp_buf = ptr;
                ptr += len*complex_elem_size;
                if( odd_real && !inv && len > 1 &&
                    !(_flags & DFT_COMPLEX_INPUT_OR_OUTPUT))
                    dptr_offset = elem_size;
            }

            if( !inv && (_flags & DFT_COMPLEX_INPUT_OR_OUTPUT) )
                dst_full_len += (len & 1) ? elem_size : complex_elem_size;

            //dft_func = dft_tbl[(!real_transform ? 0 : !inv ? 1 : 2) + (depth == CV_64F)*3];
			dft_func = dft_tbl[0];

            if( count > 1 && !(flags & DFT_ROWS) && (!inv || !real_transform) )
                stage = 1;
            else if( flags & MAT_DXT_SCALE )
                scale = 1./(len * (flags & DFT_ROWS ? 1 : count));

            if( nonzero_rows <= 0 || nonzero_rows > count )
                nonzero_rows = count;

            for( i = 0; i < nonzero_rows; i++ )
            {
                uchar* sptr = src.data + i*src.step;
                uchar* dptr0 = dst.data + i*dst.step;
                uchar* dptr = dptr0;

                if( tmp_buf )
                    dptr = tmp_buf;

                dft_func( sptr, dptr, len, nf, factors, itab, wave, len, spec, ptr, _flags, scale );
                if( dptr != dptr0 )
                    memcpy( dptr0, dptr + dptr_offset, dst_full_len );
            }

            for( ; i < count; i++ )
            {
                uchar* dptr0 = dst.data + i*dst.step;
                memset( dptr0, 0, dst_full_len );
            }

			if (stage != 1)
			{
				free(buff);
				break;
			}
                
            src = dst;
        }
        else
        {
            int a = 0, b = count;
            uchar *buf0, *buf1, *dbuf0, *dbuf1;
            uchar* sptr0 = src.data;
            uchar* dptr0 = dst.data;
            buf0 = ptr;
            ptr += len*complex_elem_size;
            buf1 = ptr;
            ptr += len*complex_elem_size;
            dbuf0 = buf0, dbuf1 = buf1;

            if( use_buf )
            {
                dbuf1 = ptr;
                dbuf0 = buf1;
                ptr += len*complex_elem_size;
            }

            //dft_func = dft_tbl[(depth == CV_64F)*3];
			dft_func = dft_tbl[0];

            if( real_transform && inv && src.cols > 1 )
                stage = 0;
            else if( flags & MAT_DXT_SCALE )
                scale = 1./(len * count);

            if( real_transform )
            {
                int even;
                a = 1;
                even = (count & 1) == 0;
                b = (count+1)/2;
                if( !inv )
                {
                    memset( buf0, 0, len*complex_elem_size );
                    CopyColumn( sptr0, src.step, buf0, complex_elem_size, len, elem_size );
                    sptr0 += dst.channels()*elem_size;
                    if( even )
                    {
                        memset( buf1, 0, len*complex_elem_size );
                        CopyColumn( sptr0 + (count-2)*elem_size, src.step,
                                    buf1, complex_elem_size, len, elem_size );
                    }
                }
                else if( src.channels() == 1 )
                {
                    CopyColumn( sptr0, src.step, buf0, elem_size, len, elem_size );
                    ExpandCCS( buf0, len, elem_size );
                    if( even )
                    {
                        CopyColumn( sptr0 + (count-1)*elem_size, src.step,
                                    buf1, elem_size, len, elem_size );
                        ExpandCCS( buf1, len, elem_size );
                    }
                    sptr0 += elem_size;
                }
                else
                {
                    CopyColumn( sptr0, src.step, buf0, complex_elem_size, len, complex_elem_size );
                    if( even )
                    {
                        CopyColumn( sptr0 + b*complex_elem_size, src.step,
                                       buf1, complex_elem_size, len, complex_elem_size );
                    }
                    sptr0 += complex_elem_size;
                }

                if( even )
                    dft_func( buf1, dbuf1, len, nf, factors, itab,
                              wave, len, spec, ptr, inv, scale );
                dft_func( buf0, dbuf0, len, nf, factors, itab,
                          wave, len, spec, ptr, inv, scale );

                if( dst.channels() == 1 )
                {
                    if( !inv )
                    {
                        // copy the half of output vector to the first/last column.
                        // before doing that, defgragment the vector
                        memcpy( dbuf0 + elem_size, dbuf0, elem_size );
                        CopyColumn( dbuf0 + elem_size, elem_size, dptr0,
                                       dst.step, len, elem_size );
                        if( even )
                        {
                            memcpy( dbuf1 + elem_size, dbuf1, elem_size );
                            CopyColumn( dbuf1 + elem_size, elem_size,
                                           dptr0 + (count-1)*elem_size,
                                           dst.step, len, elem_size );
                        }
                        dptr0 += elem_size;
                    }
                    else
                    {
                        // copy the real part of the complex vector to the first/last column
                        CopyColumn( dbuf0, complex_elem_size, dptr0, dst.step, len, elem_size );
                        if( even )
                            CopyColumn( dbuf1, complex_elem_size, dptr0 + (count-1)*elem_size,
                                           dst.step, len, elem_size );
                        dptr0 += elem_size;
                    }
                }
                else
                {
                    assert( !inv );
                    CopyColumn( dbuf0, complex_elem_size, dptr0,
                                   dst.step, len, complex_elem_size );
                    if( even )
                        CopyColumn( dbuf1, complex_elem_size,
                                       dptr0 + b*complex_elem_size,
                                       dst.step, len, complex_elem_size );
                    dptr0 += complex_elem_size;
                }
            }

            for( i = a; i < b; i += 2 )
            {
                if( i+1 < b )
                {
                    CopyFrom2Columns( sptr0, src.step, buf0, buf1, len, complex_elem_size );
                    dft_func( buf1, dbuf1, len, nf, factors, itab,
                              wave, len, spec, ptr, inv, scale );
                }
                else
                    CopyColumn( sptr0, src.step, buf0, complex_elem_size, len, complex_elem_size );

                dft_func( buf0, dbuf0, len, nf, factors, itab,
                          wave, len, spec, ptr, inv, scale );

                if( i+1 < b )
                    CopyTo2Columns( dbuf0, dbuf1, dptr0, dst.step, len, complex_elem_size );
                else
                    CopyColumn( dbuf0, complex_elem_size, dptr0, dst.step, len, complex_elem_size );
                sptr0 += 2*complex_elem_size;
                dptr0 += 2*complex_elem_size;
            }

            if( stage != 0 )
            {
                if( !inv && real_transform && dst.channels() == 2 && len > 1 )
                {
                    int n = dst.cols;
                    if( elem_size == (int)sizeof(float) )
                    {
                        float* p0 = (float*)dst.data;
                        size_t dstep = dst.step/sizeof(p0[0]);
                        for( i = 0; i < len; i++ )
                        {
                            float* p = p0 + dstep*i;
                            float* q = i == 0 || i*2 == len ? p : p0 + dstep*(len-i);

                            for( int j = 1; j < (n+1)/2; j++ )
                            {
                                p[(n-j)*2] = q[j*2];
                                p[(n-j)*2+1] = -q[j*2+1];
                            }
                        }
                    }
                    else
                    {
                        double* p0 = (double*)dst.data;
                        size_t dstep = dst.step/sizeof(p0[0]);
                        for( i = 0; i < len; i++ )
                        {
                            double* p = p0 + dstep*i;
                            double* q = i == 0 || i*2 == len ? p : p0 + dstep*(len-i);

                            for( int j = 1; j < (n+1)/2; j++ )
                            {
                                p[(n-j)*2] = q[j*2];
                                p[(n-j)*2+1] = -q[j*2+1];
                            }
                        }
                    }
                }
				free(buff);
                break;
            }
            src = dst;
        }
		free(buff);
    }
}
/* End of file. */
