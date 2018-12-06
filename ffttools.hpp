/* 
Author: Christian Bailer
Contact address: Christian.Bailer@dfki.de 
Department Augmented Vision DFKI 

                          License Agreement
               For Open Source Computer Vision Library
                       (3-clause BSD License)

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the names of the copyright holders nor the names of the contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

This software is provided by the copyright holders and contributors "as is" and
any express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are disclaimed.
In no event shall copyright holders or contributors be liable for any direct,
indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused
and on any theory of liability, whether in contract, strict liability,
or tort (including negligence or otherwise) arising in any way out of
the use of this software, even if advised of the possibility of such damage.
*/

#pragma once

//#include <cv.h>

//#ifndef _OPENCV_FFTTOOLS_HPP_
//#define _OPENCV_FFTTOOLS_HPP_
//#endif
#include <assert.h>
#include "basetype.h"
#include "dxt.h"

namespace FFTTools
{
// Previous declarations, to avoid warnings
CMat fftd(CMat img, bool backwards = false);
CMat real(CMat img);
CMat imag(CMat img);
CMat complexMultiplication(CMat a, CMat b, bool bconj = false);
CMat complexDivision(CMat a, CMat b);
void rearrange(CMat &img);


CMat fftd(CMat img, bool backwards)
{
	//if (img.channels() == 1)	//channels() has to be 2
	//{
	//    cv::Mat planes[] = {cv::Mat_<float> (img), cv::Mat_<float>::zeros(img.size())};
	//    //cv::Mat planes[] = {cv::Mat_<double> (img), cv::Mat_<double>::zeros(img.size())};
	//    cv::merge(planes, 2, img);
	//}
	assert(img.channels() == 2);
	CMat res;
	dxt::dft(img, res, backwards ? (dxt::DFT_INVERSE | dxt::DFT_SCALE) : 0);
	//for (int i = 0; i < res.rows; i++)
	//{
	//	float* pres = (float*)res.ptr(i);
	//	for (int j = 0; j < res.cols; j++)
	//	{
	//		pres[2 * j + 1] *= (-1);
	//	}
	//}
	return res;
}

CMat real(CMat img)
{
	assert(img.channels() == 2);
	CMat realimg = CMat(img.rows, img.cols, MAT_32FC1);
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			//printf("(%f,%f)", ((float*)(img.ptr(i, j)))[0], ((float*)(img.ptr(i, j)))[1]);
			((float*)(realimg.ptr(i, j)))[0] = ((float*)(img.ptr(i, j)))[0];
		}
		//printf("\n");
	}
	return realimg;
}

CMat imag(CMat img)
{
	assert(img.channels() == 2);
	CMat imagimg = CMat(img.rows, img.cols, MAT_32FC1);
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			((float*)(imagimg.ptr(i, j)))[0] = ((float*)(img.ptr(i, j)))[1];
		}
	}
	return imagimg;
}


CMat complexMultiplication(CMat a, CMat b, bool bconj)
{
	assert(a.cols == b.cols && a.rows == b.rows && a.channels() == b.channels());
	CMat res = CMat(a.rows, a.cols, MAT_32FC2);

	for (int i = 0; i < a.rows; i++)
	{
		float *pres = (float*)res.ptr(i);
		float *pa = (float*)a.ptr(i);
		float *pb = (float*)b.ptr(i);
		for (int j = 0; j < a.cols; j++)
		{
			int k = 2 * j;
			if (!bconj)
			{
				pres[k] = (pa[k] * pb[k] - pa[k + 1] * pb[k + 1]);
				pres[k + 1] = (pa[k + 1] * pb[k] + pa[k] * pb[k + 1]);
			}
			else{
				pres[k] = (pa[k] * pb[k] + pa[k + 1] * pb[k + 1]);
				pres[k + 1] = (pa[k + 1] * pb[k] - pa[k] * pb[k + 1]);
				//pres[k + 1] = (pa[k] * pb[k + 1] - pa[k + 1] * pb[k]);
			}
		}
	}

	return res;
}

CMat complexDivision(CMat a, CMat b)
{
	assert(a.cols == b.cols && a.rows == b.rows && a.channels() == b.channels());
	CMat res = CMat(a.rows,a.cols,MAT_32FC2);

	for (int i = 0; i < a.rows; i++)
	{
		float *pres = (float*)res.ptr(i);
		float *pa = (float*)a.ptr(i);
		float *pb = (float*)b.ptr(i);
		for (int j = 0; j < a.cols; j++)
		{
			int k = 2 * j;
			float divisor = 1.0 / (pb[k] * pb[k] + pb[k + 1] * pb[k + 1]);
			pres[k] = (pa[k] * pb[k] + pa[k + 1] * pb[k + 1])*divisor;
			pres[k + 1] = (pa[k + 1] * pb[k] - pa[k] * pb[k + 1])*divisor;
		}
	}

    return res;
}

void rearrange(CMat &img)
{
    // img = img(cv::Rect(0, 0, img.cols & -2, img.rows & -2));
    int cx = img.cols / 2;
    int cy = img.rows / 2;

	int cpy_bytes = cx*img.step[1];
	void* ptemp = malloc(cpy_bytes);
	// swap quadrants (Top-Left with Bottom-Right)
	for (int i = 0; i < cy; i++)
	{
		void* ptl = (void*)img.ptr(i);
		void* pbr = (void*)img.ptr(cy+i, cx);
		memcpy(ptemp, ptl, cpy_bytes);
		memcpy(ptl, pbr, cpy_bytes);
		memcpy(pbr, ptemp, cpy_bytes);
	}
	// swap quadrant (Top-Right with Bottom-Left)
	for (int i = 0; i < cy; i++)
	{
		void* ptr = (void*)img.ptr(i,cx);
		void* pbl = (void*)img.ptr(cy + i);
		memcpy(ptemp, ptr, cpy_bytes);
		memcpy(ptr, pbl, cpy_bytes);
		memcpy(pbl, ptemp, cpy_bytes);
	}
	free(ptemp);
    //cv::Mat q0(img, cv::Rect(0, 0, cx, cy)); // Top-Left - Create a ROI per quadrant
    //cv::Mat q1(img, cv::Rect(cx, 0, cx, cy)); // Top-Right
    //cv::Mat q2(img, cv::Rect(0, cy, cx, cy)); // Bottom-Left
    //cv::Mat q3(img, cv::Rect(cx, cy, cx, cy)); // Bottom-Right
    //cv::Mat tmp; // swap quadrants (Top-Left with Bottom-Right)
    //q0.copyTo(tmp);
    //q3.copyTo(q0);
    //tmp.copyTo(q3);
    //q1.copyTo(tmp); // swap quadrant (Top-Right with Bottom-Left)
    //q2.copyTo(q1);
    //tmp.copyTo(q2);
}
}
