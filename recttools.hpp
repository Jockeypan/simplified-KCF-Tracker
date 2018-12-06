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
#include <math.h>
#include <string.h>
#include <assert.h>
#include "basetype.h"

namespace RectTools
{
inline CMat subwindow(const CMat &in, const CRect & window)
{
	assert(in.channels() == 1);
	CMat submat = CMat(window.height, window.width, MAT_8UC1);

	CRect comm, diff;
	int win_max_x = window.x + window.width;
	int win_max_y = window.y + window.height;
	assert(window.x < in.cols && window.y<in.rows && win_max_x>0&&win_max_y>0);
	comm.x = window.x >= 0 ? window.x : 0;
	comm.y = window.y >= 0 ? window.y : 0;
	diff.x = comm.x - window.x;
	diff.y = comm.y - window.y;
	if (win_max_x > in.cols)
	{
		comm.width = in.cols - comm.x;
		diff.width = win_max_x - in.cols;
	}
	else{
		comm.width = win_max_x - comm.x;
		diff.width = 0;
	}
	if (win_max_y > in.rows)
	{
		comm.height = in.rows - comm.y;
		diff.height = win_max_y - in.rows;
	}
	else{
		comm.height = win_max_y - comm.y;
		diff.height = 0;
	}

	uchar *psub = NULL;
	const uchar	*pin = NULL;
	//common area
	psub = submat.ptr(diff.y, diff.x);
	pin = in.ptr(comm.y, comm.x);
	for (int i = 0; i < comm.height; i++)
	{
		memcpy(psub, pin, comm.width);
		psub += submat.step;
		pin += in.step;
	}
	//left
	if (diff.x>0)
	{
		for (int i = 0; i < comm.height; i++)
		{
			psub = submat.ptr(diff.y+i);
			uchar padval = psub[diff.x];
			memset(psub, padval, diff.x);
		}
		//psub += submat.step;
	}
	//right
	if (diff.width > 0)
	{
		for (int i = 0; i < comm.height; i++)
		{
			psub = submat.ptr(diff.y + i, diff.x);
			uchar padval = psub[comm.width-1];
			memset(psub+comm.width, padval, diff.width);
		}
	}

	uchar* psub_src = NULL;
	//up
	if (diff.y > 0)
	{
		psub_src = submat.ptr(diff.y, 0);
		for (int i = 0; i < diff.y; i++)
		{
			psub = submat.ptr(i);
			memcpy(psub, psub_src, submat.cols);
		}
	}
	//down
	if (diff.height > 0)
	{
		psub_src = submat.ptr(diff.y + comm.height - 1);
		for (int i = 0; i < diff.height; i++)
		{
			psub = submat.ptr(diff.y + comm.height + i);
			memcpy(psub, psub_src, submat.cols);
		}
	}

	return submat;
}
}



