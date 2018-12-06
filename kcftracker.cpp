#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "kcftracker.hpp"
#include "basetype.h"
#include "fhog.hpp"
#include "ffttools.hpp"
#include "recttools.hpp"

//#include "opencv2/opencv.hpp"	//-! added for test
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

void __resizeNN_gray(const unsigned char *srcData,
	unsigned char *dstData,
	const int srcDataWidth,
	const int srcDataHeight,
	const int srcDataStep,
	const int dstDataWidth,
	const int dstDataHeight,
	const int dstDataStep) {
	int *x_ofs = (int *)malloc(sizeof(int)* dstDataWidth);
	double fx = (double)dstDataWidth / srcDataWidth;
	double fy = (double)dstDataHeight / srcDataHeight;
	double ifx = 1. / fx, ify = 1. / fy;
	int x;
	for (x = 0; x < dstDataWidth; x++) {
		int sx = floor(x * ifx);
		x_ofs[x] = min(sx, srcDataWidth - 1);
	}

	int y;
	for (y = 0; y < dstDataHeight; y++) {
		unsigned char *D = dstData + dstDataStep * y;
		int sy = min(floor(y * ify), srcDataHeight - 1);
		const unsigned char *S = srcData + srcDataStep * sy;
		for (x = 0; x <= dstDataWidth - 2; x += 2) {
			unsigned char t0 = S[x_ofs[x]];
			unsigned char t1 = S[x_ofs[x + 1]];
			D[x] = t0;
			D[x + 1] = t1;
		}
		for (; x < dstDataWidth; x++)
			D[x] = S[x_ofs[x]];
	}
	free(x_ofs);
	x_ofs = NULL;
}

// Constructor
KCFTracker::KCFTracker(int tmpl_sz, float pad_sz)
{
	// Parameters equal in all cases
	lambda = 0.0001;
	padding = pad_sz;//2.5;
	//output_sigma_factor = 0.1;
	output_sigma_factor = 0.125;

	// HOG
	interp_factor = 0.012;
	sigma = 0.6;

	cell_size = 4;

	// fit correction without multiscale
	template_size = tmpl_sz;
	scale_step = 1;
}

void KCFTracker::init(CRect &roi, CMat image)
{
	_roi.x = roi.x;
	_roi.y = roi.y;
	_roi.width = roi.width;
	_roi.height = roi.height;
	assert(roi.width >= 0 && roi.height >= 0);
	_tmpl = getFeatures(image, 1);
	_prob = createGaussianPeak(size_patch[0], size_patch[1]);
	//printCMat2f(_prob);
	_alphaf = CMat(size_patch[0], size_patch[1], MAT_32FC2);
	memset(_alphaf.data,0,_alphaf.rows*_alphaf.step);	//init memory, may be wrong
	//printCMat2f(_alphaf);
	train(_tmpl, 1.0); // train with initial frame
}

// Update position based on the new frame
CRect KCFTracker::update(CMat image)
{
	if (_roi.x + _roi.width <= 0) _roi.x = -_roi.width + 1;
	if (_roi.y + _roi.height <= 0) _roi.y = -_roi.height + 1;
	if (_roi.x >= image.cols - 1) _roi.x = image.cols - 2;
	if (_roi.y >= image.rows - 1) _roi.y = image.rows - 2;

	float cx = _roi.x + _roi.width / 2.0f;
	float cy = _roi.y + _roi.height / 2.0f;


	float peak_value;
	CPoint2f res = detect(_tmpl, getFeatures(image, 0, 1.0f), peak_value);

	//if (scale_step != 1) {
	//	// Test at a smaller _scale
	//	float new_peak_value;
	//	CPoint2f new_res = detect(_tmpl, getFeatures(image, 0, 1.0f / scale_step), new_peak_value);

	//	if (scale_weight * new_peak_value > peak_value) {
	//		res = new_res;
	//		peak_value = new_peak_value;
	//		_scale /= scale_step;
	//		_roi.width /= scale_step;
	//		_roi.height /= scale_step;
	//	}

	//	// Test at a bigger _scale
	//	new_res = detect(_tmpl, getFeatures(image, 0, scale_step), new_peak_value);

	//	if (scale_weight * new_peak_value > peak_value) {
	//		res = new_res;
	//		peak_value = new_peak_value;
	//		_scale *= scale_step;
	//		_roi.width *= scale_step;
	//		_roi.height *= scale_step;
	//	}
	//}

	// Adjust by cell size and _scale
	_roi.x = cx - _roi.width / 2.0f + ((float)res.x * cell_size * _scale);
	_roi.y = cy - _roi.height / 2.0f + ((float)res.y * cell_size * _scale);

	if (_roi.x >= image.cols - 1) _roi.x = image.cols - 1;
	if (_roi.y >= image.rows - 1) _roi.y = image.rows - 1;
	if (_roi.x + _roi.width <= 0) _roi.x = -_roi.width + 2;
	if (_roi.y + _roi.height <= 0) _roi.y = -_roi.height + 2;

	assert(_roi.width >= 0 && _roi.height >= 0);
	CMat x = getFeatures(image, 0);
	train(x, interp_factor);

	CRect roi_ret(round(_roi.x), round(_roi.y), (int)_roi.width, (int)_roi.height);
	return roi_ret;
}

// Detect object in the current frame.
CPoint2f KCFTracker::detect(CMat z, CMat x, float &peak_value)
{
	using namespace FFTTools;

	CMat k = gaussianCorrelation(x, z, false);
	CMat res = (real(fftd(complexMultiplication(_alphaf, fftd(k)), true)));
	//printCMat2f(complexMultiplication(_alphaf, fftd(k)));
	//printCMat2f(res);
	CPoint pi(0, 0);
	peak_value = 0.0;
	for (int i = 0; i < res.rows; i++)
	{
		float *pres = (float*)res.ptr(i);
		for (int j = 0; j < res.cols; j++)
		{
			//printf("%f ", pres[j]);
			if (pres[j] > peak_value)
			{
				peak_value = pres[j];
				pi.x = j;
				pi.y = i;
			}
		}
		//printf("\n");
	}
	//minMaxLoc only accepts doubles for the peak, and integer points for the coordinates
	/*cv::Point2i pi;
	double pv;
	cv::minMaxLoc(res, NULL, &pv, NULL, &pi);
	peak_value = (float)pv;*/

	//subpixel peak estimation, coordinates will be non-integer
	CPoint2f p((float)pi.x, (float)pi.y);

	if (pi.x > 0 && pi.x < res.cols - 1) {
		p.x += subPixelPeak(((float*)res.ptr(pi.y, pi.x - 1))[0], peak_value, ((float*)res.ptr(pi.y, pi.x + 1))[0]);
		//p.x += subPixelPeak(res.at<float>(pi.y, pi.x - 1), peak_value, res.at<float>(pi.y, pi.x + 1));
	}

	if (pi.y > 0 && pi.y < res.rows - 1) {
		p.y += subPixelPeak(((float*)res.ptr(pi.y - 1, pi.x))[0], peak_value, ((float*)res.ptr(pi.y +1 , pi.x))[0]);
		//p.y += subPixelPeak(res.at<float>(pi.y - 1, pi.x), peak_value, res.at<float>(pi.y + 1, pi.x));
	}

	p.x -= (res.cols) / 2;
	p.y -= (res.rows) / 2;

	return p;
}

// train tracker with a single image
void KCFTracker::train(CMat x, float train_interp_factor)
{
	using namespace FFTTools;

	CMat k = gaussianCorrelation(x, x, true);
	CMat alphaf = complexDivision(_prob, (fftd(k) + lambda));

	_tmpl = _tmpl * (1 - train_interp_factor) + x * (train_interp_factor);
	_alphaf = _alphaf * (1 - train_interp_factor) + alphaf * (train_interp_factor);

	//printCMat2f(fftd(k)+lambda);
	//printCMat2f(_prob);
	//printCMat2f(alphaf);
	//printCMat2f(_tmpl);
	//printCMat2f(_alphaf);
}

// Evaluates a Gaussian kernel with bandwidth SIGMA for all relative shifts between input images X and Y, which must both be MxN. They must    also be periodic (ie., pre-processed with a cosine window).
CMat KCFTracker::gaussianCorrelation(CMat x1, CMat x2, bool self_correlation)
{
	using namespace FFTTools;
	CMat c = CMat(CSize(size_patch[1], size_patch[0]), MAT_32F);
	memset(c.data, 0, c.rows*c.step);	//init memory, may be wrong

	//printCMat2f(c);
	// HOG features

	if (!self_correlation) {
		//printCMat2f(x1);
		//printCMat2f(x2);
		for (int i = 0; i < size_patch[2]; i++) {
			float *px1 = (float*)x1.ptr(i);
			float *px2 = (float*)x2.ptr(i);
			CMat x1row = CMat(CSize(size_patch[1], size_patch[0]), MAT_32FC2);
			CMat x2row = CMat(CSize(size_patch[1], size_patch[0]), MAT_32FC2);
			for (int j = 0; j < size_patch[0]; j++){
				float *px1row = (float*)x1row.ptr(j);
				float *px2row = (float*)x2row.ptr(j);
				for (int k = 0; k < size_patch[1]; k++){
					px1row[2 * k] = px1[j * size_patch[1] + k];
					px1row[2 * k + 1] = 0;
					px2row[2 * k] = px2[j * size_patch[1] + k];
					px2row[2 * k + 1] = 0;
				}
			}
			CMat x1aux = fftd(x1row);
			CMat x2aux = fftd(x2row);
			CMat caux = fftd(complexMultiplication(x1aux, x2aux, true), true);
			//printCMat2f(x1aux);
			//printCMat2f(x2aux);
			//printCMat2f(complexMultiplication(x1aux, x2aux, true));
			//printCMat2f(caux);
			rearrange(caux);
			c = c + real(caux);
		}
	}
	else{
		for (int i = 0; i < size_patch[2]; i++) {
			float *px1 = (float*)x1.ptr(i);
			CMat x1row = CMat(CSize(size_patch[1], size_patch[0]), MAT_32FC2);
			for (int j = 0; j < size_patch[0]; j++){
				float *px1row = (float*)x1row.ptr(j);
				for (int k = 0; k < size_patch[1];k++){
					px1row[2 * k] = px1[j * size_patch[1] + k];
					px1row[2 * k + 1] = 0;
				}
			}
			CMat x1aux = fftd(x1row);
			//printCMat2f(x1aux);
			CMat caux = fftd(complexMultiplication(x1aux, x1aux, true),true);
			//printCMat2f(complexMultiplication(x1aux, x1aux, true));
			//printCMat2f(caux);
			rearrange(caux);
			c = c + real(caux);
			//printCMat2f(c);
		}
	}

	//printCMat2f(c);
	float sum_x1 = 0.0;
	float sum_x2 = 0.0;
	for (int i = 0; i < x1.rows; i++)
	{
		float *pdata_x1 = (float *)x1.ptr(i);
		for (int j = 0; j < x1.cols; j++)
		{
			sum_x1 += pdata_x1[j] * pdata_x1[j];
		}
		//pdata_x1 += x1.step / sizeof(float);
	}
	if (self_correlation){
		sum_x2 = sum_x1;
	}
	else{
		for (int i = 0; i < x1.rows; i++)
		{
			float *pdata_x2 = (float *)x2.ptr(i);
			for (int j = 0; j < x1.cols; j++)
			{
				sum_x2 += pdata_x2[j] * pdata_x2[j];
			}
		}
	}
	
	CMat d = c*(-2.0) + (sum_x1 + sum_x2);
	CMat k = CMat(CSize(size_patch[1], size_patch[0]), MAT_32FC2);
	int sum_size = size_patch[0] * size_patch[1] * size_patch[2];
	float denominator = sigma*sigma*sum_size;
	for (int i = 0; i < d.rows; i++)
	{
		float *pd = (float*)d.ptr(i);
		float *pk = (float*)k.ptr(i);
		for (int j = 0; j < d.cols; j++)
		{
			float tmp = pd[j];
			if (tmp > 0){
				pk[2 * j] = exp(tmp*(-1) / denominator);
			}
			else{
				pk[2 * j] = 1.0;
			}
			pk[2 * j + 1] = 0.0;
		}
	}
	return k;
	//cv::Mat d;
	//cv::Scalar sum_x1 = cv::sum(x1.mul(x1));
	//cv::Scalar sum_x2 = cv::sum(x2.mul(x2));
	//cv::max(((sum_x1[0] + sum_x2[0]) - 2. * c) / (size_patch[0] * size_patch[1] * size_patch[2]), 0, d);
	////cv::max(((cv::sum(x1.mul(x1))[0] + cv::sum(x2.mul(x2))[0]) - 2. * c) / (size_patch[0] * size_patch[1] * size_patch[2]), 0, d);

	//cv::Mat k;
	//cv::exp((-d / (sigma * sigma)), k);
	//return k;	
}

// Create Gaussian Peak. Function called only in the first frame.
CMat KCFTracker::createGaussianPeak(int sizey, int sizex)
{
	CMat res(sizey, sizex, MAT_32FC2);

	int syh = (sizey) / 2;
	int sxh = (sizex) / 2;

	float output_sigma = sqrt((float)sizex * sizey) / padding * output_sigma_factor;
	float mult = -0.5 / (output_sigma * output_sigma);

	
	for (int i = 0; i < sizey; i++)
	{
		float* pres_data = (float*)(res.data+i*res.step);
		for (int j = 0; j < sizex; j++)
		{
			int ih = i - syh;
			int jh = j - sxh;
			pres_data[2 * j] = exp(mult * (float)(ih * ih + jh * jh));
			pres_data[2 * j + 1] = 0.0;
		}
	}
	//printCMat2f(FFTTools::real(res));
	return FFTTools::fftd(res);
}

// Obtain sub-window from image, with replication-padding and extract features
CMat KCFTracker::getFeatures(const CMat & image, bool inithann, float scale_adjust)
{
	CRect extracted_roi;

	float cx = _roi.x + _roi.width / 2;
	float cy = _roi.y + _roi.height / 2;

	if (inithann) {
		int padded_w = _roi.width * padding;
		int padded_h = _roi.height * padding;

		if (template_size > 1) {  // Fit largest dimension to the given template size
			if (padded_w >= padded_h)  //fit to width
				_scale = padded_w / (float)template_size;
			else
				_scale = padded_h / (float)template_size;

			_tmpl_sz.width = padded_w / _scale;
			_tmpl_sz.height = padded_h / _scale;
		}
		else {  //No template size given, use ROI size
			_tmpl_sz.width = padded_w;
			_tmpl_sz.height = padded_h;
			_scale = 1;
		}

		if (true) {
			// Round to cell size and also make it even
			_tmpl_sz.width = (((int)(_tmpl_sz.width / (2 * cell_size))) * 2 * cell_size) + cell_size * 2;
			_tmpl_sz.height = (((int)(_tmpl_sz.height / (2 * cell_size))) * 2 * cell_size) + cell_size * 2;
		}
		else {  //Make number of pixels even (helps with some logic involving half-dimensions)
			_tmpl_sz.width = (_tmpl_sz.width / 2) * 2;
			_tmpl_sz.height = (_tmpl_sz.height / 2) * 2;
		}
	}

	extracted_roi.width = scale_adjust * _scale * _tmpl_sz.width;
	extracted_roi.height = scale_adjust * _scale * _tmpl_sz.height;

	// center roi with new size
	extracted_roi.x = cx - extracted_roi.width / 2;
	extracted_roi.y = cy - extracted_roi.height / 2;

	CMat FeaturesMap;
	CMat z = RectTools::subwindow(image, extracted_roi);

	//cv::Mat cv_z = cv::Mat(z.rows, z.cols, CV_8U, z.data);
	//cv::imshow("test_z", cv_z);
	//cv::waitKey(0);

	CMat resized_z;// = CMat(_tmpl_sz, MAT_8U);
	if (z.cols != _tmpl_sz.width || z.rows != _tmpl_sz.height) {
		//cv::resize(z, z, _tmpl_sz, CV_INTER_NN);
		resized_z.create(_tmpl_sz, MAT_8U);
		__resizeNN_gray(z.data, resized_z.data, z.cols, z.rows, z.step, resized_z.cols, resized_z.rows, resized_z.step);
	}
	else{
		resized_z = z;
	}
	//printCMat2i(resized_z);
	//cv::Mat cv_z = cv::Mat(resized_z.rows, resized_z.cols, CV_8U, resized_z.data);
	//cv::imshow("test_z", cv_z);
	//cv::waitKey(0);

	// HOG features
	if (true) {
		//IplImage z_ipl = z;
		CvLSVMFeatureMapCaskade *map;
		getFeatureMaps(resized_z, cell_size, &map);
		normalizeAndTruncate(map, 0.2f);
		PCAFeatureMaps(map);
		size_patch[0] = map->sizeY;
		size_patch[1] = map->sizeX;
		size_patch[2] = map->numFeatures;
		//for (int i = 0; i < size_patch[0] * size_patch[1] * size_patch[2]; i++)
		//{
		//	printf("%f ", map->map[i]);
		//}
		//FeaturesMap = cv::Mat(cv::Size(map->numFeatures, map->sizeX*map->sizeY), CV_32F, map->map);  // Procedure do deal with cv::Mat multichannel bug
		//FeaturesMap = FeaturesMap.t();
		FeaturesMap.create(CSize(map->sizeX*map->sizeY, map->numFeatures), MAT_32F);
		for (int i = 0; i < FeaturesMap.rows; i++)
		{
			float *pfm = (float*)FeaturesMap.ptr(i);
			for (int j = 0; j < FeaturesMap.cols; j++)
			{
				pfm[j] = map->map[j*FeaturesMap.rows + i];
			}
		}
		//printCMat2f(FeaturesMap);
		freeFeatureMapObject(&map);
	}	

	if (inithann) {
		createHanningMats();
	}
	FeaturesMap = hann.mul(FeaturesMap);
	//printCMat2f(FeaturesMap);
	//printCMat2f(hann);
	return FeaturesMap;
}

// Initialize Hanning window. Function called only in the first frame.
void KCFTracker::createHanningMats()
{
	CMat hann1t = CMat(CSize(size_patch[1], 1), MAT_32F);
	CMat hann2t = CMat(CSize(size_patch[0], 1), MAT_32F);

	float *ph1 = (float*)hann1t.data;
	float *ph2 = (float*)hann2t.data;
	for (int i = 0; i < hann1t.cols; i++)
		ph1[i] = 0.5 * (1 - cos(2 * 3.14159265358979323846 * i / (hann1t.cols - 1)));
		//hann1t.at<float >(0, i) = 0.5 * (1 - std::cos(2 * 3.14159265358979323846 * i / (hann1t.cols - 1)));
	for (int i = 0; i < hann2t.cols; i++)
		ph2[i] = 0.5 * (1 - cos(2 * 3.14159265358979323846 * i / (hann2t.cols - 1)));
		//hann2t.at<float >(i, 0) = 0.5 * (1 - std::cos(2 * 3.14159265358979323846 * i / (hann2t.rows - 1)));

	CMat hann2d = CMat(CSize(size_patch[0]*size_patch[1], 1), MAT_32F);
	float *pha2d = (float*)hann2d.data;
	for (int i = 0; i < size_patch[0]; i++)
	{
		for (int j = 0; j < size_patch[1]; j++)
		{
			pha2d[j] = ph1[j] * ph2[i];
		}
		pha2d += size_patch[1];
	}

	// HOG features
	hann.create(CSize(size_patch[0] * size_patch[1], size_patch[2]), MAT_32F);
	pha2d = (float*)hann2d.data;	
	for (int i = 0; i < size_patch[2]; i++)
	{
		float *phann = (float*)hann.ptr(i);
		memcpy(phann, pha2d, hann2d.cols*sizeof(float));
	}
	
	//if (true) {
	//	cv::Mat hann1d = hann2d.reshape(1, 1); // Procedure do deal with cv::Mat multichannel bug
	//	hann = cv::Mat(cv::Size(size_patch[0] * size_patch[1], size_patch[2]), CV_32F, cv::Scalar(0));
	//	for (int i = 0; i < size_patch[2]; i++) {
	//		for (int j = 0; j<size_patch[0] * size_patch[1]; j++) {
	//			hann.at<float>(i, j) = hann1d.at<float>(0, j);
	//		}
	//	}
	//}
	//// Gray features
	//else {
	//	hann = hann2d;
	//}
}

// Calculate sub-pixel peak for one dimension
float KCFTracker::subPixelPeak(float left, float center, float right)
{
	float divisor = 2 * center - right - left;

	if (divisor == 0)
		return 0;

	return 0.5 * (right - left) / divisor;
}
