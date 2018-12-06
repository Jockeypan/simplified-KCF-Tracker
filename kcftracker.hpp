
#ifndef _KCFTRACKER_HPP_
#define _KCFTRACKER_HPP_


#include "basetype.h"

namespace kcf
{
	struct fBBox
	{
		float x, y, width, height;
	};
}
class KCFTracker
{
public:
public:
	// Constructor
	KCFTracker(int tmpl_sz = 96, float pad_sz = 2.0);

	// Initialize tracker 
	void init(CRect &roi, CMat image);

	// Update position based on the new frame
	CRect update(CMat image);

protected:
	// Detect object in the current frame.
	CPoint2f detect(CMat z, CMat x, float &peak_value);

	// train tracker with a single image
	void train(CMat x, float train_interp_factor);

	// Evaluates a Gaussian kernel with bandwidth SIGMA for all relative shifts between input images X and Y, which must both be MxN. They must    also be periodic (ie., pre-processed with a cosine window).
	CMat gaussianCorrelation(CMat x1, CMat x2, bool self_correlation = false);

	// Create Gaussian Peak. Function called only in the first frame.
	CMat createGaussianPeak(int sizey, int sizex);

	// Obtain sub-window from image, with replication-padding and extract features
	CMat getFeatures(const CMat & image, bool inithann, float scale_adjust = 1.0f);

	// Initialize Hanning window. Function called only in the first frame.
	void createHanningMats();

	// Calculate sub-pixel peak for one dimension
	float subPixelPeak(float left, float center, float right);
private:
	float interp_factor; // linear interpolation factor for adaptation
	float sigma; // gaussian kernel bandwidth
	float lambda; // regularization
	int cell_size; // HOG cell size
	float padding; // extra area surrounding the target
	float output_sigma_factor; // bandwidth of gaussian target
	int template_size; // template size
	float scale_step; // scale step for multi-scale estimation

	//int _width, _height; // image width & height
	int size_patch[3];
	CMat hann;
	float _scale;

	CRect2f _roi;
	CSize _tmpl_sz;
	CMat _alphaf;
	CMat _prob;
	CMat _tmpl;
};
#endif