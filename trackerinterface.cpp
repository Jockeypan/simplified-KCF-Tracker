#include "kcftracker.hpp"
#include "trackerinterface.h"

TrackHandle CreateTracker(int template_size, float padding_size)
{
    KCFTracker* tracker = new KCFTracker(template_size,padding_size);
    return (TrackHandle)tracker;
}

int InitTracker(TrackHandle th, BBox roi, unsigned char* gray, const int width, const int height, const int step)
{
    KCFTracker* tracker = (KCFTracker*)th;
    CMat img = CMat(height,width,MAT_8UC1,(void*)gray,step);
    //img.data = gray;
    CRect box(roi.x,roi.y,roi.w,roi.h);
    //cv::imshow("image", img );
    //cv::waitKey(0);
    tracker->init( box, img );
    
    return 0;
}

BBox UpdateTracker(TrackHandle th,unsigned char* gray, const int width, const int height,const int step)
{
    KCFTracker* tracker = (KCFTracker*)th;
    CMat img = CMat(height,width,MAT_8UC1,(void*)gray,step);
    //img.data = gray;
    CRect bb = tracker->update(img);
    BBox box;
    box.x = bb.x;
    box.y = bb.y;
    box.w = bb.width;
    box.h = bb.height;

    return box;
}

bool ReleaseTracker(TrackHandle* pth)
{
    KCFTracker** ptracker = (KCFTracker**)pth;
    if (ptracker!=NULL && *ptracker != NULL)
    {
        delete *ptracker;
        *ptracker = NULL;
        ptracker = NULL;
    }

    return true;
}
