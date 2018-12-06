#ifndef _TRACK_INTERFACE_H_
#define _TRACK_INTERFACE_H_
struct BBox
{
    int x, y, w, h;
};

typedef void* TrackHandle;
TrackHandle CreateTracker(int template_size,float padding_size);
int InitTracker(TrackHandle th, BBox Roi, unsigned char* gray, const int width, const int height,const int step);
BBox UpdateTracker(TrackHandle th,unsigned char* gray, const int width, const int height,const int step);
bool ReleaseTracker(TrackHandle* pth);
#endif
