#include <stdlib.h>

#include "kcftracker.hpp"
#include "basetype.h"
#include <iostream>
#include <string>
//#include "trackerinterface.h"
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;
cv::Rect box;
bool drawing_box = false;
bool selected = false;
void create_mouse_callback(int event, int x, int y, int flag, void* param)
{
	cv::Mat *image = (cv::Mat*) param;
	switch (event){
	case CV_EVENT_MOUSEMOVE:
		if (drawing_box){
			box.width = x - box.x;
			box.height = y - box.y;
		}
		break;

	case CV_EVENT_LBUTTONDOWN:
		drawing_box = true;
		box = cv::Rect(x, y, 0, 0);
		break;

	case CV_EVENT_LBUTTONUP:
		drawing_box = false;
		if (box.width < 0){
			box.x += box.width;
			box.width *= -1;
		}
		if (box.height < 0){
			box.y += box.height;
			box.height *= -1;
		}
		cv::rectangle(*image, box, cv::Scalar(0), 2);
		selected = true;
		break;
	}

}
int main()
{

	KCFTracker tracker(256);
	//TrackHandle tracker = CreateTracker(72);
	cv::Mat frame;

	string img_root = "./DragonBaby/";
	string img_path = img_root + "/img/";
	string first_frame_name = img_path + "0001.jpg";
	frame = imread(first_frame_name, 1);
	cv::namedWindow("original image");
	Mat temp_img = frame.clone();
	/*box = cv::Rect(400, 47, 93, 174);

	cv::Mat temp;
	temp_img.copyTo(temp);
	cv::rectangle(temp, box, cv::Scalar(0), 2);*/
	/*cv::imshow("original image", temp);
	cv::waitKey(15);*/
	cv::setMouseCallback("original image", create_mouse_callback, (void*)&temp_img);
	cv::imshow("original image", frame);
	while (selected == false)
	{
		cv::Mat temp;

		temp_img.copyTo(temp);

		if (drawing_box)
			cv::rectangle(temp, box, cv::Scalar(0), 2);

		cv::imshow("original image", temp);

		if (cv::waitKey(15) == 27)
			break;
	}
	cv::setMouseCallback("original image", NULL, NULL);
	waitKey(0);
	if (box.width == 0 || box.height == 0)
		return 0;
	//tracker.init(frame, box);
	Mat frame_gray;
	cv::cvtColor(frame, frame_gray, CV_BGR2GRAY);

	CMat cgray = CMat(CSize(frame_gray.cols, frame_gray.rows), MAT_8UC1, frame_gray.data);
	CRect cbox(box.x, box.y, box.width, box.height);
	tracker.init(cbox, cgray);
	/*BBox initBox;
	initBox.x = box.x;
	initBox.y = box.y;
	initBox.w = box.width;
	initBox.h = box.height;
	InitTracker(tracker,initBox,frame_gray.data,frame_gray.cols,frame_gray.rows,frame_gray.step);*/

	printf("Start the tracking process\n");
	char img_name[128];
	int img_idx = 2;
	//BBox_c bb;
	for (;;) {
		// get frame from the video
		//cap >> frame;
		sprintf(img_name, "%04d.jpg", img_idx++);
		string img_full_name = img_path + img_name;
		frame = imread(img_full_name, 1);
		if (frame.empty()) break;

		// stop the program if no more images
		if (frame.rows == 0 || frame.cols == 0)
			break;

		// update the tracking result
		//! [update]
		//tracker->update(frame, box);
		Mat frame_gray;
		cv::cvtColor(frame, frame_gray, CV_BGR2GRAY);
		CMat cgray = CMat(CSize(frame_gray.cols, frame_gray.rows), MAT_8UC1, frame_gray.data);
		double time_profile_counter = (double)cvGetTickCount();
		CRect cbb = tracker.update(cgray);
		Rect bb(cbb.x, cbb.y, cbb.width, cbb.height);
		/*BBox bbox = UpdateTracker(tracker,frame_gray.data,frame_gray.cols,frame_gray.rows,frame_gray.step);
		Rect bb(bbox.x,bbox.y,bbox.w,bbox.h);*/

		time_profile_counter = (double)cvGetTickCount() - time_profile_counter;
		std::cout << "time used " << time_profile_counter / ((double)cvGetTickFrequency() * 1000) << "ms" << endl;
		//! [update]

		//! [visualization]
		// draw the tracked object
		rectangle(frame, bb, Scalar(255, 0, 0), 2, 1);

		// show image with the tracked object
		imshow("tracker", frame);
		//! [visualization]
		//quit on ESC button
		if (waitKey(0) == 27)
			break;
	}

	//ReleaseTracker(&tracker);
	return EXIT_SUCCESS;
}

int main2()
{

	KCFTracker tracker(96,2.0);
	//TrackHandle tracker = CreateTracker(72);
	cv::Mat frame;
	VideoCapture cap(
		"C:/Users/jimmypang/Desktop/testVideo/clip/ice_video_20170412-144949-5.webm"
		);
	cap >> frame;
	//string img_root = "D:/WorkSpace/kcf-master/BlurBody/";
	//string img_path = img_root + "/img/";
	//string first_frame_name = img_path + "0001.jpg";
	//frame = imread(first_frame_name, 1);
	cv::namedWindow("original image");
	Mat temp_img = frame.clone();
	/*box = cv::Rect(400, 47, 93, 174);

	cv::Mat temp;
	temp_img.copyTo(temp);
	cv::rectangle(temp, box, cv::Scalar(0), 2);*/
	/*cv::imshow("original image", temp);
	cv::waitKey(15);*/
	cv::setMouseCallback("original image", create_mouse_callback, (void*)&temp_img);
	cv::imshow("original image", frame);
	while (selected == false)
	{
		cv::Mat temp;

		temp_img.copyTo(temp);

		if (drawing_box)
			cv::rectangle(temp, box, cv::Scalar(0), 2);

		cv::imshow("original image", temp);

		if (cv::waitKey(15) == 27)
			break;
	}
	cv::setMouseCallback("original image", NULL, NULL);
	waitKey(0);
	if (box.width == 0 || box.height == 0)
		return 0;
	//tracker.init(frame, box);
	Mat frame_gray;
	cv::cvtColor(frame, frame_gray, CV_BGR2GRAY);

	CMat cgray = CMat(CSize(frame_gray.cols, frame_gray.rows), MAT_8UC1, frame_gray.data);
	CRect cbox(box.x, box.y, box.width, box.height);
	tracker.init(cbox, cgray);
	/*BBox initBox;
	initBox.x = box.x;
	initBox.y = box.y;
	initBox.w = box.width;
	initBox.h = box.height;
	InitTracker(tracker,initBox,frame_gray.data,frame_gray.cols,frame_gray.rows,frame_gray.step);*/

	printf("Start the tracking process\n");
	char img_name[128];
	int img_idx = 2;
	//BBox_c bb;
	for (;;) {
		// get frame from the video
		//cap >> frame;
		//sprintf_s(img_name, "%04d.jpg", img_idx++);
		//string img_full_name = img_path + img_name;
		//frame = imread(img_full_name, 1);
		cap >> frame;
		if (frame.empty()) break;

		// stop the program if no more images
		if (frame.rows == 0 || frame.cols == 0)
			break;

		// update the tracking result
		//! [update]
		//tracker->update(frame, box);
		Mat frame_gray;
		cv::cvtColor(frame, frame_gray, CV_BGR2GRAY);
		CMat cgray = CMat(CSize(frame_gray.cols, frame_gray.rows), MAT_8UC1, frame_gray.data);
		double time_profile_counter = (double)cvGetTickCount();
		CRect cbb = tracker.update(cgray);
		Rect bb(cbb.x, cbb.y, cbb.width, cbb.height);
		/*BBox bbox = UpdateTracker(tracker,frame_gray.data,frame_gray.cols,frame_gray.rows,frame_gray.step);
		Rect bb(bbox.x,bbox.y,bbox.w,bbox.h);*/

		time_profile_counter = (double)cvGetTickCount() - time_profile_counter;
		std::cout << "time used " << time_profile_counter / ((double)cvGetTickFrequency() * 1000) << "ms" << endl;
		//! [update]

		//! [visualization]
		// draw the tracked object
		rectangle(frame, bb, Scalar(255, 0, 0), 2, 1);

		// show image with the tracked object
		imshow("tracker", frame);
		//! [visualization]
		//quit on ESC button
		if (waitKey(1) == 27)
			break;
	}

	//ReleaseTracker(&tracker);
	return EXIT_SUCCESS;
}