//
//  HPatch.h
//  fypFirstDraft
//
//  Created by callen on 10/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#ifndef __fypFirstDraft__HPatch__
#define __fypFirstDraft__HPatch__

#include <iostream>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace std;
using namespace cv;

#endif /* defined(__fypFirstDraft__HPatch__) */

struct ground_Truth{
	float rotationMatrix[9];
	float nosePosition[3];
};

struct threeDPostCal{
	float x,y,d;
};

struct patchCenter{
	int r,c;
	threeDPostCal p;
    
};

struct boundingBox{
	int x,y;
	int width,height;
};

struct sub_patch{
    int x,y,w,h;
};



class HPatch{
public:
    
    HPatch(int w, int h){
        p_width = w;
        p_height = h;
    }
    ~HPatch(){}
    void setPatchXY(boundingBox bbox);
    int getPx();
    int getPy();
    void setPatchCenter(Mat dImage);
    patchCenter getPatchCenter(vector<threeDPostCal> dImage);
    void chooseSubPatches(vector<sub_patch> SP,Mat itegralImage);
    void loadSubPatches(const string fname);
    void setSubPatchDistance(vector<threeDPostCal> dImage);
    float getSubPatchDistance();
    float findEuclideanDistance(vector<float> gt);
    void loadGroundTruth(vector<float> groundT);
    vector<sub_patch> rectangles;
    vector<float> groundT;
    //int groundTruthIndex;

    int index;
    int p_width, p_height, p_x, p_y;
    //vector<sub_patch> f;
    patchCenter pC;
    float subPDistance;
    bool positive;
    void printRectangles();
};

class PatchSet{
    
public:
    PatchSet(){}
    PatchSet(int n){size = n;}
    ~PatchSet(){}
    void getRandomPatches(boundingBox bbox , Mat dImage, vector<float> groundT ,bool p);
    void storeGroundTruth(vector<float> groundT);
    vector<HPatch> pSet;
    vector<float> groundT;
    int size;
};