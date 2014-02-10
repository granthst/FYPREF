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
using namespace std;


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
    void setPatchCenter(vector<threeDPostCal> dImage);
    patchCenter getPatchCenter(vector<threeDPostCal> dImage);
    float integralImage(vector<threeDPostCal> dImage);
    float findEuclideanDistance(ground_Truth gt);
    
private:
    
    int p_width, p_height, p_x, p_y;
    patchCenter pC;
};