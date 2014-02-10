//
//  HPatch.cpp
//  fypFirstDraft
//
//  Created by callen on 10/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "HPatch.h"

void HPatch::setPatchXY(boundingBox bbox){
    p_x = bbox.x + rand() % bbox.width - 80;
    p_y = bbox.y + rand() % bbox.height - 80 ;
}

int HPatch::getPx(){
    return p_x;
}

int HPatch::getPy(){
    return p_y;
}

void HPatch::setPatchCenter(vector<threeDPostCal> dImage){
    int centerP = (this->p_y+40 ) * 640 + (this->p_x+40);
	//cout << dImage.size() << endl;
	//cout << centerP <<endl;
	pC.c = this->p_x+40;
	pC.r = this->p_y+40;
	pC.p.x = dImage[centerP].x;
	pC.p.y = dImage[centerP].y;
	pC.p.d = dImage[centerP].d;
}

patchCenter HPatch::getPatchCenter(vector<threeDPostCal> dImage){
    return this->pC;
}

float HPatch::integralImage(vector<threeDPostCal> dImage){
    int x,y,w,h;
    //ofstream fOut;
    //fOut.open("subpatch.txt");
    float integral = 0;
    x = this->p_x + rand() % 80;
    y = this->p_y + rand() % 80;
    w = 1 + rand() % (80 - x + this->p_x);
    h = 1 + rand() % (80 - y + this->p_y);
    //cout << "rP.x : " << rP.x << endl;
    //cout << "rP.y : " << rP.y << endl;
    //cout << "x : " << x << endl;
    //cout << "y : " << y << endl;
    //cout << "w : " << w << endl;
    //cout << "h : " << h << endl;
    //cout << y + h << endl;
    //cout << x + w << endl;
    //cout << dImage.size() << endl;
    for(int j = y; j < y + h; j++){
        for(int i = x; i < x + w; i++){
            integral = integral + dImage[j*640+i].d;
            //cout << i << ", " << j << endl;
            //fOut << dImage[j*640+i].d <<  " " ;
        }
        //fOut << endl;
    }
    
    integral = integral / (w*h);
    //cout << integral << endl;
    //fOut.close();
    return integral;
}

float HPatch::findEuclideanDistance(ground_Truth gt){
return sqrt(pow(pC.p.x - gt.nosePosition[0],2) + pow(pC.p.y - gt.nosePosition[1],2)+pow(pC.p.d - gt.nosePosition[2],2));
    
}