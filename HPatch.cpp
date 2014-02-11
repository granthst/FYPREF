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

void HPatch::chooseSubPatches(){
    sub_patch temp;
    for(int i = 0; i < 2; i++){
        temp.x = this->p_x + rand() % 80;
        temp.y = this->p_y + rand() % 80;
        temp.w = 1 + rand() % (80 - temp.x + this->p_x);
        temp.h = 1 + rand() % (80 - temp.y + this->p_y);
        f.push_back(temp);
    }
}

float HPatch::subPatchDistance(vector<threeDPostCal> dImage){
    
    float integral[2] = {};
    for(int n = 0; n < 2; n++) {
        for(int j = f[n].y; j < f[n].y + f[n].h; j++){
            for(int i = f[n].x; i < f[n].x + f[n].w; i++){
                integral[n] = integral[n] + dImage[j*640+i].d;
            }
            
        }
        integral[n] = integral[n] / (f[n].w*f[n].h);
    }
    return integral[0]-integral[1];
}

float HPatch::findEuclideanDistance(ground_Truth gt){
return sqrt(pow(pC.p.x - gt.nosePosition[0],2) + pow(pC.p.y - gt.nosePosition[1],2)+pow(pC.p.d - gt.nosePosition[2],2));
    
}