//
//  HPatch.cpp
//  fypFirstDraft
//
//  Created by callen on 10/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "HPatch.h"

void HPatch::setPatchXY(boundingBox bbox){
    cout<< "bbox x,y : " << bbox.x << ", " << bbox.y << endl;
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
        rectangles.push_back(temp);
        /*cout << "xt : "<< temp.x << " yt : "<< temp.y << endl;
        cout << "f1 x, y : " << rectangles[i].x << " , " << rectangles[i].y << endl;
        cout << "f size : " << rectangles.size()<< endl;
        cout << endl;*/
    }
    //cout << "xt : "<< temp.x << " yt : "<< temp.y << endl;
}

void HPatch::setSubPatchDistance(vector<threeDPostCal> dImage){
    
    //cout << f[1].x << ", " << f[1].y << endl;
    float integral[2] = {};
    for(int n = 0; n < 2; n++) {
        for(int j = rectangles[n].y; j < rectangles[n].y + rectangles[n].h; j++){
            for(int i = rectangles[n].x; i < rectangles[n].x + rectangles[n].w; i++){
                integral[n] = integral[n] + dImage[j*640+i].d;
            }
            
        }
        integral[n] = integral[n] / (rectangles[n].w*rectangles[n].h);
    }
    subPDistance = integral[0]-integral[1];
    //cout << subPDistance << endl;
}

float HPatch::getSubPatchDistance(){
    return subPDistance;
}

float HPatch::findEuclideanDistance(ground_Truth gt){
return sqrt(pow(pC.p.x - gt.nosePosition[0],2) + pow(pC.p.y - gt.nosePosition[1],2)+pow(pC.p.d - gt.nosePosition[2],2));
    
}

void HPatch::loadGroundTruth(vector<float> gt){
    for(int i = 0; i < 6; i++){
        this->groundT.push_back(gt[i]);
    }
    
}

void PatchSet::getRandomPatches(boundingBox bbox, vector<threeDPostCal> dImage, vector<float> groundT){
	for (int i = 0; i < size; i++){
		HPatch temp(80,80);
        temp.setPatchXY(bbox);
        //cout << "xp : "<< temp.getPx() << " yp : "<< temp.getPy() << endl;
		temp.setPatchCenter(dImage);
        temp.chooseSubPatches();
        //cout << "f1 x, y : " << temp.f[1].x << " , " << temp.f[1].y << endl;
        temp.setSubPatchDistance(dImage);
        temp.loadGroundTruth(groundT);
		pSet.push_back(temp);
        this->storeGroundTruth(groundT);
        //cout << temp.getSubPatchDistance() << endl;
	}
    
	//cout << pSet[5].f[0].x <<
}

void PatchSet::storeGroundTruth(vector<float> gt){
    for(int i = 0; i < 6; i++){
        this->groundT.push_back(gt[i]);
    }
}