//
//  HPatch.cpp
//  fypFirstDraft
//
//  Created by callen on 10/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "HPatch.h"
#define SIZE_X 10
#define SIZE_Y 10
void HPatch::setPatchXY(boundingBox bbox){
    //cout<< "bbox x,y : " << bbox.x << ", " << bbox.y << endl;
    p_x = bbox.x + rand() % bbox.width;
    p_y = bbox.y + rand() % bbox.height;
}

int HPatch::getPx(){
    return p_x;
}

int HPatch::getPy(){
    return p_y;
}

void HPatch::setPatchCenter(Mat dImage){
    //int centerP = (this->p_y+40 ) * 640 + (this->p_x+40);
	//cout << dImage.size() << endl;
	//cout << centerP <<endl;
	pC.c = this->p_x+40;
	pC.r = this->p_y+40;
	pC.p.x = dImage.at<Vec3f>(pC.r,pC.c)[0];
	pC.p.y = dImage.at<Vec3f>(pC.r,pC.c)[1];
	pC.p.d = dImage.at<Vec3f>(pC.r,pC.c)[2];
}

patchCenter HPatch::getPatchCenter(vector<threeDPostCal> dImage){
    return this->pC;
}

void HPatch::chooseSubPatches(vector<sub_patch> SP, Mat itegralImage){
        /*cout << "xt : "<< temp.x << " yt : "<< temp.y << endl;
        cout << "f1 x, y : " << rectangles[i].x << " , " << rectangles[i].y << endl;
        cout << "f size : " << rectangles.size()<< endl;
        cout << endl;*/
    sub_patch temp;
    //cout << SP[0].x << endl;
    for(int i = 0; i < 2; i++){
        temp.y = SP[i].x + this->p_x;
        temp.x = SP[i].y + this->p_y;
        temp.h = SP[i].w ;
        temp.w = SP[i].h ;
        rectangles.push_back(temp);
        //cout << "sub patch x " << rectangles[i].x << endl;
    }
//    cout << " check violation x1: " << rectangles[0].x+rectangles[0].w << " y: " << rectangles[0].y+rectangles[0].h<< endl;
//    cout << " check violation x2: " << rectangles[0].x << " y: " << rectangles[0].y<< endl;
//    cout << " check violation x3: " << rectangles[0].x+rectangles[0].w << " y: " << rectangles[0].y<< endl;
//    cout << " check violation x4: " << rectangles[0].x << " y: " << rectangles[0].y+rectangles[0].h<< endl;
    double sum1 = (itegralImage.at<int>(rectangles[0].x+rectangles[0].w,rectangles[0].y+rectangles[0].h) + itegralImage.at<int>(rectangles[0].x,rectangles[0].y) - itegralImage.at<int>(rectangles[0].x+rectangles[0].w,rectangles[0].y) - itegralImage.at<int>(rectangles[0].x,rectangles[0].y+rectangles[0].h))/(double)(rectangles[0].w*rectangles[0].h);
    //cout << sum1 << endl;
    double sum2 = (itegralImage.at<int>(rectangles[1].x+rectangles[1].w,rectangles[1].y+rectangles[1].h) + itegralImage.at<int>(rectangles[1].x,rectangles[1].y) - itegralImage.at<int>(rectangles[1].x+rectangles[1].w,rectangles[1].y) - itegralImage.at<int>(rectangles[1].x,rectangles[1].y+rectangles[1].h))/(double)(rectangles[1].w*rectangles[1].h);
    //cout << "xt : "<< temp.w << " yt : "<< temp.h << endl;
    subPDistance = sum1 - sum2;
    rectangles.clear();
}


void HPatch::loadSubPatches(const string fname){
    sub_patch temp;
    int dummy;
    ifstream fInp;
    fInp.open(fname);
    for(int i = 0; i < 2; i++){
        fInp >> temp.x;
        fInp >> temp.y;
        fInp >> temp.w;
        fInp >> temp.h;
        rectangles.push_back(temp);
        // cout << "sub patch x " << rectangles[i].x << endl;
    }
    fInp >> subPDistance;
    fInp >> dummy;
}

void HPatch::printRectangles() {
    int noRectangles = (int) rectangles.size();
    cout << "noRectangles = " << noRectangles << endl;
    for (int i = 0; noRectangles; i++) {
        cout << "rectangle: " << i << " [x, y, w, h] = [" << rectangles[i].x << ", " << rectangles[i].y << ", " << rectangles[i].w << ", " << rectangles[i].h << "]" << endl;
    }
}

void HPatch::setSubPatchDistance(vector<threeDPostCal> dImage){
    //cout << f[1].x << ", " << f[1].y << endl;
    float integral[2] = {};
    for(int n = 0; n < 2; n++) {
        //cout <<  rectangles[n].x << endl;
        for(int j = rectangles[n].y; j < rectangles[n].y + rectangles[n].h; j++){
            for(int i = rectangles[n].x; i < rectangles[n].x + rectangles[n].w; i++){
                integral[n] = integral[n] + dImage[(j*640)+i].d;
            }
            
        }
        integral[n] = integral[n] / (rectangles[n].w * rectangles[n].h);
    }
    //printRectangles();
    subPDistance = integral[0]-integral[1];
    rectangles.clear();
    //cout << subPDistance << endl;
}

float HPatch::getSubPatchDistance(){
    return subPDistance;
}

float HPatch::findEuclideanDistance(vector<float> gt){

    return sqrt(pow(pC.p.x - gt[0],2) + pow(pC.p.y - gt[1],2)+pow(pC.p.d - gt[2],2));
    
}

void HPatch::loadGroundTruth(vector<float> gt){
    for(int i = 0; i < 6; i++){
        this->groundT.push_back(gt[i]);
    }
    
    
    
}

void HPatch::setTestPatch(boundingBox bbox,int stepsizeX, int stepsizeY){
    p_x = bbox.x + stepsizeX;
    p_y = bbox.y + stepsizeY;
}

void PatchSet::getRandomPatches(boundingBox bbox, Mat dImage, vector<float> groundT, bool p){
	//for (int i = 0; i < size; i++){
    for(int i = 0; i < size/2; i++){
		HPatch temp(80,80);
        temp.setPatchXY(bbox);
        //cout << "bbox " << bbox.x <<  " " << bbox.y << " " << bbox.width << " " << bbox.height << endl;
        //cout << "xp : "<< temp.getPx() << " yp : "<< temp.getPy() << endl;
		temp.setPatchCenter(dImage);
        //temp.chooseSubPatches();
        //cout << "f1 x, y : " << temp.f[1].x << " , " << temp.f[1].y << endl;
        //temp.setSubPatchDistance(dImage);
        //cout << temp.pC.r << " " << temp.pC.c << endl;
        //cout << groundT[0] << " " << groundT[1] << " " << groundT[2] << endl;
        //vector<float> nearestPoint;
        //cout << temp.findEuclideanDistance(groundT) << endl;
        
        temp.loadGroundTruth(groundT);
        //cout << temp.pC.c*640+temp.pC.r << endl;
        temp.groundT[0] = temp.groundT[0] - temp.pC.p.x;
        temp.groundT[1] = temp.groundT[1] - temp.pC.p.y;
        temp.groundT[2] = temp.groundT[2] - temp.pC.p.d;
        temp.positive = p;
        pSet.push_back(temp);
        this->storeGroundTruth(groundT);
        
        //cout << temp.getSubPatchDistance() << endl;
	}
    
	//cout << pSet[5].f[0].x <<
}

void PatchSet::getRandomPatches2d(boundingBox bbox, Mat dImage, vector<vector<float>> groundT2d, bool p){
    int np = 30;
    if( p == 1)
        np = 10;
    for(int i = 0; i < np; i++){
        HPatch temp(80,80);
        temp.setPatchXY(bbox);
        temp.pC.c = temp.p_x + 40;
        temp.pC.r = temp.p_y + 40;
        //cout << groundT2d.size() <<endl;
        for(int j = 0; j < groundT2d.size(); j ++){
            temp.groundT2d.push_back(groundT2d[j]);
        }
        //if(p == 1)
        //cout << "p center x y " << (float)temp.pC.c/640 << " " << (float)temp.pC.r/480 << endl;
        for(int k = 0; k < 17; k++){
            
            temp.groundT2d[k][0] = temp.groundT2d[k][0] - (float)temp.pC.c/640;
            temp.groundT2d[k][1] = temp.groundT2d[k][1] - (float)temp.pC.r/480;
            //if(p == 1)
            //cout << "offset x,y for feature "<< k << " : " << temp.groundT2d[k][0] << " " << temp.groundT2d[k][1]<< endl;
        }
        temp.positive = p;
        pSet.push_back(temp);
    }
}

void PatchSet::sampleTestPatches(boundingBox bbox,Mat dImage){
    //cout << bbox.width << endl;
    //cout << bbox.height << endl;
    for(int x = 0; x < bbox.width-80; x = x + SIZE_X){
        for(int y = 0; y < bbox.height-80; y = y + SIZE_Y){
            HPatch temp(80,80);
            temp.setTestPatch(bbox, x, y);
            temp.setPatchCenter(dImage);
            
            //cout << "x " << temp.p_x << endl;
            //cout << "y " << temp.p_y << endl;
            //cout << "x coor " << temp.pC.p.x << endl;
            pSet.push_back(temp);
        }
    }
}

void PatchSet::sampleTestPatches2d(boundingBox bbox,Mat dImage){
    //cout << bbox.width << endl;
    //cout << bbox.height << endl;
    for(int x = 0; x < bbox.width-80; x = x + SIZE_X){
        for(int y = 0; y < bbox.height-80; y = y + SIZE_Y){
            HPatch temp(80,80);
            temp.setTestPatch(bbox, x, y);
            //temp.setPatchCenter(dImage);
            temp.pC.c = temp.p_x + 40;
            temp.pC.r = temp.p_y + 40;
//            cout << "x " << temp.p_x << endl;
//            cout << "y " << temp.p_y << endl;
            //cout << "x coor " << temp.pC.p.x << endl;
            pSet.push_back(temp);
        }
    }
}

void PatchSet::storeGroundTruth(vector<float> gt){
    for(int i = 0; i < 6; i++){
        this->groundT.push_back(gt[i]);
    }
}