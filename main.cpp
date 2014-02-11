//
//  main.cpp
//  fypFirstDraft
//
//  Created by callen on 10/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "HPatch.h"

int g_max_z = 1300;

int n_patch = 40;
//Ultility functions

//load depth image
int16_t* loadDepthImageCompressed( const char* fname );

//load calfile
vector<float> loadCalFile(string fname);

//load trainingset
vector<int16_t*> loadTrainingSet(string depth_path);

//load groundtruthset
vector<ground_Truth> loadGroundSet(string depth_path);

//int to string
string convertInt(int number);

//load pose ground truth
void loadPoseTxt(string fname,ground_Truth& gt);

//depth to 3D
vector<threeDPostCal> get3dFromDepth(int16_t* depthImage, vector<float> depth_intrinsic);

// get bounding box around data
boundingBox getBoundingBox(int16_t* depthImage);

//generate sample patches
vector<HPatch> getRandomPatches(boundingBox bbox, int n , vector<threeDPostCal> dImage);

int main(int argc, const char * argv[])
{
    string depth_path = "/Users/Knowhow10/Desktop/01/";
	//string gt_path = "D:/kinect_head_pose_db/01/";
	string cal_filename = depth_path + "depth.cal";
    vector<float> calM;
    calM = loadCalFile(cal_filename);
    
    vector<int16_t*> trainingSet = loadTrainingSet(depth_path);
    vector<ground_Truth> groundTruthSet = loadGroundSet(depth_path);

    boundingBox bbox = getBoundingBox(trainingSet[0]);
    //cout << "bbox " << bbox.x << ", " << bbox.y << endl;
    
    vector<threeDPostCal> depthTo3D;
    depthTo3D = get3dFromDepth(trainingSet[0],calM);
    /*for(int y = 0; y < 480; y++){
        for(int x = 0; x < 640; x++){
            if(depthTo3D[y*640+x].d != 0)
				cout << depthTo3D[y*640+x].x << " , " << depthTo3D[y*640+x].y << " , " << depthTo3D[y*640+x].d << endl;
		}
	}*/
    float lambda = 0;
    vector<HPatch> rP = getRandomPatches(bbox, n_patch, depthTo3D);
    for(int i = 0; i < n_patch; i++){
        
		/*cout << " patch number " << i << " at position x, y, " << rP[i].getPx() << ", " << rP[i].getPy() << endl;
		cout << " patch center " << rP[i].getPatchCenter(depthTo3D).c << ", " << rP[i].getPatchCenter(depthTo3D).r << endl;
		cout << " 3d offset vector " << rP[i].getPatchCenter(depthTo3D).p.x << ", " << rP[i].getPatchCenter(depthTo3D).p.y << ", " << rP[i].getPatchCenter(depthTo3D).p.d << endl;*/
        
        
            rP[i].chooseSubPatches();
            lambda = rP[i].subPatchDistance(depthTo3D);
            if ( lambda != 0)
				cout << "binary test threshold candidate for patch " << i << " is : " << lambda << endl;
            //cout << "EuclideanDistance between the center of the patch and the nose tip " << rP[i].findEuclideanDistance(groundTruthSet[0]) << endl;
        
    }
    return 0;
}







//function definitions
int16_t* loadDepthImageCompressed( const char* fname ){
    
	//now read the depth image
	FILE* pFile = fopen(fname, "rb");
	if(!pFile){
		std::cerr << "could not open file " << fname << std::endl;
		return NULL;
	}
    
	int im_width = 0;
	int im_height = 0;
	bool success = true;
    
	success &= ( fread(&im_width,sizeof(int),1,pFile) == 1 ); // read width of depthmap
	success &= ( fread(&im_height,sizeof(int),1,pFile) == 1 ); // read height of depthmap
    
	int16_t* depth_img = new int16_t[im_width*im_height];
	
	int numempty;
	int numfull;
	int p = 0;
    
	while(p < im_width*im_height ){
        
		success &= ( fread( &numempty,sizeof(int),1,pFile) == 1 );
        
		for(int i = 0; i < numempty; i++)
			depth_img[ p + i ] = 0;
        
		success &= ( fread( &numfull,sizeof(int), 1, pFile) == 1 );
		success &= ( fread( &depth_img[ p + numempty ], sizeof(int16_t), numfull, pFile) == (unsigned int) numfull );
		p += numempty+numfull;
        
	}
    
	fclose(pFile);
    //cout << "depth file " << fname << " load completed " << endl;
	if(success)
		return depth_img;
	else{
		delete [] depth_img;
		return NULL;
	}
    
}

vector<float> loadCalFile(string fname){
    ifstream fInp;
	fInp.open(fname.c_str());
	if (!fInp.is_open()){
		cerr << "depth.cal file not found " << endl;
		
	}
	//read intrinsics only
    vector<float> v;
	float depth_intrinsic;
	for(int i =0; i<9; ++i){
        
		fInp >> depth_intrinsic;
		//cout << depth_intrinsic[i] << endl;
        v.push_back(depth_intrinsic);
	}
	fInp.close();
    
    return v;
}


//load trainingset
vector<int16_t*> loadTrainingSet(string depth_path){
    
    vector<int16_t*> trainingSet;
    int16_t* img = new int16_t[640*480];
    string subS = "";
    string depthFname = depth_path;
    for(int i = 3; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
		depthFname = depth_path  + "frame_00"+subS+"_depth.bin";
		
		img = loadDepthImageCompressed( depthFname.c_str());
		
		trainingSet.push_back(img);
		
	}
    return trainingSet;
}

//load groundtruthset
vector<ground_Truth> loadGroundSet(string depth_path){
    
    vector<ground_Truth> groundTruthSet;
    ground_Truth tempGt;
    string subS = "";
    string poseFname = depth_path;
    for(int i = 3; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
        poseFname = depth_path  + "frame_00"+subS+"_pose.txt";
        loadPoseTxt(poseFname,tempGt);
        groundTruthSet.push_back(tempGt);
    }
    return groundTruthSet;
}


string convertInt(int number)
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

void loadPoseTxt(string fname,ground_Truth& gt){
	ifstream fInp;
	fInp.open(fname.c_str());
	if (!fInp.is_open()){
		cerr << "pose.txt can not be found" << endl;
		
	}
	for(int i =0; i<9; ++i){
        
		fInp >> gt.rotationMatrix[i];
		//cout << depth_intrinsic[i] << endl;
	}
	for(int j =0; j<3; ++j){
        
		fInp >> gt.nosePosition[j];
		//cout << depth_intrinsic[i] << endl;
	}
}

vector<threeDPostCal> get3dFromDepth(int16_t* depthImage, vector<float> depth_intrinsic){
	
    vector<threeDPostCal> depthTo3D;
	threeDPostCal temp;
	for(int y = 0; y < 480; y++)
	{
		
        
		for(int x = 0; x < 640; x++){
            
			float d = (float)depthImage[y*640+x];
            
			if ( d < g_max_z && d > 0 ){
                
				temp.x = d * (float(x) - depth_intrinsic[2])/depth_intrinsic[0];
				temp.y = d * (float(y) - depth_intrinsic[5])/depth_intrinsic[4];
				temp.d = d;
				//cout << img3Di[x][2] << endl;
			}
			else{
                
				temp.x = 0;
				temp.y = 0;
				temp.d = 0;
			}
			depthTo3D.push_back(temp);
		}
	}
	return depthTo3D;
}

vector<HPatch> getRandomPatches(boundingBox bbox, int n , vector<threeDPostCal> dImage){
	
	vector<HPatch> rP;
    HPatch temp(80,80);
	for (int i = 0; i < n; i++){
		//cout << "xp : "<< temp.x << " yp : "<< temp.y << endl;
        temp.setPatchXY(bbox);
		temp.setPatchCenter(dImage);
		rP.push_back(temp);
	}
	return rP;
}

boundingBox getBoundingBox(int16_t* depthImage){
    
	boundingBox bbox;
	int min_x = 640;
	int min_y = 480;
	int max_x = 0;
	int max_y = 0;
    
    
	for(int y = 0; y < 480; y++)
	{
	    
	    for(int x = 0; x < 640; x++){
            
	    	if( depthImage[y*640+x] > 0) {
                
				min_x = min(min_x,x); min_y = min(min_y,y);
				max_x = max(max_x,x); max_y = max(max_y,y);
			}
            
	    }
	}
	bbox.x = min_x;
	bbox.y = min_y;
	bbox.width = max_x - min_x;
	bbox.height = max_y - min_y;
    
	return bbox;
    
}