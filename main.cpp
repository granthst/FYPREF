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
#include <ctime>
#include "DTree.h"

#define PATCH_SET_SIZE 40
#define DATA_SET_SIZE 10
#define n 6
#define THESHOLD_MAX 1000
using namespace cv;
string preProcessedDataFile = "randomPatches.txt";
string treefile = "tree.txt";
bool swaped=false;

int g_max_z = 1300;

int n_patch = 40;
//Ultility functions

//load depth image
int16_t* loadDepthImageCompressed( const char* fname );

//load calfile
vector<float> loadCalFile(string fname);

//load pose bin
vector<float> loadPoseBin(string fname);

//load trainingset
vector<int16_t*> loadTrainingSet(string depth_path);

//load groundtruthset txt
vector<ground_Truth> loadGroundSet(string depth_path);

//load groundtruthset bin
vector<vector<float>> loadGroundSetBin(string bin_path);

//int to string
string convertInt(int number);

//load pose ground truth
void loadPoseTxt(string fname,ground_Truth& gt);

//depth to 3D
vector<threeDPostCal> get3dFromDepth(int16_t* depthImage, vector<float> depth_intrinsic);

// get bounding box around data
boundingBox getBoundingBox(int16_t* depthImage);

//generate sample patches
vector<HPatch> getRandomPatches(boundingBox bbox, int n1 , vector<threeDPostCal> dImage);

void generateRandomSubPatches(vector<vector<sub_patch>>& rectangles);

void generateRandomThreshold(vector<int>& rT);

void writePreProcessedData(vector<HPatch> PS, vector<vector<sub_patch>>& rectangles, vector<int>& rT, string fname);
//compute the meanVector
vector<float> computeMeanVector(const vector<vector<float>>& groundTruthSet);

//compute the cv
vector<vector<float>> computeCovariance(const vector<vector<float>>& groundTruthSet);


void getIntegralImageDepthChannel(Mat& depthI,int16_t* trainingSet,Mat* channels,vector<float> calM,Mat& img3D);


int main(int argc, const char * argv[])
{
    //srand(time(NULL));
    
    string depth_path = "/Users/Knowhow10/Desktop/01/";
    string bin_path = "/Users/Knowhow10/Desktop/db_annotations/01/";
	string cal_filename = depth_path + "depth.cal";
    vector<float> calM;
    calM = loadCalFile(cal_filename);
    /***vector<threeDPostCal> depthTo3D,integralImage;
    vector<int16_t*> trainingSet = loadTrainingSet(depth_path);
    vector<ground_Truth> groundTruthSet = loadGroundSet(depth_path);
    vector<vector<float>> groundTruthSetBin = loadGroundSetBin(bin_path);
    boundingBox bbox;
    depthTo3D = get3dFromDepth(trainingSet[0],calM);
    vector<Mat> depthIntegral;
    //vector<vector<Mat>> channelsList;
    vector<Mat> img3DList;
    for(int i = 0; i < trainingSet.size(); i++){
        
        Mat* channels = new Mat[3];
        Mat depthI = Mat(480, 640, CV_32FC3);
        Mat sumI = Mat(481, 641, CV_64FC3);
        Mat img3D;
        getIntegralImageDepthChannel(depthI,trainingSet[i],channels,calM,img3D);
        
        integral(channels[2],sumI);
        depthIntegral.push_back(sumI);
        img3DList.push_back(img3D);
        cout << "image " << i << " done" << endl;
    }


    
    
    cout << trainingSet.size() << endl;
    
    //vector<PatchSet> pSList;
    vector<HPatch> wholeDataSet;
    vector<vector<threeDPostCal>> depthTo3DSet;
    for( int i = 0; i < DATA_SET_SIZE; i++ ){
        cout << "#####     i = " << i << endl;
        PatchSet pS(PATCH_SET_SIZE);
        bbox = getBoundingBox(trainingSet[i]);
        //depthTo3D = get3dFromDepth(trainingSet[i],calM);
        //depthTo3DSet.push_back(depthTo3D);
        pS.getRandomPatches(bbox, img3DList[i], groundTruthSetBin[i]);
        cout << " gt : " << groundTruthSetBin[i][0] << " " << groundTruthSetBin[i][1] << " " << groundTruthSetBin[i][2] << " " << groundTruthSetBin[i][3] << " " << groundTruthSetBin[i][4] << " " << groundTruthSetBin[i][5] << endl;
        cout << " meanVector : " << pS.pSet[0].groundT[0] << " " << pS.pSet[0].groundT[1] << " " << pS.pSet[0].groundT[2] << " " << pS.pSet[0].groundT[3] << " " << pS.pSet[0].groundT[4] << " " << pS.pSet[0].groundT[5] << endl;
        //pSList.push_back(pS);
        for (int j = 0; j <  PATCH_SET_SIZE; j++){
            pS.pSet[j].index = i;
            wholeDataSet.push_back(pS.pSet[j]);
        }
    }
    
    
    //cout << "depthTo3D size : " << depthTo3D.size() << endl;
    //cout << "depthTo3Dset size : " << depthTo3DSet.size() << endl;
    

    
    /*string outputTreefile = "tree1.txt";
    DTree tree2(10);

    tree2.growTree(wholeDataSet, depthIntegral);
    tree2.write_tree(outputTreefile);*/
    
    
    
    DTree testTree;
    testTree.read_tree(treefile);
    //cout << testTree.noNodes << endl;
    //cout << testTree.treeTable[5].bestF[0].x << endl;
    
    string testImagefile = "/Users/Knowhow10/Desktop/01/frame_00025_depth.bin";
    string testImageGTfile = "/Users/Knowhow10/Desktop/db_annotations/01/frame_00025_depth.bin";
    int16_t* testImage = loadDepthImageCompressed(testImagefile.c_str());
    boundingBox testBbox = getBoundingBox(testImage);
    vector<vector<float>> groundTruthSetBin = loadGroundSetBin(bin_path);
    Mat* channels = new Mat[3];
    Mat depthI = Mat(480, 640, CV_32FC3);
    Mat sumI = Mat(481, 641, CV_64FC3);
    Mat img3D;
    getIntegralImageDepthChannel(depthI,testImage,channels,calM,img3D);
    integral(channels[2],sumI);
    vector<float> testGt = groundTruthSetBin[21];
        cout << "test " << testGt[0] << " " << testGt[1] << " " << testGt[2] << " " << testGt[3] << " " << testGt[4] << " " << testGt[5] << endl;
    vector<vector<float>> estimatedMean;
    testTree.regressionEstimation(sumI, testBbox, testGt, estimatedMean,img3D);
    cout << " size mean " << estimatedMean.size() << endl;
    vector<float> estimatedMean2 = computeMeanVector(estimatedMean);
    
        
    cout << "estimated " << estimatedMean2[0] << " " << estimatedMean2[1] << " " << estimatedMean2[2] << " " << estimatedMean2[3] << " " << estimatedMean2[4] << " " << estimatedMean2[5] << endl;
    cout << "test " << testGt[0] << " " << testGt[1] << " " << testGt[2] << " " << testGt[3] << " " << testGt[4] << " " << testGt[5] << endl;
    //PatchSet testPS(PATCH_SET_SIZE);
    //testPS.getRandomPatches(testBbox, test3D, testGt);
    //testPS.pSet[0].rectangles = testTree.treeTable[0].bestF;
    //cout << testPS.pSet[0].rectangles[0].x << endl;
    
    
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
    for(int i = 4; i < 500; i++){
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
    for(int i = 4; i < 500; i++){
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

vector<vector<float>> loadGroundSetBin(string bin_path){
    vector<vector<float>> groundSetBin;
    vector<float> tempG;
    string subS = "";
    string poseFname = bin_path;
    for(int i = 4; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
        poseFname = bin_path  + "frame_00"+subS+"_pose.bin";
        tempG = loadPoseBin(poseFname);
        //cout << tempG[0] << endl;
        groundSetBin.push_back(tempG);
    }
    
    
    return groundSetBin;
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

vector<float> loadPoseBin(string fname){
    bool have_gt = false;
    vector<float> gt1;
	float *gt = new float[6] ;
    FILE* pFile = fopen(fname.c_str(), "rb");
	if(pFile){
        //cout << " bin file read correctly " << endl;
		have_gt = true;
		have_gt &= ( fread( &gt[0], sizeof(float),6, pFile) == 6 );
		fclose(pFile);
	}
    for (int i = 0; i < 6; i ++)
        gt1.push_back(gt[i]);
   // cout << gt[0] << " " << gt[1] << " " << gt[2] << endl;
    return gt1;
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

/*vector<HPatch> getRandomPatches(boundingBox bbox, int n1 , vector<threeDPostCal> dImage){
	
	vector<HPatch> rP;
    HPatch temp(80,80);
	for (int i = 0; i < n; i++){
		//cout << "xp : "<< temp.x << " yp : "<< temp.y << endl;
        temp.setPatchXY(bbox);
		temp.setPatchCenter(dImage);
		rP.push_back(temp);
	}
	return rP;
}*/

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

void generateRandomSubPatches(vector<vector<sub_patch>>& rectangles){
    sub_patch temp;
    vector<sub_patch> tempSub;
    for (int j = 0; j < PATCH_SET_SIZE * DATA_SET_SIZE; j++){
        for(int i = 0; i < 2; i++){
            temp.x = rand() % 80;
            temp.y = rand() % 80;
            temp.w = 1 + rand() % (80 - temp.x);
            temp.h = 1 + rand() % (80 - temp.y);
            tempSub.push_back(temp);
            //cout << tempSub[i].x << endl;
        }
        rectangles.push_back(tempSub);
        //cout << " rec " << rectangles[j][0].x << endl;
        tempSub.clear();
        //cout << " rec " << rectangles[j][0].x << endl;
    }
}

void generateRandomThreshold(vector<int>& rT){
    
    
    float rt = 0;
    for (int i = 0; i < PATCH_SET_SIZE * DATA_SET_SIZE; i++){
        rt = -THESHOLD_MAX + (rand() %(2*THESHOLD_MAX));
        rT.push_back((int) rt);
    }
    
}

/*void writePreProcessedData(vector<HPatch> PS, vector<vector<sub_patch>>& rectangles, vector<int>& rT, string fname){
    generateRandomSubPatches(rectangles);
    generateRandomThreshold(rT);
    ofstream fOut;
    fOut.open(fname);
    for(int j = 0; j < rT.size(); j++){
        for( int i = 0; i < PS.size(); i++){
            
            PS[i].chooseSubPatches(rectangles[j]);
            //cout << " f1 " << i << "x " << PS[i].rectangles[0].x << endl;
            fOut << PS[i].p_x << " " << PS[i].p_y << endl;
            fOut << rectangles[j][0].x << " " << rectangles[j][0].y << " " << rectangles[j][0].w << " " << rectangles[j][0].h << endl;
            fOut << rectangles[j][1].x << " " << rectangles[j][1].y << " " << rectangles[j][1].w << " " << rectangles[j][1].h << endl;
            fOut  << rT[j] << endl;
        }
        cout << " rt " << j << " " << rT[j]<< " done" << endl;
    }
    fOut.close();
    
}*/


vector<float> computeMeanVector(const vector<vector<float>>& groundTruthSet){
    vector<float> mean;
    for(int i = 0; i < 6; i++){
        float temp = 0.0;
        for(int j = 0; j < groundTruthSet.size(); j++){
            temp = groundTruthSet[j][i] + temp;
        }
        temp = temp/groundTruthSet.size();
        mean.push_back(temp);
    }
    return mean;
}

vector<vector<float>> computeCovariance(const vector<vector<float>>& groundTruthSet){
    vector<vector<float>> Covariance;
    vector<float> meanVector = computeMeanVector(groundTruthSet);
    for(int i = 0; i < 6; i++){
        vector<float> col;
        for(int j = 0; j < 6; j++){
            float temp = 0.0;
            for(int k = 0; k < groundTruthSet.size(); k++){
                temp = temp + groundTruthSet[k][i]*groundTruthSet[k][j];
            }
            float cov = temp - (meanVector[i]*meanVector[j]);
            col.push_back(cov);
        }
        Covariance.push_back(col);
    }
    return Covariance;
}



void getIntegralImageDepthChannel(Mat& depthI,int16_t* trainingSet,Mat* channels,vector<float> calM,Mat& img3D){
    
    
    depthI = Mat(480, 640, CV_32FC3);
    img3D.create( depthI.rows, depthI.cols, CV_32FC3 );
    for(int row = 0; row < 480; row++){
        for(int col = 0; col < 640; col++){
            depthI.at<int16_t>(row, col) = trainingSet[640*row + col];
            //cout << depthI.at<int16_t>(row,col) << endl;
            
        }
    }
    
    
    
    for(int y = 0; y < img3D.rows; y++)
	{
		Vec3f* img3Di = img3D.ptr<Vec3f>(y);
        //cout << *img3Di << endl;
        //cout << *img3D.ptr<Vec3f>(y) << endl;
		const int16_t* depthImgi = depthI.ptr<int16_t>(y);
        //cout << *depthImgi << endl;
		for(int x = 0; x < img3D.cols; x++){
            
			float d = (float)depthImgi[x];
            //cout << depthImgi[x] << endl;
			if ( d < g_max_z && d > 0 ){
                
				img3Di[x][0] = d * (float(x) - calM[2])/calM[0];
				img3Di[x][1] = d * (float(y) - calM[5])/calM[4];
				img3Di[x][2] = d;
                
			}
			else{
                
				img3Di[x] = 0;
			}
            //cout << x << endl;
            //cout << img3Di[x] << " after cal " << endl;
		}
	}
    
    
	split(img3D, channels);
}
