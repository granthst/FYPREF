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
#include <opencv2/objdetect/objdetect.hpp>
#include "DForest.h"

#define PATCH_SET_SIZE 80
#define DATA_SET_SIZE 496
#define NO_TEST_SUBJECTS 100
#define n 6
#define THESHOLD_MAX 800
#define NUMBER_TREES 10
#define NO_MEAN 1
#define ERROR_TH 80
#define ERROR_TH_2D 0.015
#define NUM_TEST 500
#define NUM_LANDMARKS 17
#define NUMBER_POSE 19
using namespace cv;

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
void loadTrainingSet(string depth_path,vector<int16_t*>& trainingSet);

//load groundtruthset txt
vector<ground_Truth> loadGroundSet(string depth_path);

//load groundtruthset bin
void loadGroundSetBin(string bin_path,vector<vector<float>>& groundSetBin );

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

vector<double> Error(vector<float> estimated, vector<float> gt);

void loadRGBImages(vector<Mat>& rgb, string filepath);

void loadMaskImages(vector<Mat>& rgb, string filepath);

void saveBBox(vector<boundingBox> masks, string fname);

void loadBBox(vector<boundingBox>& bboxSet, int size);

void findMasks(vector<Mat> masks, vector<boundingBox>& headMasks);

boundingBox findMask(Mat mask);

void eliminateOutliars(vector<vector<float>> estimated,vector<float>& newMean);

void runTestImage(int i,DForest forest,vector<float> calM,bool& noseA,bool& poseA);

vector<vector<float>> readLandMarks(string fname);

vector<Point> getdots(vector<vector<float>> landmarks);

void displydots(vector<Point> dots,Mat image, bool t);

void loadPoses( const string& fileName,vector<vector<int>>& mDepthList );

vector<int> loadLandmarks(string fname);

boundingBox findBbox2d(vector<Point> dots);

void load2dTrainingset(vector<Mat>& trainingset2d, vector<vector<vector<float>>>& landmarks, vector<boundingBox>& bboxs);

void pickWantedLandmarks(vector<vector<float>> landmarks, vector<vector<float>>& newlandmarks);

vector<vector<float>> computeMeanVector2d(vector<vector<float>> data);

void removeoutliars(vector<vector<float>>& newMean,vector<vector<float>> data);

void writeAccuracyvsNoTree(int trees,vector<float> noseAccuracy, vector<float> poseAccuracy);

void writeAngleError(vector<float> estimated, vector<float> ground, int angle);

vector<float> estimatedPitch,estimatedYaw,estimatedRoll,groundPitch,groundYaw,groundRoll;

void detectAndDisplay( Mat frame, vector<Rect>& faces );

/** Global variables */
String face_cascade_name = "haarcascade_frontalface_alt.xml";
//String eyes_cascade_name = "haarcascade_eye_tree_eyeglasses.xml";
CascadeClassifier face_cascade;
//CascadeClassifier eyes_cascade;
string window_name = "Capture - Face detection";
RNG rng(12345);
int main(int argc, const char * argv[])
{
    //srand(time(NULL));
    

    
//    vector<threeDPostCal> integralImage;
//    vector<int16_t*> trainingSet;
//    vector<vector<float>> groundTruthSetBin;
//    vector<HPatch> wholeDataSet;
//    vector<Mat> img3DList;
//    vector<Mat> depthIntegral;
//    vector<Mat> rgbs;
//    vector<Mat> rgbIntegral;
//    //read the images
//    for(int j = 1; j <= NO_TEST_SUBJECTS; j++){
//        string depth_path;
//        string bin_path;
//        if(j!= 6){
//            if(j < 10){
//                depth_path = "/Users/Knowhow10/Downloads/kinect_head_pose_db/0"+convertInt(j)+"/";
//                bin_path = "/Users/Knowhow10/Desktop/db_annotations/0"+convertInt(j)+"/";
//            }
//            else{
//                depth_path = "/Users/Knowhow10/Downloads/kinect_head_pose_db/"+convertInt(j)+"/";
//                bin_path = "/Users/Knowhow10/Desktop/db_annotations/"+convertInt(j)+"/";
//            }
//            //loadTrainingSet(depth_path,trainingSet);
//            loadRGBImages(rgbs, depth_path);
//            loadGroundSetBin(bin_path,groundTruthSetBin);
//            cout << rgbs.size() << endl;
//            cout << groundTruthSetBin.size() << endl;
//        }
//    }
//    for(int i = 0; i < rgbs.size(); i++){
//        Mat Sumimage = Mat(481,641,CV_64FC3);
//        integral(rgbs[i], Sumimage);
//        rgbIntegral.push_back(Sumimage);
//    }
    //cout << rgbIntegral[0] << endl;
//    for(int i = 1; i <= NO_TEST_SUBJECTS; i++){
//        if(i != 6){
//            cout << i << "th person " << endl;
//            string depth_path;
//            if(i < 10){
//                depth_path = "/Users/Knowhow10/Downloads/kinect_head_pose_db/0"+convertInt(i)+"/";
//            }
//            else{
//                depth_path = "/Users/Knowhow10/Downloads/kinect_head_pose_db/"+convertInt(i)+"/";
//            }
//            string cal_filename = depth_path + "depth.cal";
//            vector<float> calM;
//            calM = loadCalFile(cal_filename);
//            for(int k = 0; k < DATA_SET_SIZE; k++ ){
//                //cout << k << "th image " <<endl;
//                Mat* channels = new Mat[3];
//                Mat depthI = Mat(480, 640, CV_32FC3);
//                Mat sumI = Mat(481, 641, CV_64FC3);
//                Mat img3D;
//                int offset = 0;
//                if(i<7)
//                    offset = (i-1)*DATA_SET_SIZE+k;
//                else
//                    offset = (i-2)*DATA_SET_SIZE+k;
//                getIntegralImageDepthChannel(depthI,trainingSet[offset],channels,calM,img3D);
//                integral(channels[2],sumI);
//                depthIntegral.push_back(sumI);
//                
//                img3DList.push_back(img3D);
//            }
//        }
//        cout << depthIntegral.size() << endl;
//    }
    //cout << depthIntegral.size() << endl;
//    vector<boundingBox> headmasks;
//    loadBBox(headmasks, NO_TEST_SUBJECTS);
//    //generating the patches
//    boundingBox bbox;
//    for( int i = 0; i < trainingSet.size(); i++ ){
//        cout << "#####     i = " << i << endl;
//        PatchSet pS(PATCH_SET_SIZE);
//        bbox = getBoundingBox(trainingSet[i]);
//        //cout << "bbox big" << bbox.x <<  " " << bbox.y << endl;
//        //depthTo3D = get3dFromDepth(trainingSet[i],calM);
//        pS.getRandomPatches(headmasks[i], img3DList[i], groundTruthSetBin[i],0);
//        //cout << "number of positive patches" << pS.pSet.size() << endl;
//        pS.getRandomPatches(bbox, img3DList[i], groundTruthSetBin[i],1);
//        //cout << " gt : " << groundTruthSetBin[i][0] << " " << groundTruthSetBin[i][1] << " " << groundTruthSetBin[i][2] << " " << groundTruthSetBin[i][3] << " " << groundTruthSetBin[i][4] << " " << groundTruthSetBin[i][5] << endl;
//        //cout << "number of patches" << pS.pSet.size() << endl;;
//        for (int j = 0; j <  PATCH_SET_SIZE; j++){
//            pS.pSet[j].index = i;
//            wholeDataSet.push_back(pS.pSet[j]);
//        }
//        cout << "number of patches " << wholeDataSet.size() << endl;
//    }
    
//    trainingSet.clear();
//    groundTruthSetBin.clear();
//    //build the forest
//    DForest forest(NUMBER_TREES);
//    forest.growForest(wholeDataSet, depthIntegral);
//    forest.writeForest();
    
    
    //DTree testTree;
    //testTree.read_tree(treefile);
//    vector<float> noseA,poseA;
//    
//    DForest forest(NUMBER_TREES);
//    forest.loadTree();
//    string cal_filename = "/Users/Knowhow10/Downloads/kinect_head_pose_db/09/depth.cal";
//    vector<float> calM;
//    calM = loadCalFile(cal_filename);
//    
//    bool nose,pose;
//    int a = 0;
//    int b = 0;
//    for(int j = 4; j < NUM_TEST; j++){
//        nose = 0;
//        pose = 0;
//        runTestImage(j, forest, calM,nose,pose);
//        if(nose ==1)
//            a++;
//        if(pose == 1)
//            b++;
//    }
//    cout << estimatedPitch.size() << endl;
//    writeAngleError(estimatedPitch, groundPitch, 0);
//    writeAngleError(estimatedYaw, groundYaw, 1);
//    writeAngleError(estimatedRoll, groundRoll, 2);
    //writeAccuracyvsNoTree(NUMBER_TREES,noseA,poseA);
//    vector<double> error = Error(trueMean,testGt);
//    cout << "nose error : " << error[0] << " mm " << endl;
//    cout << "angle error : " << error[1] << " degrees " << endl;

//    string landmarkfilepath = "/Users/Knowhow10/Desktop/FaceWarehouse_Data_0.part1/Tester_1/TrainingPose";
//    vector<vector<float>> landmarks;
//    string fname = "landmarks.txt";
//    vector<int> lms = loadLandmarks(fname);
//    landmarks = readLandMarks(landmarkfilepath+"/pose_1.land");
//    vector<vector<vector<float>>> trainLandmarks;
//    vector<Mat> trainSet2d;
//    vector<boundingBox> bboxs;
//    boundingBox bigBbox;
//    bigBbox.height = 400;
//    bigBbox.width = 560;
//    bigBbox.x = 0;
//    bigBbox.y = 0;
//    face_cascade.load(face_cascade_name);
//    load2dTrainingset(trainSet2d,trainLandmarks,bboxs);
//    //cout << trainLandmarks[0][0][0] << " " << trainLandmarks[0][0][1]<< endl;
//    cout << trainSet2d.size() << endl;
//    vector<HPatch> wholeDataSet;
//    for( int i = 0; i < trainSet2d.size(); i++ ){
//        cout << "#####     i = " << i << endl;
//        PatchSet pS(PATCH_SET_SIZE);
//        
//        //cout << "bbox big" << bbox.x <<  " " << bbox.y << endl;
//        //depthTo3D = get3dFromDepth(trainingSet[i],calM);
//        pS.getRandomPatches2d(bboxs[i], trainSet2d[i], trainLandmarks[i],0);
//        //cout << "number of positive patches" << pS.pSet.size() << endl;
//        pS.getRandomPatches2d(bigBbox, trainSet2d[i], trainLandmarks[i],1);
//        //cout << " gt : " << groundTruthSetBin[i][0] << " " << groundTruthSetBin[i][1] << " " << groundTruthSetBin[i][2] << " " << groundTruthSetBin[i][3] << " " << groundTruthSetBin[i][4] << " " << groundTruthSetBin[i][5] << endl;
//        cout << "number of patches" << pS.pSet.size() << endl;;
//        for (int j = 0; j <  PATCH_SET_SIZE; j++){
//            pS.pSet[j].index = i;
//            wholeDataSet.push_back(pS.pSet[j]);
//        }
//        cout << "number of patches " << wholeDataSet.size() << endl;
//    }
//    //cout << wholeDataSet[123].index <<endl;
//    //imshow("Display window", trainSet2d[156] );
//    //cout << gray_image << endl;
//    //waitKey(0);
//    Mat *trainSet2dArray = new Mat[trainSet2d.size()];
//    for(int tr = 0; tr < trainSet2d.size(); tr++){
//        trainSet2dArray[tr] = trainSet2d[tr];
//    }
//    trainSet2d.clear();
//    DForest forest(NUMBER_TREES);
//    forest.growForest(wholeDataSet, trainSet2dArray);
//    forest.writeForest();
//    vector<vector<float>> accuracy1,accuracy2,accuracy3;
//    vector<float> timeT;
//    //for(int noft = 5; noft <= 10; noft++){
//        DForest forest(NUMBER_TREES);
//        forest.loadTree();
//        
//        
//        Mat testImage,testImage2;
//        int count[17];
//        int count2[17];
//        int count3[17];
//        for(int c = 0; c< 17; c++){
//            count[c] = 0;
//            count2[c] = 0;
//            count3[c] = 0;
//        }
//        
//        int countI = 0;
//        vector<float> times;
//        //for(int p = 1; p <= 10; p++){
//           // for(int p2 = 0; p2 <20 ; p2++){
//                //string testimagepath = "/Users/Knowhow10/Desktop/projectFacedata/Tester_"+convertInt(p)+"/TrainingPose/pose_"+convertInt(p2);
//                string testimagepath = "/Users/Knowhow10/Desktop/projectFacedata/Tester_105/TrainingPose/pose_15";
//                testImage = imread(testimagepath+".png",CV_LOAD_IMAGE_GRAYSCALE);
//                testImage2 = imread(testimagepath+".png",CV_LOAD_IMAGE_COLOR);
//                Mat testSum = Mat(481,641,CV_64FC3);
//                integral(testImage, testSum);
//                vector<vector<float>> templms1,templms2;
//                vector<Point> dots;
//                templms1 = readLandMarks(testimagepath+".land");
//                //cout << templms1.size() <<endl;
//                pickWantedLandmarks(templms1, templms2);
//                dots = getdots(templms1);
//                face_cascade.load( face_cascade_name );
//                vector<Rect> face;
//                detectAndDisplay( testImage2,face );
//                //if(face.size() == 1){
//        //            cout << "face x :" << face.x- face.width<<" face y: " << face.y <<endl;
//        //            cout << "face w :" << face.width << " face h: " << face.height <<endl;
//                    boundingBox testBbox ;
//        //            cout << testBbox.x << " " << testBbox.y << endl;
//        //            cout << testBbox.width << " " << testBbox.height << endl;
//                    countI++;
//                    testBbox.x = face[0].x;
//                    testBbox.y = face[0].y;
//                    testBbox.width = face[0].width;
//                    testBbox.height =face[0].height;
//                    clock_t t = clock();
//                    forest.regressionEstimation2d(testSum, testBbox, templms2);
//                    //cout << forest.estimatedMean.size() << endl;
//                    vector<vector<float>> mean ;//= computeMeanVector2d(forest.estimatedMean);
//                    cout << forest.estimatedMean.size() <<endl;
//                    removeoutliars(mean,forest.estimatedMean);
//                    t = clock() - t;
//                    times.push_back((float)t/CLOCKS_PER_SEC);
//                    cout << "time " << (float)t/CLOCKS_PER_SEC << endl;
//            //        for(int j = 0; j < forest.estimatedMean.size() ; j++){
//            //            for(int i =0; i < 17; i= i+2){
//            //                cout << "feature " << i << " x: " << forest.estimatedMean[j][i] << " y: " << forest.estimatedMean[j][i+1] << endl;
//            //            }
//            //        }
//                    for(int j = 0; j < mean.size(); j++){
//                        cout << "feature " << j << endl;
//                        cout << "estimated x : " << mean[j][0] << " y : " << mean[j][1] << endl;
//                        cout << "true x : " << templms2[j][0] << " true y: " << templms2[j][1] <<endl;
//                        float err = sqrt((mean[j][0] - templms2[j][0])*(mean[j][0] - templms2[j][0]) + ((mean[j][1]-templms2[j][1])*(mean[j][1]-templms2[j][1])));
//                        cout << " error " << err <<endl;
//                        if(err <= 0.018)
//                            count[j] ++;
//                        if(err <= 0.015)
//                            count2[j] ++;
//                        if(err <= 0.013)
//                            count3[j] ++;
//                    }
                    //forest.estimatedMean.clear();
                //}
            //}
            
        //}
//        vector<float> temp1;
//        for(int e = 0; e < 17; e++){
//            //cout <<(float)count[e]/countI << endl;
//            temp1.push_back((float)count[e]/countI);
//            
//            
//        }
//        accuracy1.push_back(temp1);
//        
//        vector<float> temp2;
//        for(int e = 0; e < 17; e++){
//            
//            
//            //cout <<(float)count[e]/countI << endl;
//            temp2.push_back((float)count2[e]/countI);
//            
//            
//            
//        }
//        accuracy2.push_back(temp2);
//        
//        vector<float> temp3;
//        for(int e = 0; e < 17; e++){
//            
//            
//            cout <<(float)count[e]/countI << endl;
//            temp3.push_back((float)count3[e]/countI);
//            
//            
//            
//        }
//        accuracy3.push_back(temp3);
//        float atime = 0;
//        for(int i =0; i < times.size(); i++){
//            atime = atime + times[i];
//        }
//        cout << atime / times.size() <<endl;
//        timeT.push_back(atime / times.size());
//    }
//    cout << accuracy1.size() <<endl;
//    ofstream fOut,fOut2,fOut3;
//    fOut.open("2daccuracy0.018.txt");
//    fOut2.open("2daccuracy0.015.txt");
//    fOut3.open("2daccuracy0.013.txt");
//    for(int j =0; j < accuracy1.size(); j++){
//        fOut << j+1 << " " << timeT[j] << " ";
//        fOut2 << j+1 << " " << timeT[j] << " ";
//        fOut3 << j+1 << " " << timeT[j] << " ";
//        for(int i =0; i <17; i++){
//            //cout << accuracy1[j][i] << " " << j+1 << " " << timeT[j] << endl;
//            fOut << accuracy1[j][i] << " " ;
//            fOut2 << accuracy2[j][i] << " " ;
//            fOut3 << accuracy3[j][i] << " " ;
//        }
//        fOut2 << endl;
//        fOut3 << endl;
//        fOut << endl;
//    }
//    
//    fOut.close();
//    fOut2.close();
//    fOut3.close();
//    vector<Point> estimatedDots = getdots(mean);
//    vector<Point> trueDots = getdots(templms2);
//    //displydots(trueDots, testImage2,1);
//    //displydots(estimatedDots, testImage2,0);
//    
//    imshow("Display window", testImage2 );
//    waitKey(0);
    DForest forest(NUMBER_TREES);
    forest.loadTree();
    VideoCapture stream1(0);   //0 is the id of video device.0 if you have only one camera.
    
    if (!stream1.isOpened()) { //check if video device has been initialised
        cout << "cannot open camera";
    }
    
    //unconditional loop
    while (true) {
        Mat cameraFrame;
        stream1.read(cameraFrame);
        face_cascade.load( face_cascade_name );
        vector<Rect> face;
        boundingBox testBbox;
        vector<vector<float>> templms2;
        //vector<vector<float>>
        //cout << cameraFrame.cols << endl;
        detectAndDisplay( cameraFrame,face );
        testBbox.x = face[0].x;
        testBbox.y = face[0].y;
        testBbox.width = face[0].width;
        testBbox.height =face[0].height;
        Mat testSum = Mat(481,641,CV_64FC3);
        integral(cameraFrame, testSum);
        vector<vector<float>> mean ;
        forest.regressionEstimation2d(testSum, testBbox, templms2);
        removeoutliars(mean,forest.estimatedMean);
        vector<Point> estimatedDots = getdots(mean);
        displydots(estimatedDots, cameraFrame,0);
        imshow("cam", cameraFrame);
        if (waitKey(30) >= 0)
            break;
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
void loadTrainingSet(string depth_path,vector<int16_t*>& trainingSet){
    
    
    int16_t* img = new int16_t[640*480];
    string subS = "";
    string depthFname = depth_path;
    ifstream fInp;
    for(int i = 4; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
		depthFname = depth_path  + "frame_00"+subS+"_depth.bin";
		fInp.open(depthFname);
        if(fInp.is_open()){
            fInp.close();
            img = loadDepthImageCompressed( depthFname.c_str());
            trainingSet.push_back(img);
        }
		
	}
    
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

void loadGroundSetBin(string bin_path,vector<vector<float>>& groundSetBin){
    
    vector<float> tempG;
    string subS = "";
    string poseFname = bin_path;
    ifstream fInp;
    for(int i = 4; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
        poseFname = bin_path  + "frame_00"+subS+"_pose.bin";
        fInp.open(poseFname);
        if(fInp.is_open()){
            fInp.close();
            tempG = loadPoseBin(poseFname);
            groundSetBin.push_back(tempG);
        }
    }
    
    
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
    //cout << gt[0] << " " << gt[1] << " " << gt[2] << endl;
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

vector<double> Error(vector<float> estimated, vector<float> gt){
    vector<double> err;
    double noseError = 0;
    double angleError = 0;
    for(int i = 0; i < 3; i++){
        noseError = noseError + (gt[i]-estimated[i])*(gt[i]-estimated[i]);
        angleError = angleError + (gt[i+3]-estimated[i+3])*(gt[i+3]-estimated[i+3]);
    }
    err.push_back(sqrt(noseError));
    err.push_back(sqrt(angleError));
    return err;
}

void loadRGBImages(vector<Mat>& rgb, string filepath){
    Mat temp;
    string subS = "";
    string filename = "";
    for(int i = 4; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
        filename = filepath  + "frame_00"+subS+"_rgb.png";
        temp = imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
        rgb.push_back(temp);
    }
}

void loadMaskImages(vector<Mat>& rgb, string filepath){
    Mat temp;
    string subS = "";
    string filename = "";
    ifstream fInp;
    for(int i = 4; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
        filename = filepath  + "frame_00"+subS+"_depth_mask.png";
        fInp.open(filename);
        if(fInp.is_open()){
            fInp.close();
            temp = imread(filename, CV_LOAD_IMAGE_COLOR);
            rgb.push_back(temp);
        }
    }
}

void saveBBox(vector<boundingBox> masks,string fname){
    ofstream fOut;
    fOut.open(fname);
    for(int i = 0; i < masks.size(); i++){
        fOut << masks[i].x << " " << masks[i].y << " " << masks[i].width << " " << masks[i].height << endl;
    }
    fOut.close();
}

void loadBBox(vector<boundingBox>& bboxSet, int size){
    boundingBox temp;
    for(int j = 0; j < size; j++){
        if(j!=5){
            string fname = "headmask"+convertInt(j+1)+".txt";
            ifstream fInp;
            fInp.open(fname);
            for(int i = 0; i < DATA_SET_SIZE; i ++){
                fInp >> temp.x >> temp.y >> temp.width >> temp.height;
                bboxSet.push_back(temp);
                //cout << bboxSet[i].x << " " << bboxSet[i].y << " " << bboxSet[i].width << " " << bboxSet[i].height << endl;
            }
            fInp.close();
        }
    }
    
}

boundingBox findMask(Mat mask){
    boundingBox bbox;
    int min_x = 640;
	int min_y = 480;
	int max_x = 0;
	int max_y = 0;
    
    
	for(int y = 0; y < 480; y++)
	{
	    
	    for(int x = 0; x < 640; x++){
            
	    	if( mask.at<int>(y,x) > 0) {
                
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

void findMasks(vector<Mat> masks, vector<boundingBox>& headMasks){
    boundingBox temp;
    for(int i = 0; i < masks.size(); i++){
        temp = findMask(masks[i]);
        headMasks.push_back(temp);
    }
}

void eliminateOutliars(vector<vector<float>> estimated,vector<float>& newMean){
    vector<float> tempMean;
    //vector<float> err;
    vector<vector<float>> afterElimination;
    tempMean = computeMeanVector(estimated);
    for(int i =0; i < estimated.size(); i++){
        vector<double> err = Error(estimated[i], tempMean);
        //cout << err[0] << " mm " << endl;
        if(err[0] <= ERROR_TH)
            afterElimination.push_back(estimated[i]);
            
    }
    cout << "after size " << afterElimination.size() << endl;
    if(afterElimination.size() >0)
        newMean = computeMeanVector(afterElimination);

}

void runTestImage(int i,DForest forest,vector<float> calM,bool& noseA,bool& poseA){
    string subS = "";
    string testImagefile = "";
    string testImageGTfile = "";
    if(i<10)
        subS = "00"+convertInt(i);
    else if (i>=10 && i < 100)
        subS = "0"+convertInt(i);
    else
        subS = convertInt(i);
    testImagefile = "/Users/Knowhow10/Downloads/kinect_head_pose_db/09/frame_00"+subS+ "_depth.bin";
    testImageGTfile = "/Users/Knowhow10/Desktop/db_annotations/09/frame_00"+subS+"_pose.bin";

    int16_t* testImage = loadDepthImageCompressed(testImagefile.c_str());
    boundingBox testBbox = getBoundingBox(testImage);
    //vector<vector<float>> groundTruthSetBin = loadGroundSetBin(bin_path);
    Mat* channels = new Mat[3];
    Mat depthI = Mat(480, 640, CV_32FC3);
    Mat sumI = Mat(481, 641, CV_64FC3);
    Mat img3D;
    getIntegralImageDepthChannel(depthI,testImage,channels,calM,img3D);
    integral(channels[2],sumI);
    vector<float> testGt = loadPoseBin(testImageGTfile);
    vector<float> trueMean;
    clock_t t = clock();
    forest.regressionEstimation(sumI, testBbox, testGt, img3D);
    
    vector<float> after;
    eliminateOutliars(forest.estimatedMean,after);
    t = clock() -t;
    cout << "time " << (float)t/CLOCKS_PER_SEC <<endl;
    if(after.size() > 0){
        cout << "estimated " << after[0] << " " << after[1] << " " << after[2] << " " << after[3] << " " << after[4] << " " << after[5] << endl;
        cout << "ground_T " << testGt[0] << " " << testGt[1] << " " << testGt[2] << " " << testGt[3] << " " << testGt[4] << " " << testGt[5] << endl;
        estimatedPitch.push_back(after[3]); groundPitch.push_back(testGt[3]);
        estimatedYaw.push_back(after[4]);   groundYaw.push_back(testGt[4]);
        estimatedRoll.push_back(after[5]);  groundRoll.push_back(testGt[5]);
        vector<double> error = Error(after,testGt);
        cout << "nose error : " << error[0] << " mm " << endl;
        cout << "angle error : " << error[1] << " degrees " << endl;
        if(error[0] <= 20)
            noseA = 1;
        if(error[1] <= 15)
            poseA = 1;
    }
    forest.estimatedMean.clear();
    
}

vector<vector<float>> readLandMarks(string fname){
    vector<vector<float>> landmarks;
    vector<float> temp;
    float tempf;
    ifstream fInp;
    fInp.open(fname);
    if(fInp.is_open()){
        fInp>> tempf;
        for(int i = 0; i < 74; i++){
            fInp >> tempf;
            //cout << tempf;
            temp.push_back(tempf);
            fInp >> tempf;
            temp.push_back(1-tempf);
            landmarks.push_back(temp);
            temp.clear();
        }
    }
    fInp.close();

    return landmarks;
}

vector<Point> getdots(vector<vector<float>> landmarks){
    vector<Point> dots;
    for(int i = 0; i < NUM_LANDMARKS; i++){
        //cout << "x : " << 640*landmarks[i][0] << " y " << 480 - 480*landmarks[i][1] <<endl;
        Point temp = Point(640*landmarks[i][0],480*landmarks[i][1]);
        dots.push_back(temp);
    }
    return dots;
}
void displydots(vector<Point> dots,Mat image, bool t){
    int r = 0;
    int g = 0;
    if(t == 0)
        g = 255;
    else
        r = 255;
    for(int i = 0; i < NUM_LANDMARKS; i=i+1){
        if(t == 1)
        circle(image, dots[i],1,CV_RGB(255,0,0),8);
        string s = convertInt(i);
        if(t == 0)
        cv::putText(image,s, dots[i], CV_FONT_HERSHEY_PLAIN, 0.4,CV_RGB(r,g,0));
    }
}

vector<int> loadLandmarks(string fname){
    vector<int> lms;
    int temp = 0;
    ifstream fInp;
    fInp.open(fname);
    while(!fInp.eof()){
        fInp >> temp;
        lms.push_back(temp);
    }
    return lms;
}

boundingBox findBbox2d(vector<Point> dots){
    boundingBox bbox;
    bbox.x = dots[1].x-10;
    //cout << " dot7 y " << dots[7].y << endl;
    bbox.y = dots[7].y-200;
    //cout << dots[13].x - dots[1].x << endl;
    bbox.width = dots[13].x - dots[1].x + 20;
    bbox.height = dots[7].y - dots[7].y+200 + 10;
    return bbox;
}

void load2dTrainingset(vector<Mat>& trainingset2d, vector<vector<vector<float>>>& landmarks, vector<boundingBox>& bboxs){
    Mat temp,temp2;
    vector<vector<float>> templms1,templms2;
    vector<Point> tempdots;
    boundingBox tempbbox;
    string folderpath = "/Users/Knowhow10/Desktop/projectFacedata/Tester_";
    string pngfilename = "pose_";
    string landmarkfilename = "";
    int count = 0;
    for(int i = 1; i <= NO_TEST_SUBJECTS; i++){
        for(int j = 0; j < NUMBER_POSE; j++){
            Mat Sumimage = Mat(481,641,CV_64FC3);
            string subs1 = convertInt(i) + "/TrainingPose/";
            temp = imread(folderpath+subs1+pngfilename+convertInt(j)+".png",CV_LOAD_IMAGE_GRAYSCALE);
            temp2 = imread(folderpath+subs1+pngfilename+convertInt(j)+".png",CV_LOAD_IMAGE_COLOR);
            //cout << folderpath+subs1+pngfilename+convertInt(j)+".png" << endl;
//            CascadeClassifier tempFace_cascade;
//            tempFace_cascade.load(face_cascade_name);
            vector<Rect> faces;
            detectAndDisplay( temp2,faces );
            if(faces.size() == 1){
                count++;
                //cout << faces[0].x << " " << faces[0].width <<endl;
                tempbbox.x = faces[0].x - 0.1*faces[0].width;
                tempbbox.y = faces[0].y - 0.1*faces[0].height;
                tempbbox.width = faces[0].width*1.2;
                tempbbox.height = faces[0].height*1.2;
                //cout << tempbbox.x << " " << tempbbox.width <<endl;
                cout << " image " << count << " read ok " << endl;
                integral(temp,Sumimage);
                templms1 = readLandMarks(folderpath+subs1+pngfilename+convertInt(j)+".land");
                pickWantedLandmarks(templms1,templms2);
                tempdots = getdots(templms1);
                //tempbbox = findBbox2d(tempdots);
                trainingset2d.push_back(Sumimage);
                bboxs.push_back(tempbbox);
                landmarks.push_back(templms2);
                templms1.clear();
                tempdots.clear();
                templms2.clear();
            }
           // cout << i << " " << j << " done " << endl;
        }
        
    }
    cout << trainingset2d.size() << endl;
    
}

void pickWantedLandmarks(vector<vector<float>> landmarks, vector<vector<float>>& newlandmarks){
    
    string fname = "landmarks.txt";
    vector<int> lms = loadLandmarks(fname);
    for(int i = 0; i < lms.size(); i++){
        newlandmarks.push_back(landmarks[lms[i]]);
    }
}

vector<vector<float>> computeMeanVector2d(vector<vector<float>> data){
    vector<vector<float>> bigmean;
    for(int i = 0; i < 34; i= i+2){
        vector<float> smallmean;
        float tempx = 0;
        float tempy = 0;
        for(int j = 0; j < data.size(); j++){
            tempx = tempx + data[j][i];
            tempy = tempy + data[j][i+1];
        }
        tempx = tempx / data.size();
        tempy = tempy / data.size();
        smallmean.push_back(tempx); smallmean.push_back(tempy);
        bigmean.push_back(smallmean);
        
    }
    //cout << bigmean.size() <<endl;
    return bigmean;
}

void removeoutliars(vector<vector<float>>& newMean,vector<vector<float>> data){
    vector<vector<float>> oldmean = computeMeanVector2d(data);
    //cout << "data size " << data[18].size() <<endl;
    for(int i = 0; i < 17; i++){
        //cout << " i " << i <<endl;
        float tempx = 0;
        float tempy = 0;
        vector<float> smallmean;
        int countx = 0;
        int county = 0;
        for(int j = 0; j < data.size(); j++){
//            cout << "j " << j << endl;
//            cout << " 2i " << i + i <<endl;
//            cout << " 2i+1 " << i+i+1 <<endl;
            float errx = oldmean[i][0] - data[j][i+i];
            float erry = oldmean[i][1] - data[j][i+i+1];
//            cout << "errx " << errx << endl;
//            cout << "erry " << erry << endl;
            if(errx <= ERROR_TH_2D){
                tempx = tempx + data[j][i+i];
                
                countx++;
            }
            if(erry <= ERROR_TH_2D){
                tempy = tempy + data[j][i+i+1];
                county++;
            }
        }
        tempx = tempx /countx;
        tempy = tempy /county;
        smallmean.push_back(tempx); smallmean.push_back(tempy);
        newMean.push_back(smallmean);
    }
    //cout << " new mean size " << newMean.size() <<endl;
    
}

void writeAccuracyvsNoTree(int trees,vector<float> noseAccuracy, vector<float> poseAccuracy){
    ofstream fOut;
    fOut.open("NotreesVSaccuracy.txt");
    for(int i = 0; i < trees; i++){
        fOut << noseAccuracy[i] << " " << poseAccuracy[i] << " " << i+1 <<endl;
    }
    fOut.close();
}

void writeAngleError(vector<float> estimated, vector<float> ground, int angle){
    ofstream fOut;
    string fname = "";
    if(angle == 0)
        fname= "pitchError.txt";
    if(angle == 1)
        fname= "yawError.txt";
    if(angle == 2)
        fname= "rollError.txt";
    fOut.open(fname);
    for(int i = 0; i < estimated.size(); i++){
        fOut << estimated[i] << " " << ground[i] << " " << i << endl;
        
    }
    fOut.close();
}

void detectAndDisplay( Mat frame, vector<Rect>& faces )
{
    //std::vector<Rect> faces;
    Mat frame_gray;
    
    cvtColor( frame, frame_gray, CV_BGR2GRAY );
    equalizeHist( frame_gray, frame_gray );
    
    //-- Detect faces
    face_cascade.detectMultiScale( frame_gray, faces, 1.1, 2, 0|CV_HAAR_SCALE_IMAGE, Size(30, 30) );
    //cout << faces.size() << endl;
    rectangle(frame, faces[0], Scalar( 255, 0, 255 ));
    //face = faces;
//    for( size_t i = 0; i < faces.size(); i++ )
//    {
//        Point center( faces[i].x + faces[i].width*0.5, faces[i].y + faces[i].height*0.5 );
//        ellipse( frame, center, Size( faces[i].width*0.5, faces[i].height*0.5), 0, 0, 360, Scalar( 255, 0, 255 ), 4, 8, 0 );
//        
//        Mat faceROI = frame_gray( faces[i] );
//
//    }
    
    //-- Show what you got
    //imshow( window_name, frame );
}