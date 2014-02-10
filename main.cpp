#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <sstream>
#include <cmath>
using namespace std;
struct ground_Truth{
	float rotationMatrix[9];
	float nosePosition[3];
};


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

	if(success)
		return depth_img;
	else{
		delete [] depth_img;
		return NULL;
	}
}

float* matrixToEuler(float* m){
	float euler[3] = {};
	euler[0] = atan2 (m[7],m[8])*180/3.1415926;
	euler[1] = atan2 (-1*m[6],sqrt(m[7]*m[7]+m[8]*m[8]) )*180/3.1415926;
	euler[2] = atan2 (m[3],m[0])*180/3.1415926;
	return euler;
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





int main(){

	string depth_path = "D:/kinect_head_pose_db/01/";
	string gt_path = "D:/kinect_head_pose_db/01/";
	string cal_filename = depth_path + "depth.cal";
	
	ifstream fInp;
	fInp.open(cal_filename.c_str());
	if (!fInp.is_open()){
		cerr << "depth.cal file not found in the same folder as the depth image! " << endl;
		return -1;
	}
	//read intrinsics only
	float depth_intrinsic[9];	
	for(int i =0; i<9; ++i){

		fInp >> depth_intrinsic[i];
		cout << depth_intrinsic[i] << endl;
	}
	fInp.close();

	vector<int16_t*> img2;
	vector<ground_Truth> gt;
	int16_t* img = new int16_t[640*480];
	ground_Truth gt1;
	string filename = depth_path ;
	string pose_filename = gt_path ;
	string subS = "";
	for(int i = 3; i < 500; i++){
		if(i<10)
			subS = "00"+convertInt(i);
		else if (i>=10 && i < 100)
			subS = "0"+convertInt(i);
		else
			subS = convertInt(i);
		filename = depth_path  + "frame_00"+subS+"_depth.bin";
		pose_filename = gt_path  + "frame_00"+subS+"_pose.txt";
		img = loadDepthImageCompressed( filename.c_str());
		loadPoseTxt(pose_filename,gt1);
		img2.push_back(img);
		gt.push_back(gt1);
	}
	int p_x = 0;
	int p_y = 0;
	for (int i = 0; i < 20; i++){
		//p_x = rand() % 400;
		//p_y = 200 + rand() % 180 ;
		//cout << "xp : "<< p_x << " yp : "<< p_y << endl;
		//cout << "depth : " << img2.size() << endl;
		//cout << "depth : " << gt.size() << endl;
		//cout << "rotation matrix " << endl << gt[i].rotationMatrix[0] << " " << gt[i].rotationMatrix[1] << " "<< gt[i].rotationMatrix[2] << " " << endl;
		//cout << gt[i].rotationMatrix[3] << " " << gt[i].rotationMatrix[4] << " "<< gt[i].rotationMatrix[5] << " " << endl;
		//cout << gt[i].rotationMatrix[6] << " " << gt[i].rotationMatrix[7] << " "<< gt[i].rotationMatrix[8] << " " << endl;
		//cout << "pose : " << gt[i].nosePosition[0] << " " << gt[i].nosePosition[1] << " "<< gt[i].nosePosition[2] << " " << endl;
	}
	
	//cout << "rotation matrix " << endl << gt[1].rotationMatrix[0] << " " << gt[1].rotationMatrix[1] << " "<< gt[1].rotationMatrix[2] << " " << endl;
	
	//cout << gt[1].rotationMatrix[3] << " " << gt[1].rotationMatrix[4] << " "<< gt[1].rotationMatrix[5] << " " << endl;
	
	//cout << gt[1].rotationMatrix[6] << " " << gt[1].rotationMatrix[7] << " "<< gt[1].rotationMatrix[8] << " " << endl;
	float *m = matrixToEuler(gt[1].rotationMatrix);
	cout << m[0] << " " << m[1] << " " << m[2] << endl;
	bool have_gt = false;
	float *gt2;
	//try to read in the ground truth from a binary file
	string pose_filename1 = "D:/db_annotations/01/frame_00004_pose.bin";
	FILE* pFile = fopen(pose_filename1.c_str(), "rb");
	if(pFile){

		have_gt = true;
		have_gt &= ( fread( &gt2[0], sizeof(float),6, pFile) == 6 );
		fclose(pFile);
	}
	cout << gt2[3] << " " << gt2[4] << " " << gt2[5] << endl;
	return 0;

}

