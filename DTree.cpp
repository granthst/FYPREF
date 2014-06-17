//
//  DTree.cpp
//  fypFirstDraft
//
//  Created by callen on 11/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "DTree.h"
#include <ctime>
#define n_F 100
#define n_T 10
#define THESHOLD_MAX 200
#define DEPTH_TREE 20
#define TRACE_MAX 300
#define SAMPLE_PATCHES 1000
#define P_TH 1
#define T_P 0.7
#define TD_LAMBDA 0.125
#define step_x 5
#define step_y 5



void Node::generateRandomSubPatches(const vector<HPatch>& PS){
    sub_patch temp;
    vector<sub_patch> tempSub;
    for (int j = 0; j < n_F; j++){
        for(int i = 0; i < 2; i++){
            temp.x = rand() % 80;
            temp.y = rand() % 80;
            if(temp.x <= 40)
                temp.w = 1 + rand()% 40;
            else
                temp.w = 1 + rand() % (80 - temp.x);
            if(temp.y <= 40)
                temp.h = 1 + rand() % 40;
            else
                temp.h = 1 + rand() % (80 - temp.y);
            tempSub.push_back(temp);
           // cout << tempSub[i].x << endl;
        }
        rectangles.push_back(tempSub);
        //cout << " rec " << rectangles[j][0].x << endl;
        tempSub.clear();
        //cout << " rec " << rectangles[j][0].x << endl;
    }
    
}

//generate a set of random Thresholds for binary test
void Node::generateRandomThreshold(const vector<HPatch>& PS){
    
    int rt = 0;
    for (int i = 0; i < n_T; i++){
        rt = -THESHOLD_MAX + rand() % (2*THESHOLD_MAX);
        rT.push_back(rt);
    }
    
}

void Node::setPatchSetBeforeSplit( vector<HPatch> PS){

    beforeSplit = PS;
}

void Node::findBestT(vector<HPatch> PS, Mat* itegralImage){
    
    //cout << PS.size() << endl;
    float infoGainTemp = 0;//numeric_limits<int>::min();
    //cout << "int min " << infoGainTemp << endl;
    vector<HPatch> tempLeft, tempRight;
    int x1,y1,w1,h1,x2,y2,w2,h2;
    float tempSB = 0;
    float tempRT = 0;
    sub_patch tempS1,tempS2;
    generateRandomThreshold(PS);
    generateRandomSubPatches(PS);
    //sub_patch temp;
    //ifstream fInp;
    //ofstream fOut;
    //fOut.open(treeFilename,fstream::app);
    //fOut.open(integralImageFilename);
    //fInp.open(integralImageFilename);
    /*for( int i = 0; i < PS.size(); i++){
        PS[i].chooseSubPatches(rectangles[i]);
        PS[i].setSubPatchDistance(dImage);
    }*/
    for(int k = 0; k < rT.size(); k++ ){
        
        tempRT = rT[k];
        //cout << "threshold " << tempRT << endl;
        for(int j = 0; j < rectangles.size(); j++){
            for( int i = 0; i < PS.size(); i++){
                //cout << itegralImage[PS[i].index].at<int>(rectangles[j][0].x+rectangles[j][0].w,rectangles[j][0].y+rectangles[j][0].h) << endl;
                PS[i].chooseSubPatches(rectangles[j],itegralImage[PS[i].index]);
                x1 = rectangles[j][0].x; x2 = rectangles[j][1].x;
                y1 = rectangles[j][0].y; y2 = rectangles[j][1].y;
                w1 = rectangles[j][0].w; w2 = rectangles[j][1].w;
                h1 = rectangles[j][0].h; h2 = rectangles[j][1].h;
                tempSB = PS[i].subPDistance;
                //cout << "sb difference " << tempSB << endl;
                if(tempSB > tempRT)
                    tempLeft.push_back(PS[i]);
                else
                    tempRight.push_back(PS[i]);

                //PS[i].rectangles.clear();
            }
            //rT.push_back(tempRT);
            //leftSplit = tempLeft;
            //rightSplit = tempRight;
           //cout << j << " left child patch size : " << tempLeft.size() << " right child patch size : " << tempRight.size() << endl;
            //cout << "info gain " << infoGain2d(beforeSplit, tempLeft, tempRight) << endl;
            
           //cout << "rt " << tempRT << endl;
            //cout << " tempsb " << tempSB << endl;

                if (infoGainTemp < infoGain2d(beforeSplit, tempLeft, tempRight)){
                    leftSplit = tempLeft;
                    rightSplit = tempRight;
                    infoGainTemp = infoGain2d(beforeSplit, tempLeft, tempRight);
                    tempS1.x = x1;   tempS1.y = y1;
                    tempS1.w = w1;   tempS1.h = h1;
                    tempS2.x = x2;   tempS2.y = y2;
                    tempS2.w = w2;   tempS2.h = h2;
                    best_T = tempRT;
                }
            
            
            tempLeft.clear();
            tempRight.clear();

        }

    if(((leftSplit.size()==0) && (rightSplit.size()==0)) || depth == DEPTH_TREE){
        //cout << " leaf node " << endl;
        computePositiveP();
        isleaf = 1;
    }
    //cout << "depth : " << depth << endl;
    //cout << " isleaf " << isleaf << endl;
        //cout << best_T << endl;
        //cout << "ig " << infoGainTemp << endl;
    }

   bestF.push_back(tempS1);  bestF.push_back(tempS2);
    

    //cout << "at depth " << depth << " beforesplit patch size " << beforeSplit.size() << endl;
    //cout << "at depth " << depth << " left child patch size : " << leftSplit.size() << " right child patch size : " << rightSplit.size() << endl;


}

vector<Node> Node::makeTreeNoRecursion(vector<HPatch> PS, Mat* integralImage) {
    
    int iterationNo = 0;
    vector<Node> pointers;
    Node initialPointer = *this;
    initialPointer.setPatchSetBeforeSplit(PS);
    initialPointer.findBestT(PS,integralImage);
    pointers.push_back(initialPointer);
    int prevNoPointers = 0;
    while (iterationNo < DEPTH_TREE) {
        int noPointers = (int) pointers.size();
        
        cout << "iterationNo = " << iterationNo << ", noPointers = " << noPointers << endl;
        
        for (int j = prevNoPointers; j < noPointers; j++) {
            Node parentPointer = pointers[j];
            if(parentPointer.isleaf == 0){
                Node leftPointer = Node(iterationNo);
                leftPointer.setPatchSetBeforeSplit(parentPointer.leftSplit);/***/
                leftPointer.findBestT(leftPointer.beforeSplit, integralImage);
                leftPointer.depth = iterationNo+1;
                pointers.push_back(leftPointer);
                //pointers[j].left_child = leftPointer;
                
                Node rightPointer = Node(iterationNo);
                rightPointer.setPatchSetBeforeSplit(parentPointer.rightSplit);/***/
                rightPointer.depth = iterationNo+1;
                rightPointer.findBestT(rightPointer.beforeSplit, integralImage);
                pointers.push_back(rightPointer);
                //pointers[j]->right_child = rightPointer;
            }
        }
        prevNoPointers = noPointers;
        iterationNo++;
    }
    return pointers;
}

void Node::splitPatchSet(int bin_test){
    
    vector<HPatch> lS,rL;
    
    
    
    leftSplit = lS;
    rightSplit = rL;
}

vector<float> Node::computeMeanVector(vector<HPatch> sPvector){
    vector<float> mean;
    for(int i = 0; i < 6; i++){
        float temp = 0.0;
        for(int j = 0; j < sPvector.size(); j++){
            temp = sPvector[j].groundT[i] + temp;
        }
        temp = temp/sPvector.size();
        mean.push_back(temp);
    }
    return mean;
}

vector<float> Node::computeMeanVector2d(vector<HPatch> sPvector, int lm){
    vector<float> mean;
    for(int i = 0; i < 2; i++){
        float temp = 0.0;
        for(int j = 0; j < sPvector.size(); j++){
            temp = sPvector[j].groundT2d[lm][i] + temp;
        }
        temp = temp/sPvector.size();
        mean.push_back(temp);
    }
    return mean;
}

vector<vector<float>> Node::computeCovariance(vector<HPatch> sPvector, bool angleOrNose){
    vector<vector<float>> Covariance;
    vector<float> meanVector = computeMeanVector(sPvector);
    int offset = 0;
    if(angleOrNose == 0)
        offset = 3;
    for(int i = 0; i < 3; i++){
        vector<float> col;
        for(int j = 0; j < 3; j++){
            float temp = 0.0;
            for(int k = 0; k < sPvector.size(); k++){
                temp = temp + (sPvector[k].groundT[offset+i] - meanVector[offset+i]) *(sPvector[k].groundT[offset+j] - meanVector[offset+j]);
            }
            float cov;
            if(sPvector.size() != 0)
                cov = temp / sPvector.size();
            else
                cov = 0;
            col.push_back(cov);
        }
        Covariance.push_back(col);
    }
    //cout << "size " << sPvector.size() << endl;
    return  Covariance;

    //cv = Covariance;
}

vector<vector<float>> Node::computeCovariance2d(vector<HPatch> sPvector, int lm){
    vector<vector<float>> Covariance;
    vector<float> meanVector = computeMeanVector2d(sPvector,lm);
    for(int i = 0; i < 2; i++){
        vector<float> col;
        for(int j = 0; j < 2; j++){
            float temp = 0.0;
            for(int k = 0; k < sPvector.size(); k++){
               temp = temp + (sPvector[k].groundT2d[lm][i] - meanVector[i]) *(sPvector[k].groundT2d[lm][j] - meanVector[j]);
            }
            float cov;
            if(sPvector.size() != 0)
                cov = temp / sPvector.size();
            else
                cov = 0;
            col.push_back(cov);
        }
        Covariance.push_back(col);
    }
    return Covariance;
}


float Node::computeDeterminant(vector<vector<float>> cvM){
    TMatrix tM;
    //bool allzero = 1;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            tM.m[j*3+i] = cvM[i][j];
            //cout << cvM[i][j] << " ";
        }
        //cout << endl;
    }
    //cout << "covarianceM det " << tM.determinant() << endl;
    return tM.determinant();
    
}

float Node::computeEntropy(vector<vector<float>> cvM){
    
    
    //cout << "det cvm " << computeDeterminant(cvM) << endl;
    if(computeDeterminant(cvM) != 0)
        return log2 (abs(computeDeterminant(cvM)));
    else
        return 0;
}

float Node::computeClassEntropy(vector<HPatch> patches){
    float p1,p2;
    p1 = classProbability(patches);
    p2 = 1 - p1;
    if(p1 == 0)
        return - p2*log2(p2);
    else if(p2 == 0)
        return - p1*log2(p1);
    else
        return - (p1*log2(p1)+p2*log2(p2));
    
}
//compute infomation gain
float Node::infoGain(vector<HPatch> parent, vector<HPatch> left, vector<HPatch> right){

    float IG = 0 ;
    IG = OkadaEntropy(parent) - float(left.size())/parent.size()*OkadaEntropy(left) - float(right.size())/parent.size()*OkadaEntropy(right);
    //cout << IG << endl;
    return IG;
}

float Node::OkadaEntropy(vector<HPatch> patches){
    vector<vector<float>> p1 = computeCovariance(patches,1);
    vector<vector<float>> p2 = computeCovariance(patches,0);
    //cout << "regression " << computeEntropy(p1)+computeEntropy(p2) << endl;
    //cout << "classification " << computeClassEntropy(patches) << endl;
    return computeEntropy(p1)+computeEntropy(p2)+max(classProbability(patches)-T_P,0.0)*computeClassEntropy(patches);
}


float traceCovariance(vector<vector<float>> cv){
    float t = 0;
    for(int i = 0; i < cv.size(); i++){
        t = t + cv[i][i];
        //cout << "cv i i " << cv[i][i] << endl;
    }
    return t;
}

void Node::computePositiveP(){
    float count = 0;
    for(int i = 0; i < beforeSplit.size(); i++ ){
        if(beforeSplit[i].positive == 0)
            count++;
    }
    positiveP = count/beforeSplit.size();
    
}

float Node::classProbability(vector<HPatch> patches){
    float count = 0;
    for(int i = 0; i < patches.size(); i++ ){
        if(patches[i].positive == 0)
            count++;
    }
    if(patches.size() == 0)
        return 0.0;
    return count/patches.size();
}

float Node::classUncertainty(vector<HPatch> patches){
    float entropy2d = 0;
    for(int i = 0; i < 17; i ++){
        float sum = 0;
        for(int j =0; j < patches.size(); j++){
            float c = exp(-1/TD_LAMBDA * findModulas(patches[j].groundT2d[i]));
            sum = sum + c;
        }
        entropy2d = entropy2d + sum/patches.size()*log2(sum/patches.size());
    }
    return -entropy2d;
}

float Node::infoGain2d(vector<HPatch> parent, vector<HPatch> left, vector<HPatch> right){
    float IG = 0;
    IG = classUncertainty(parent) - float(left.size())/parent.size()*classUncertainty(left) - float(right.size())/parent.size()*classUncertainty(right);
    return IG;
}
float Node::findModulas(vector<float> f){
    float x = 0;
    for(int i = 0; i < f.size(); i++ ){
        x= x + f[i]*f[i];
    }
    return sqrt(x);
}
//load tree
void DTree::read_tree(const string& fname){
    ifstream fInp;
    
    int tempDepth;
    vector<sub_patch> tempBestF;
    sub_patch tempSubpatch;
    int tempBest_T;
    bool tempIsLeaf;
    vector<float> tempMean;
    float tempData,tempTrace,tempP;
    fInp.open(fname);
    int nodeCounter = 0;
    int d = -1;
    while (fInp.eof()!=1){
        Node tempNode;
        fInp >> tempDepth >> tempIsLeaf;
        if (d == tempDepth) {
            nodeCounter++;
            nodesAtEachLevel[d] = nodeCounter;
        }
        else {
            nodeCounter = 1;
            nodesAtEachLevel.push_back(nodeCounter);
            d = tempDepth;
        }
        if(tempIsLeaf == 0){
            for(int i = 0; i < 2; i ++){
                fInp >> tempSubpatch.x >> tempSubpatch.y >> tempSubpatch.w >> tempSubpatch.h;
                //cout << tempSubpatch.x << endl;
                tempBestF.push_back(tempSubpatch);
            }
            fInp >> tempBest_T ;
            tempNode.best_T = tempBest_T;
            //cout << tempBestF[0].x << endl;
            //cout << tempBestF.size() << endl;
            tempNode.bestF = tempBestF;
            //cout << tempNode.bestF[0].x << endl;
            //cout << tempNode.bestF.size() << endl;
            tempNode.bestFArray = new sub_patch[tempBestF.size()];
            for(int bf =0; bf < tempBestF.size(); bf++){
                tempNode.bestFArray[bf] = tempBestF[bf];
            }
            tempBestF.clear();
        }
        else{
//            for(int i = 0; i < 6; i++){
//                fInp >> tempData;
//                //cout << tempData <<endl;
//                tempMean.push_back(tempData);
//            }
//            fInp >> tempTrace >> tempP;
//            tempNode.meanVector = tempMean;
//            tempNode.trace = tempTrace;
//            tempNode.positiveP = tempP;
//            tempMean.clear();
            for(int i = 0; i < 17; i++){
                fInp >> tempData;
                tempMean.push_back(tempData);
                fInp >> tempData;
                tempMean.push_back(tempData);
                fInp >> tempTrace;
                tempNode.trace2d.push_back(tempTrace);
                
            }
            fInp >> tempP;
            tempNode.positiveP = tempP;
            tempNode.meanVector = tempMean;
            tempMean.clear();
        }
        tempNode.depth = tempDepth; tempNode.isleaf = tempIsLeaf;
        treeTable.push_back(tempNode);
       // cout << "temp depth inside function: " << tempNode.depth << endl;
        //cout << treeTable.size() << endl;
        //cout << "table depth inside function: " << treeTable[t].bestF[0].x << endl;
        
        
    }
    
    fInp.close();
    noNodes = treeTable.size();
    treeTableArray = new Node[treeTable.size()];
    for(int i = 0; i < treeTable.size(); i++){
        treeTableArray[i] = treeTable[i];
    }
    
    nodesAtEachLevelArray = new int[nodesAtEachLevel.size()];
    for(int j = 0; j < nodesAtEachLevel.size(); j++){
        nodesAtEachLevelArray[j] = nodesAtEachLevel[j];
        
    }
}

//store tree
void DTree::write_tree(const string& fname){
    ofstream fOut;
    fOut.open(fname);
    for (int j = 0; j < treeTable.size(); j++){
        fOut << treeTable[j].depth << " " << treeTable[j].isleaf << " ";
        if(treeTable[j].isleaf == 0){
            for(int i = 0; i < 2; i++){
                fOut << treeTable[j].bestF[i].x << " " << treeTable[j].bestF[i].y << " " << treeTable[j].bestF[i].w << " " << treeTable[j].bestF[i].h << " " ;
            }
            fOut<< treeTable[j].best_T << endl;
        }else {
            for(int num = 0; num < 17; num++){
                vector<float> mean = treeTable[j].computeMeanVector2d(treeTable[j].beforeSplit,num);
                vector<vector<float>> cv = treeTable[j].computeCovariance2d(treeTable[j].beforeSplit,num);
                fOut<< mean[0] << " " << mean[1] <<" " << traceCovariance(cv) << " " ;
            }
            fOut << treeTable[j].positiveP << endl;
            //cout<< mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << " " << mean[4] << " " << mean[5] << " " << traceCovariance(cv) << " " << treeTable[j].positiveP <<  endl;
        }
    }
    fOut.close();
}



//build the tree
void DTree::growTree(vector<HPatch> PS, Mat* dImage){
    treeTable = m_root.makeTreeNoRecursion(PS,dImage);
}




//void DTree::regressionEstimation(Mat test3D,boundingBox testBbox,vector<float> testGt,vector<vector<float>>& estimatedMean,Mat img3D,vector<Vote>& votes,PatchSet testPS){
////    PatchSet testPS(SAMPLE_PATCHES);
////    testPS.sampleTestPatches(testBbox, img3D);
//    //testPS.getRandomPatches(testBbox, img3D, testGt, 1);
//    //cout << testPS.pSet.size() << endl;
//    vector<int> accumlativeNodesAtEachLevel;
//    int accumlativeSum = 0;
//    for(int i = 0;  i < nodesAtEachLevel.size(); i++) {
//        
//        
//        accumlativeNodesAtEachLevel.push_back(accumlativeSum);
//        accumlativeSum += nodesAtEachLevelArray[i];
//        //cout << " nodesAtEachLevel " << nodesAtEachLevel[i] << endl;
//        //cout << " accumlativeNodesAtEachLevel " << accumlativeNodesAtEachLevel[i] << endl;
//    }
//    
//    for(int i = 0; i < testPS.pSet.size(); i++){
//        int chosenNode = 0;
//        vector<float> tempMean;
//        //cout << endl << "considering patch " << i << "..." << endl;
//        for(int j = 0; j < nodesAtEachLevel.size(); j++){
//            //cout << " j " << j << endl;
//            //bool chosen = 0;
//            int noLeafNodes = 0;
//            int noNonLeafNodes = 0;
//            //cout << nodesAtEachLevel[j] << endl;
//            for(int k = 0; k < nodesAtEachLevelArray[j]; k++){
//                bool outputMeanVector = false;
//                //cout << " k " << k << endl;
//                int nodePositionInTreeTable = accumlativeNodesAtEachLevel[j] + k;
//                if (treeTableArray[nodePositionInTreeTable].isleaf) {
//                    noLeafNodes++;
//                }
//                else {
//                    noNonLeafNodes++;
//                }
//                //cout << " noLeafNodes " << noLeafNodes << endl;
//                //cout << " noNonLeafNodes " << noNonLeafNodes <<  endl;
//                if (k == chosenNode) {
//                    //string direction = "";
//                    if (treeTableArray[nodePositionInTreeTable].isleaf) {
//                        outputMeanVector = true;
//                        chosenNode = -1;
//                        //direction = "STOP";
//                    }
//                    else {
//                        testPS.pSetArray[i].chooseSubPatches(treeTableArray[nodePositionInTreeTable].bestF,test3D);
//                        //testPS.pSet[i].setSubPatchDistance(test3D);
//                        if (testPS.pSetArray[i].subPDistance > treeTableArray[nodePositionInTreeTable].best_T) {
//                            chosenNode = 2*(noNonLeafNodes-1);
//                            //direction = "LEFT";
//                            //chosen = 1;
//                            k = nodesAtEachLevel[j];
//                        }
//                        else {
//                            chosenNode = 2*noNonLeafNodes - 1;
//                            //direction = "RIGHT";
//                            //chosen = 1;
//                            k = nodesAtEachLevel[j];
//                        }
//                    }
//                    //cout << "at level: " << j << ", subPatch Distance = " << testPS.pSet[i].subPDistance << ", threshold = " << treeTable[nodePositionInTreeTable].best_T << ", chosenNode = " << chosenNode << " (" << direction << ")" << endl;
//                    if (outputMeanVector) {
//                        //cout << "stopped node position in tree table = " << nodePositionInTreeTable << endl;
//                        //for(int a = 0; a < treeTable[nodePositionInTreeTable].meanVector.size(); a++)
//                        //cout << treeTable[nodePositionInTreeTable].meanVector[a] << " ";
//                        //cout << endl;
//                        //cout << testPS.pSet[i].pC.p.d <<endl;
//                        //if(treeTable[nodePositionInTreeTable].trace < 400 ){
//                        
//                        if(testPS.pSetArray[i].pC.p.d != 0 && treeTableArray[nodePositionInTreeTable].trace < TRACE_MAX && testPS.pSetArray[i].subPDistance != 0 && treeTableArray[nodePositionInTreeTable].positiveP >= P_TH ){
//                            Vote v;
//                            //cout << "trace " << treeTable[nodePositionInTreeTable].trace << endl;
//                            tempMean = treeTableArray[nodePositionInTreeTable].meanVector;
//                            //cout << testPS.pSet[i].pC.p.x<< endl;
//                            tempMean[0] = tempMean[0] + testPS.pSetArray[i].pC.p.x;
//                            tempMean[1] = tempMean[1] + testPS.pSetArray[i].pC.p.y;
//                            tempMean[2] = tempMean[2] + testPS.pSetArray[i].pC.p.d;
////                            v.vote[0] = tempMean[0];
////                            v.vote[1] = tempMean[1];
////                            v.vote[2] = tempMean[2];
////                            v.vote[3] = tempMean[3];
////                            v.vote[4] = tempMean[4];
////                            v.vote[5] = tempMean[5];
//                            //v.trace
//                            //cout << "vote used: " << tempMean[0] << " " << tempMean[1] << " " << tempMean[2] << endl;
//                            //votes.push_back(v);
//                            estimatedMean.push_back(tempMean);
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//}

void DTree::regressionEstimation(Mat test3D,boundingBox testBbox,vector<float> testGt,vector<vector<float>>& estimatedMean,Mat img3D,vector<Vote>& votes,PatchSet testPS){
    //    PatchSet testPS(SAMPLE_PATCHES);
    //    testPS.sampleTestPatches(testBbox, img3D);
    //testPS.getRandomPatches(testBbox, img3D, testGt, 1);
    //cout << testPS.pSet.size() << endl;
    vector<int> accumlativeNodesAtEachLevel;
    int accumlativeSum = 0;
    for(int i = 0;  i < nodesAtEachLevel.size(); i++) {
        
        
        accumlativeNodesAtEachLevel.push_back(accumlativeSum);
        accumlativeSum += nodesAtEachLevelArray[i];
        //cout << " nodesAtEachLevel " << nodesAtEachLevel[i] << endl;
        //cout << " accumlativeNodesAtEachLevel " << accumlativeNodesAtEachLevel[i] << endl;
    }
    
    for(int x = 0; x < testBbox.width-80; x = x + step_x){
        for(int y = 0; y < testBbox.height-80; y = y + step_y){
            int p_x = testBbox.x + x;
            int p_y = testBbox.y + y;
//            cout << x << " " << y << endl;
//            cout << x << " " << y << endl;
            int pc_x = p_x + 40;
            int pc_y = p_y + 40;
            float subpatchDistance = 0;
            int chosenNode = 0;
            vector<float> tempMean;
            //cout << endl << "considering patch " << i << "..." << endl;
            for(int j = 0; j < nodesAtEachLevel.size(); j++){
                //cout << " j " << j << endl;
                //bool chosen = 0;
                int noLeafNodes = 0;
                int noNonLeafNodes = 0;
                //cout << nodesAtEachLevel[j] << endl;
                for(int k = 0; k < nodesAtEachLevelArray[j]; k++){
                    bool outputMeanVector = false;
                    //cout << " k " << k << endl;
                    int nodePositionInTreeTable = accumlativeNodesAtEachLevel[j] + k;
                    if (treeTableArray[nodePositionInTreeTable].isleaf) {
                        noLeafNodes++;
                    }
                    else {
                        noNonLeafNodes++;
                    }
                    //cout << " noLeafNodes " << noLeafNodes << endl;
                    //cout << " noNonLeafNodes " << noNonLeafNodes <<  endl;
                    if (k == chosenNode) {
                        //string direction = "";
                        if (treeTableArray[nodePositionInTreeTable].isleaf) {
                            outputMeanVector = true;
                            chosenNode = -1;
                            //direction = "STOP";
                        }
                        else {

                            sub_patch  rectangles[2];
//                            rectangles[0] = treeTableArray[nodePositionInTreeTable].bestF[0];
//                            rectangles[1] = treeTableArray[nodePositionInTreeTable].bestF[1];
                            for(int r = 0; r < 2; r++){
                                rectangles[r].x = treeTableArray[nodePositionInTreeTable].bestF[r].y + p_y;
                                rectangles[r].y = treeTableArray[nodePositionInTreeTable].bestF[r].x + p_x;
                                rectangles[r].w = treeTableArray[nodePositionInTreeTable].bestF[r].h;
                                rectangles[r].h = treeTableArray[nodePositionInTreeTable].bestF[r].w;                            }
                            

//                            cout << " 1 " << rectangles[0].x << " " << rectangles[0].y << " " << rectangles[0].w << " " << rectangles[0].h << endl;
//                            cout << " 2 " << rectangles[1].x << " " << rectangles[1].y << " " << rectangles[1].w << " " << rectangles[1].h << endl;
                            double sum1 = (test3D.at<double>(rectangles[0].x+rectangles[0].w,rectangles[0].y+rectangles[0].h) + test3D.at<double>(rectangles[0].x,rectangles[0].y) - test3D.at<double>(rectangles[0].x+rectangles[0].w,rectangles[0].y) - test3D.at<double>(rectangles[0].x,rectangles[0].y+rectangles[0].h))/(double)(rectangles[0].w*rectangles[0].h);
                            //cout << sum1 << endl;
                            double sum2 = (test3D.at<double>(rectangles[1].x+rectangles[1].w,rectangles[1].y+rectangles[1].h) + test3D.at<double>(rectangles[1].x,rectangles[1].y) - test3D.at<double>(rectangles[1].x+rectangles[1].w,rectangles[1].y) - test3D.at<double>(rectangles[1].x,rectangles[1].y+rectangles[1].h))/(double)(rectangles[1].w*rectangles[1].h);
                            subpatchDistance = sum1 - sum2;
                            if ( subpatchDistance > treeTableArray[nodePositionInTreeTable].best_T) {
                                chosenNode = 2*(noNonLeafNodes-1);
                                //direction = "LEFT";
                                //chosen = 1;
                                k = nodesAtEachLevel[j];
                            }
                            else {
                                chosenNode = 2*noNonLeafNodes - 1;
                                //direction = "RIGHT";
                                //chosen = 1;
                                k = nodesAtEachLevel[j];
                            }
                        }
                        //cout << "at level: " << j << ", subPatch Distance = " << testPS.pSet[i].subPDistance << ", threshold = " << treeTable[nodePositionInTreeTable].best_T << ", chosenNode = " << chosenNode << " (" << direction << ")" << endl;
                        if (outputMeanVector) {
                            //cout << "stopped node position in tree table = " << nodePositionInTreeTable << endl;
                            //for(int a = 0; a < treeTable[nodePositionInTreeTable].meanVector.size(); a++)
                            //cout << treeTable[nodePositionInTreeTable].meanVector[a] << " ";
                            //cout << endl;
                            //cout << testPS.pSet[i].pC.p.d <<endl;
                            //if(treeTable[nodePositionInTreeTable].trace < 400 ){
//                            cout << img3D.at<Vec3f>(pc_x,pc_y)[0] << endl;
//                            cout << img3D.at<Vec3f>(pc_x,pc_y)[1] << endl;
//                            cout << img3D.at<Vec3f>(pc_x,pc_y)[2] << endl;
                            if( img3D.at<Vec3f>(pc_y,pc_x)[2] != 0 && treeTableArray[nodePositionInTreeTable].trace < TRACE_MAX && subpatchDistance != 0 && treeTableArray[nodePositionInTreeTable].positiveP >= P_TH ){
                                //Vote v;
                                //cout << "trace " << treeTable[nodePositionInTreeTable].trace << endl;
                                tempMean = treeTableArray[nodePositionInTreeTable].meanVector;
                                //cout << testPS.pSet[i].pC.p.x<< endl;
                                tempMean[0] = tempMean[0] + img3D.at<Vec3f>(pc_y,pc_x)[0];
                                tempMean[1] = tempMean[1] + img3D.at<Vec3f>(pc_y,pc_x)[1];
                                tempMean[2] = tempMean[2] + img3D.at<Vec3f>(pc_y,pc_x)[2];
                                //                            v.vote[0] = tempMean[0];
                                //                            v.vote[1] = tempMean[1];
                                //                            v.vote[2] = tempMean[2];
                                //                            v.vote[3] = tempMean[3];
                                //                            v.vote[4] = tempMean[4];
                                //                            v.vote[5] = tempMean[5];
                                //v.trace
                                //cout << "vote used: " << tempMean[0] << " " << tempMean[1] << " " << tempMean[2] << endl;
                                //votes.push_back(v);
                                estimatedMean.push_back(tempMean);
                            }
                        }
                    }
                }
            }
        }
    }
    
}






void DTree::regressionEstimation2d(Mat test2d,boundingBox testBbox,vector<vector<float>> testGt,vector<vector<float>>& estimatedMean,vector<Vote>& votes, PatchSet testPS){
    //PatchSet testPS(SAMPLE_PATCHES);
    //testPS.sampleTestPatches2d(testBbox, test2d);
    //testPS.getRandomPatches(testBbox, img3D, testGt, 1);
    //cout << testPS.pSet.size() << endl;
    vector<int> accumlativeNodesAtEachLevel;
    int accumlativeSum = 0;
    for(int i = 0;  i < nodesAtEachLevel.size(); i++) {
        
        
        accumlativeNodesAtEachLevel.push_back(accumlativeSum);
        accumlativeSum += nodesAtEachLevelArray[i];
        //cout << " nodesAtEachLevel " << nodesAtEachLevel[i] << endl;
        //cout << " accumlativeNodesAtEachLevel " << accumlativeNodesAtEachLevel[i] << endl;
    }
    
    
    for(int i = 0; i < testPS.pSet.size(); i++){
        int chosenNode = 0;
        vector<float> tempMean;
        //cout << endl << "considering patch " << i << "..." << endl;
        for(int j = 0; j < nodesAtEachLevel.size(); j++){
            //cout << " j " << j << endl;
            //bool chosen = 0;
            int noLeafNodes = 0;
            int noNonLeafNodes = 0;
            //cout << nodesAtEachLevel[j] << endl;
            for(int k = 0; k < nodesAtEachLevelArray[j]; k++){
                bool outputMeanVector = false;
                //cout << " k " << k << endl;
                int nodePositionInTreeTable = accumlativeNodesAtEachLevel[j] + k;
                if (treeTableArray[nodePositionInTreeTable].isleaf) {
                    noLeafNodes++;
                }
                else {
                    noNonLeafNodes++;
                }
                //cout << " noLeafNodes " << noLeafNodes << endl;
                //cout << " noNonLeafNodes " << noNonLeafNodes <<  endl;
                if (k == chosenNode) {
                    //string direction = "";
                    if (treeTableArray[nodePositionInTreeTable].isleaf) {
                        outputMeanVector = true;
                        chosenNode = -1;
                        //direction = "STOP";
                    }
                    else {
                        //clock_t t = clock();
                        testPS.pSet[i].chooseSubPatches(treeTableArray[nodePositionInTreeTable].bestF,test2d);
                        //testPS.pSet[i].chooseSubPatches(treeTable[nodePositionInTreeTable].bestF,test2d);
                        //t = clock() - t;
                        //cout << "time " << (float)t/CLOCKS_PER_SEC << endl;
                        //testPS.pSet[i].setSubPatchDistance(test3D);
                        if (testPS.pSet[i].subPDistance > treeTableArray[nodePositionInTreeTable].best_T) {
                            chosenNode = 2*(noNonLeafNodes-1);
                            //direction = "LEFT";
                            //chosen = 1;
                            k = nodesAtEachLevelArray[j];
                        }
                        else {
                            chosenNode = 2*noNonLeafNodes - 1;
                            //direction = "RIGHT";
                            //chosen = 1;
                            k = nodesAtEachLevelArray[j];
                        }
                    }
                    //cout << "at level: " << j << ", subPatch Distance = " << testPS.pSet[i].subPDistance << ", threshold = " << treeTable[nodePositionInTreeTable].best_T << ", chosenNode = " << chosenNode << " (" << direction << ")" << endl;
                    if (outputMeanVector) {
                        //cout << "stopped node position in tree table = " << nodePositionInTreeTable << endl;
                        //for(int a = 0; a < treeTable[nodePositionInTreeTable].meanVector.size(); a++)
                        //cout << treeTable[nodePositionInTreeTable].meanVector[a] << " ";
                        //cout << endl;
                        //cout << testPS.pSet[i].pC.p.d <<endl;
                        //if(treeTable[nodePositionInTreeTable].trace < 400 ){
                        
                        if(treeTableArray[nodePositionInTreeTable].positiveP >= P_TH ){
                            //Vote v;
                            //cout << "trace " << treeTable[nodePositionInTreeTable].trace << endl;
                            tempMean = treeTableArray[nodePositionInTreeTable].meanVector;
                            //cout << testPS.pSet[i].pC.p.x<< endl;
                            for(int lm = 0; lm < tempMean.size(); lm = lm +2){
//                                cout << tempMean[lm] << endl;
//                                cout << tempMean[lm+1] << endl;
//                                cout << testPS.pSet[i].pC.c<<endl;
//                                cout << testPS.pSet[i].pC.r <<endl;
                                tempMean[lm] = tempMean[lm] + (float)testPS.pSet[i].pC.c/640;
                                tempMean[lm+1] = tempMean[lm+1] + (float)testPS.pSet[i].pC.r/480;
                            }

                            //v.trace
                            //cout << "vote used: " << tempMean[0] << " " << tempMean[1] << " " << tempMean[2] << endl;
                            //votes.push_back(v);
                            estimatedMean.push_back(tempMean);
                        }
                    }
                }
            }
        }
    }
    
}

