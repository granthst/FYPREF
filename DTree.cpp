//
//  DTree.cpp
//  fypFirstDraft
//
//  Created by callen on 11/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "DTree.h"
#define n_F 50
#define n_T 10
#define THESHOLD_MAX 1000
#define DEPTH_TREE 15
#define TRACE_MAX 400
#define SAMPLE_PATCHES 10000
#define P_TH 0.5
string integralImageFilename = "integral3.txt";
string treeFilename = "tree.txt";
//read tree file
bool Node::read(ifstream& fInp){
    
    
    
    return 0;
}

//write tree file
void Node::write(ofstream& fInp){
    
}

void Node::generateRandomSubPatches(const vector<HPatch>& PS){
    sub_patch temp;
    vector<sub_patch> tempSub;
    for (int j = 0; j < n_F; j++){
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

void Node::findBestT(vector<HPatch> PS, vector<Mat> itegralImage){
    
    //cout << PS.size() << endl;
    float infoGainTemp = numeric_limits<int>::min();
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
        for(int j = 0; j < rectangles.size(); j++){
            
            //cout << "threshold : " << k << " : " << rT[k] << " for test " << j << endl;
            
            for( int i = 0; i < PS.size(); i++){
                //cout << " f1 x " << rectangles[j][0].x << endl << endl;
                //int rf = rand() % 399;
                //int rt = rand() % 399;
                PS[i].chooseSubPatches(rectangles[j],itegralImage[PS[i].index]);
                x1 = rectangles[j][0].x; x2 = rectangles[j][1].x;
                y1 = rectangles[j][0].y; y2 = rectangles[j][1].y;
                w1 = rectangles[j][0].w; w2 = rectangles[j][1].w;
                h1 = rectangles[j][0].h; h2 = rectangles[j][1].h;
                
                /*if(w1 == 0)
                    w1 = 1;
                if(h1 == 0)
                    h1 = 1;
                if(h2 == 0)
                    h2 = 1;
                if(w2 == 0)
                    w2 = 1;*/
                /*x1 = f1[j*n_P+PS[i].index].x ; y1 = f1[j*n_P+PS[i].index].y ; w1 = f1[j*n_P+PS[i].index].w; h1 = f1[j*n_P+PS[i].index].h;
                x2 = f2[j*n_P+PS[i].index].x; y2 = f2[j*n_P+PS[i].index].y; w2 = f2[j*n_P+PS[i].index].w; h2 = f2[j*n_P+PS[i].index].h;
                tempSB = subD[j*n_P+PS[i].index];*/
                
                //cout << PS[i].index << endl;
                tempSB = PS[i].subPDistance;
                //cout << "tempSB : " << tempSB << endl;
                //cout << "tempRT : " << tempRT << endl;
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
            //cout << "info gain " << infoGain(beforeSplit, tempLeft, tempRight) << endl;
            
           //cout << "rt " << tempRT << endl;
            //cout << " tempsb " << tempSB << endl;

                if (infoGainTemp < infoGain(beforeSplit, tempLeft, tempRight)){
                    leftSplit = tempLeft;
                    rightSplit = tempRight;
                    infoGainTemp = infoGain(beforeSplit, tempLeft, tempRight);
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

vector<Node> Node::makeTreeNoRecursion(vector<HPatch> PS, vector<Mat> integralImage) {
    
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
            float cov = temp / sPvector.size();
            col.push_back(cov);
        }
        Covariance.push_back(col);
    }

    return  Covariance;
    //cv = Covariance;
}

float Node::computeDeterminant(vector<vector<float>> cvM){
    TMatrix tM;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            tM.m[j*3+i] = cvM[i][j];
        }
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

//compute infomation gain
float Node::infoGain(vector<HPatch> parent, vector<HPatch> left, vector<HPatch> right){
    vector<vector<float>> p1 = computeCovariance(parent,1);
    vector<vector<float>> l1 = computeCovariance(left,1);
    vector<vector<float>> r1 = computeCovariance(right,1);
    vector<vector<float>> p2 = computeCovariance(parent,0);
    vector<vector<float>> l2 = computeCovariance(left,0);
    vector<vector<float>> r2 = computeCovariance(right,0);
    float IG = 0 ;
    IG = (computeEntropy(p1) + computeEntropy(p2)) - float(left.size())/parent.size()*(computeEntropy(l1)+computeEntropy(l2)) - float(right.size())/parent.size()*(computeEntropy(r1)+computeEntropy(r2));
    //cout << IG << endl;
    return IG;
}


float traceCovariance(vector<vector<float>> cv){
    float t = 0;
    for(int i = 0; i < 3; i++){
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
            tempBestF.clear();
        }
        else{
            for(int i = 0; i < 6; i++){
                fInp >> tempData;
                tempMean.push_back(tempData);
            }
            fInp >> tempTrace >> tempP;
            tempNode.meanVector = tempMean;
            tempNode.trace = tempTrace;
            tempNode.positiveP = tempP;
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
            vector<float> mean = treeTable[j].computeMeanVector(treeTable[j].beforeSplit);
            vector<vector<float>> cv = treeTable[j].computeCovariance(treeTable[j].beforeSplit,1);
            fOut<< mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << " " << mean[4] << " " << mean[5] << " " << traceCovariance(cv) << " " << treeTable[j].positiveP <<  endl;
        }
    }
    fOut.close();
}



//build the tree
void DTree::growTree(vector<HPatch> PS, vector<Mat> dImage){
    treeTable = m_root.makeTreeNoRecursion(PS,dImage);
}


void DTree::regressionEstimation(Mat test3D,boundingBox testBbox,vector<float> testGt,vector<vector<float>>& estimatedMean,Mat img3D){
    PatchSet testPS(SAMPLE_PATCHES);
    testPS.getRandomPatches(testBbox, img3D, testGt,0);
    //cout << testPS.pSet.size() << endl;
    vector<int> accumlativeNodesAtEachLevel;
    int accumlativeSum = 0;
    for(int i = 0;  i < nodesAtEachLevel.size(); i++) {
        
        
        accumlativeNodesAtEachLevel.push_back(accumlativeSum);
        accumlativeSum += nodesAtEachLevel[i];
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
            for(int k = 0; k < nodesAtEachLevel[j]; k++){
                bool outputMeanVector = false;
                //cout << " k " << k << endl;
                int nodePositionInTreeTable = accumlativeNodesAtEachLevel[j] + k;
                if (treeTable[nodePositionInTreeTable].isleaf) {
                    noLeafNodes++;
                }
                else {
                    noNonLeafNodes++;
                }
                //cout << " noLeafNodes " << noLeafNodes << endl;
                //cout << " noNonLeafNodes " << noNonLeafNodes <<  endl;
                if (k == chosenNode) {
                    string direction = "";
                    if (treeTable[nodePositionInTreeTable].isleaf) {
                        outputMeanVector = true;
                        chosenNode = -1;
                        direction = "STOP";
                    }
                    else {
                        testPS.pSet[i].chooseSubPatches(treeTable[nodePositionInTreeTable].bestF,test3D);
                        //testPS.pSet[i].setSubPatchDistance(test3D);
                        if (testPS.pSet[i].subPDistance > treeTable[nodePositionInTreeTable].best_T) {
                            chosenNode = 2*(noNonLeafNodes-1);
                            direction = "LEFT";
                            //chosen = 1;
                            k = nodesAtEachLevel[j];
                        }
                        else {
                            chosenNode = 2*noNonLeafNodes - 1;
                            direction = "RIGHT";
                            //chosen = 1;
                            k = nodesAtEachLevel[j];
                        }
                    }
                    //cout << "at level: " << j << ", subPatch Distance = " << testPS.pSet[i].subPDistance << ", threshold = " << treeTable[nodePositionInTreeTable].best_T << ", chosenNode = " << chosenNode << " (" << direction << ")" << endl;
                    if (outputMeanVector) {
                        //cout << "stopped node position in tree table = " << nodePositionInTreeTable << ", mean vector = ";
                        for(int a = 0; a < treeTable[nodePositionInTreeTable].meanVector.size(); a++)
                            //cout << treeTable[nodePositionInTreeTable].meanVector[a] << " ";
                        //cout << endl;
                        //cout << testPS.pSet[i].pC.p.d <<endl;
                        //if(treeTable[nodePositionInTreeTable].trace < 400 ){
                        
                        if(testPS.pSet[i].pC.p.d != 0 && treeTable[nodePositionInTreeTable].trace < TRACE_MAX && testPS.pSet[i].subPDistance != 0 && treeTable[nodePositionInTreeTable].positiveP >= P_TH ){
                            
                            //cout << "trace " << treeTable[nodePositionInTreeTable].trace << endl;
                            tempMean = treeTable[nodePositionInTreeTable].meanVector;
                            //cout << testPS.pSet[i].pC.p.x << endl;
                            tempMean[0] = tempMean[0] + testPS.pSet[i].pC.p.x;
                            tempMean[1] = tempMean[1] + testPS.pSet[i].pC.p.y;
                            tempMean[2] = tempMean[2] + testPS.pSet[i].pC.p.d;
                            estimatedMean.push_back(tempMean);
                        }
                    }
                }
            }
        }
    }
}