//
//  DTree.cpp
//  fypFirstDraft
//
//  Created by callen on 11/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "DTree.h"
#define n_P 800
#define THESHOLD_MAX 1000
#define DEPTH_TREE 5
string integralImageFilename = "integral1.txt";
string treeFilename = "tree.txt";
//read tree file
bool Node::read(ifstream& fInp){
    
    
    
    return 0;
}

//write tree file
void Node::write(ofstream& fInp){
    
}

void Node::generateRandomSubPatches(){
    sub_patch temp;
    vector<sub_patch> tempSub;
    for (int j = 0; j < n_P; j++){
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
    for (int i = 0; i < PS.size(); i++){
        rt = -THESHOLD_MAX + rand() % (2*THESHOLD_MAX);
        rT.push_back(rt);
    }
    
}

void Node::setPatchSetBeforeSplit( vector<HPatch> PS){

    beforeSplit = PS;
}

void Node::findBestT(vector<HPatch> PS, vector<threeDPostCal> dImage){
    
    //cout << PS.size() << endl;
    float infoGainTemp = numeric_limits<int>::min();
    //cout << "int min " << infoGainTemp << endl;
    vector<HPatch> tempLeft, tempRight;
    int x1,y1,w1,h1,x2,y2,w2,h2;
    float tempSB = 0;
    float tempRT = 0;
    sub_patch tempS1,tempS2;
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

        for(int j = 0; j < n_P; j++){
            //cout << "threshold : " << rT[j] << endl;
            
            for( int i = 0; i < PS.size(); i++){
                //cout << " f1 x " << rectangles[j][0].x << endl << endl;
                //int rf = rand() % 399;
                //int rt = rand() % 399;
                //PS[i].chooseSubPatches(rectangles[j]);
                //cout << " f1 " << i << "x " << PS[i].rectangles[0].x << endl;
                //fOut << rectangles[j][0].x << " " << rectangles[j][0].y << " " << rectangles[j][0].w << " " << rectangles[j][0].h << endl;
                //fOut << rectangles[j][1].x << " " << rectangles[j][1].y << " " << rectangles[j][1].w << " " << rectangles[j][1].h << endl;
                //PS[i].setSubPatchDistance(dImage);
                //fOut << PS[i].subPDistance << " " << rT[j] << endl;
                /*if (i <= 1) {
                    cout << "patch no: " << i;
                    cout << ", integral image " << PS[i].subPDistance;
                    cout << ", threshold " << rT[j] << endl;
                }*/
                
                //fInp >> x1 >> y1 >> w1 >> h1 >> x2 >> y2 >> w2 >> h2;
                //fInp >> tempSB >> tempRT;
                //cout << j*n_P+i << endl;
                //cout << " index " << j*n_P+PS[i].index << endl;
                
                x1 = f1[j*n_P+PS[i].index].x ; y1 = f1[j*n_P+PS[i].index].y ; w1 = f1[j*n_P+PS[i].index].w; h1 = f1[j*n_P+PS[i].index].h;
                x2 = f2[j*n_P+PS[i].index].x; y2 = f2[j*n_P+PS[i].index].y; w2 = f2[j*n_P+PS[i].index].w; h2 = f2[j*n_P+PS[i].index].h;
                tempSB = subD[j*n_P+PS[i].index];
                tempRT = rt[j];

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
                    //cout << " j : " << j << " infogain " << infoGain(beforeSplit, tempLeft, tempRight);
                     //cout << " j " << j  <<" " << tempS1.x << " " << tempS1.y << " " << tempS1.w << " " << tempS1.h << " " << tempS2.x << " " << tempS2.y << " " << tempS2.w << " " << tempS2.h << " " << best_T << endl;
                }
            
            
            tempLeft.clear();
            tempRight.clear();

        }

    if(((leftSplit.size()==0) && (rightSplit.size()==0)) || depth == DEPTH_TREE){
        //cout << " leaf node " << endl;
        isleaf = 1;
    }
    //cout << "depth : " << depth << endl;
    //cout << " isleaf " << isleaf << endl;
        //cout << best_T << endl;
        //cout << "ig " << infoGainTemp << endl;


   bestF.push_back(tempS1);  bestF.push_back(tempS2);
    
    /*for(int i = 0; i < 2; i++){
        if(leftSplit.size()!= 0)
            fOut << bestF[i].x << " " << bestF[i].y << " " << bestF[i].w << " " << bestF[i].h << " " ;
    }
    fOut << best_T;
    fOut << endl;
    
    //fInp.close();
    
    fOut.close();*/
    
    //cout << " best threshold " << best_T << endl;
    //left_child->setPatchSetBeforeSplit(leftSplit);
    //left_child->rt = rt;
    //left_child->f1 = f1;
    //left_child->f2 = f2;
    //left_child->subD = subD;

    /*for(int i = 0; i < leftSplit.size(); i++)
        cout << leftSplit[i].index << endl;*/
        //fOut.close();
    //cout << "at depth " << depth << " beforesplit patch size " << beforeSplit.size() << endl;
    //cout << "at depth " << depth << " left child patch size : " << leftSplit.size() << " right child patch size : " << rightSplit.size() << endl;
//    if( depth < 2){
//        if(leftSplit.size() > 3){
//            left_child = new Node(f1,f2,subD,rt,depth);
//            //cout << " lf child depth " << left_child.depth << endl;
//            left_child->setPatchSetBeforeSplit(leftSplit);
//            //cout << left_child.beforeSplit.size() << endl;
//            left_child->findBestT(left_child->beforeSplit, dImage);
//        }
//        if(rightSplit.size() > 3){
//            right_child = new Node(f1,f2,subD,rt,depth);
//            //cout << " right_child child depth " << right_child.depth << endl;
//            right_child->setPatchSetBeforeSplit(rightSplit);
//            //cout << right_child.beforeSplit.size() << endl;
//            right_child->findBestT(right_child->beforeSplit, dImage);
//        }
        
//    }
//    else
//        return;

}

vector<Node> Node::makeTreeNoRecursion(vector<HPatch> PS, vector<threeDPostCal> dImage) {
    
    int iterationNo = 0;
    vector<Node> pointers;
    Node initialPointer = *this;
    initialPointer.setPatchSetBeforeSplit(PS);
    initialPointer.findBestT(PS,dImage);
    pointers.push_back(initialPointer);
    int prevNoPointers = 0;
    while (iterationNo < DEPTH_TREE) {
        int noPointers = (int) pointers.size();
        
        //cout << "iterationNo = " << iterationNo << ", noPointers = " << noPointers << endl;
        
        for (int j = prevNoPointers; j < noPointers; j++) {
            Node parentPointer = pointers[j];
            if(parentPointer.isleaf == 0){
                Node leftPointer = Node(f1,f2,subD,rt,iterationNo);
                leftPointer.setPatchSetBeforeSplit(parentPointer.leftSplit);/***/
                leftPointer.findBestT(leftPointer.beforeSplit, dImage);
                leftPointer.depth = iterationNo+1;
                pointers.push_back(leftPointer);
                //pointers[j].left_child = leftPointer;
                
                Node rightPointer = Node(f1,f2,subD,rt,iterationNo);
                rightPointer.setPatchSetBeforeSplit(parentPointer.rightSplit);/***/
                rightPointer.depth = iterationNo+1;
                rightPointer.findBestT(rightPointer.beforeSplit, dImage);
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
    IG = (computeEntropy(p1)+ computeEntropy(p2)) - float(left.size())/parent.size()*(computeEntropy(l1)+computeEntropy(l2)) - float(right.size())/parent.size()*(computeEntropy(r1)+computeEntropy(r2));
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

//load tree
void DTree::read_tree(const string& fname){
    ifstream fInp;
    
    int tempDepth;
    vector<sub_patch> tempBestF;
    sub_patch tempSubpatch;
    int tempBest_T;
    bool tempIsLeaf;
    vector<float> tempMean;
    float tempData,tempTrace;
    fInp.open(fname);
    int t = 0;
    while (fInp.eof()!=1){
        Node tempNode;
        fInp >> tempDepth >> tempIsLeaf;
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
            fInp >> tempTrace;
            tempNode.meanVector = tempMean;
            tempNode.trace = tempTrace;
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
            fOut<< mean[0] << " " << mean[1] << " " << mean[2] << " " << mean[3] << " " << mean[4] << " " << mean[5] << " " << traceCovariance(cv) <<  endl;
        }
    }
    fOut.close();
}


void DTree::loadPreProcessedData(const string& fname){
    ifstream fInp;
    fInp.open(integralImageFilename);
    float tempSB = 0;
    float tempRT = 0;
    sub_patch temp;
    for(int j = 0; j < n_P; j++){
        for( int i = 0; i < n_P; i++){
            fInp >> temp.x >> temp.y >> temp.w >> temp.h;
            f1.push_back(temp);
            fInp >> temp.x >> temp.x >> temp.w >> temp.h;
            f2.push_back(temp);
            
            fInp >> tempSB >> tempRT;
            subD.push_back(tempSB);
        }
        rt.push_back(tempRT);
    }
    m_root.f1 = f1;
    m_root.f2 = f2;
    m_root.subD = subD;
    m_root.rt = rt;
}

//build the tree
void DTree::growTree(vector<HPatch> PS, vector<threeDPostCal> dImage){
    treeTable = m_root.makeTreeNoRecursion(PS,dImage);
}


void DTree::regressionEstimation(vector<threeDPostCal> test3D,boundingBox testBbox,vector<float> testGt){
    PatchSet testPS(40);
    testPS.getRandomPatches(testBbox, test3D, testGt);
    //cout << treeTable[0].best_T << endl;
    vector<float> estimatedMean;
    for(int i = 0; i < testPS.pSet.size(); i++){
        //for(int j = 0; j < treeTable.size(); j++){
            testPS.pSet[i].chooseSubPatches(treeTable[1].bestF);
            testPS.pSet[i].setSubPatchDistance(test3D);
            cout << "patch " << i << endl;
            //cout << treeTable[1].bestF[0].x << endl;
            cout << testPS.pSet[i].subPDistance << endl;
            cout << treeTable[0].best_T << endl;
        //}
    }
}