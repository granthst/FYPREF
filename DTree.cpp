//
//  DTree.cpp
//  fypFirstDraft
//
//  Created by callen on 11/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "DTree.h"

//read tree file
bool Node::read(ifstream& fInp){
    
    
    
    return 0;
}

//write tree file
void Node::write(ofstream& fInp){
    
}

//generate a set of random Thresholds for binary test
void Node::generateRandomThreshold(const vector<HPatch>& PS){
    
    
    int rt = 0;
    for (int i = 0; i < PS.size(); i++){
        rt = -800 + rand() % 1600;
        rT.push_back(rt);
    }
    
}

void Node::setPatchSetBeforeSplit(const vector<HPatch>& PS){

    beforeSplit = PS;
}

void Node::findBestT(const vector<HPatch>& PS){
    
    float infoGainTemp = numeric_limits<int>::min();
    cout << "int min " << infoGainTemp << endl;
    vector<HPatch> tempLeft, tempRight;
    for(int j = 0; j < rT.size(); j++){
        cout << "threshold : " << rT[j] << endl;
        for( int i = 0; i < PS.size(); i++){
            
            if(PS[i].subPDistance > rT[j])
                tempLeft.push_back(PS[i]);
            else
                tempRight.push_back(PS[i]);
        }
        //leftSplit = tempLeft;
        //rightSplit = tempRight;
        cout << j << " left child patch size : " << tempLeft.size() << " right child patch size : " << tempRight.size() << endl;
        //cout << "info gain " << infoGain(beforeSplit, tempLeft, tempRight, rT) << endl;
        
        if (infoGainTemp < infoGain(beforeSplit, tempLeft, tempRight, rT)){
            leftSplit = tempLeft;
            rightSplit = tempRight;
            infoGainTemp = infoGain(beforeSplit, tempLeft, tempRight, rT);
            
        }
        tempLeft.clear();
        tempRight.clear();
    }
            cout << " left child patch size : " << leftSplit.size() << " right child patch size : " << rightSplit.size() << endl;
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

vector<vector<float>> Node::computeCovariance(vector<HPatch> sPvector ){
    vector<vector<float>> Covariance;
    vector<float> meanVector = computeMeanVector(sPvector);
    for(int i = 0; i < 6; i++){
        vector<float> col;
        for(int j = 0; j < 6; j++){
            float temp = 0.0;
            for(int k = 0; k < sPvector.size(); k++){
                temp = temp + sPvector[k].groundT[i]*sPvector[k].groundT[j];
            }
            float cov = temp - (meanVector[i]*meanVector[j]);
            col.push_back(cov);
        }
        Covariance.push_back(col);
    }
    return  Covariance;
    //cv = Covariance;
}

float Node::computeDeterminant(vector<vector<float>> cvM){
    TMatrix tM;
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            tM.m[j*6+i] = cvM[i][j];
        }
    }
    //cout << "covarianceM det " << tM.determinant() << endl;
    return tM.determinant();
    
}

float Node::computeEntropy(vector<vector<float>> cvM){
    
    float entropy = log2 (computeDeterminant(cvM));
    //cout << "entropy " << entropy << endl;
    return entropy;
}

//compute infomation gain
float Node::infoGain(vector<HPatch> parent, vector<HPatch> left, vector<HPatch> right, const vector<int>& randomThreshold){
    vector<vector<float>> p = computeCovariance(parent);
    vector<vector<float>> l = computeCovariance(left);
    vector<vector<float>> r = computeCovariance(right);
    float IG = computeEntropy(p) - float(left.size())/parent.size()*computeEntropy(l) - float(right.size())/parent.size()*computeEntropy(r);
    //cout << IG << endl;
    return IG;
}

//load tree
bool DTree::read_tree(const string& fname){
    return 0;
}

//store tree
bool DTree::write_tree(const string& fname){
    return 0;
}











//build the tree
void DTree::growTree(){
    
}