//
//  DTree.h
//  fypFirstDraft
//
//  Created by callen on 11/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#ifndef __fypFirstDraft__DTree__
#define __fypFirstDraft__DTree__

#include <iostream>
#include <limits>
#include "HPatch.h"
#include "TMatrix.h"
#include <cmath>
#include <fstream>
#endif /* defined(__fypFirstDraft__DTree__) */

using namespace std;

class leafNode{

public:
    
    leafNode(){}
    ~leafNode(){}
    
    float pfg;
    
    vector<float> mean;
    
    float trace;
    
};


class Node{
    
public:
    
    Node(){
        //left_child = 0;
        //right_child = 0;
        depth = 0;
        isleaf = 0;
    }
    Node(int d){
        //cout << d << endl;
        depth = d + 1;
        isleaf = 0;
        //cout << depth << endl;
    }
    bool read(ifstream& fInp);
    void write(ofstream& fInp);
    void generateRandomSubPatches(const vector<HPatch>& PS);
    void generateRandomThreshold(const vector<HPatch>& PS);
    void setPatchSetBeforeSplit(vector<HPatch> PS);
    void findBestT(vector<HPatch> PS,vector<Mat> integralImage);
    void splitPatchSet(int bin_test);
    vector<float> computeMeanVector(vector<HPatch> sPvector);
    vector<vector<float>> computeCovariance(vector<HPatch> sPvector, bool angleOrNose);
    float computeDeterminant(vector<vector<float>> cvM);
    float computeEntropy(vector<vector<float>> cvM);
    float infoGain(vector<HPatch> total, vector<HPatch> l, vector<HPatch> r);
    float traceCovariance(vector<vector<float>>  cv);
    
    
    bool isleaf;
    int depth;
    Node *left_child, *right_child;
    int best_T;
    vector<HPatch> beforeSplit;
    vector<HPatch> leftSplit;
    vector<HPatch> rightSplit;
    vector<int> lIndex,rIndex;
    vector<vector<sub_patch>> rectangles;
    vector<sub_patch> bestF;
    vector<int> rT;
    float detConvariance;
    vector<float> meanVector;
    vector<Node> makeTreeNoRecursion(vector<HPatch> PS,vector<Mat> dImage);
    float trace;
};


class DTree{

public:
    DTree(){};
    DTree(int d){max_depth = d;};
    ~DTree() { };
    
    Node m_root;
    
    void read_tree(const string& fname);
    
    void write_tree(const string& fname);
    
    void loadPreProcessedData(const string& fname,vector<HPatch>& PS);
    
    void growTree(vector<HPatch> PS, vector<Mat> dImage);
    
    void regressionEstimation(vector<threeDPostCal> test3D,boundingBox testBbox,vector<float> testGt);
    
    vector<sub_patch> f2;
    vector<sub_patch> f1;
    vector<float> subD;
    vector<int> rt;
    int max_depth;
    
    vector<float> mean;
    vector<Node> treeTable;
    vector<int> nodesAtEachLevel;
    int noNodes;
};