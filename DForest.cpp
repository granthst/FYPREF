//
//  DForest.cpp
//  fypFirstDraft
//
//  Created by callen on 31/3/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "DForest.h"
#define SIZE 100

string convertIntToString(int number)
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

void DForest::growForest(vector<HPatch> wholeDataSet, vector<Mat> depthIntegral){
    for(int i = 0; i < noTrees; i++){
        cout << "growing tree" << i << endl;
        DTree tree(10);
        tree.growTree(generateSubSet(wholeDataSet,SIZE),depthIntegral);
        trees.push_back(tree);
        cout << "tree" << i << " done" << endl;
    }
}

void DForest::writeForest(){
    string treeFilename = "";
    for(int i = 0; i < noTrees; i++){
        treeFilename = "treeNew"+ convertIntToString(i) + ".txt";
        trees[i].write_tree(treeFilename);
    }
}

void DForest::loadTree(){
    string treeFilename = "";
    for(int i = 0; i < noTrees; i++){
        DTree tree(10);
        treeFilename = "treeNew"+ convertIntToString(i) + ".txt";
        tree.read_tree(treeFilename);
        trees.push_back(tree);
        cout << "tree" << i << " load" << endl;
    }
}

void DForest::regressionEstimation(Mat test3D,boundingBox testBbox,vector<float> testGt,Mat img3D){
    for(int i = 0; i < noTrees; i++){
        trees[i].regressionEstimation(test3D, testBbox, testGt, estimatedMean, img3D);
        //cout << "size mean after tree" << i+1 << " " << estimatedMean.size() << endl;
    }
    
}

vector<HPatch> DForest::generateSubSet(vector<HPatch> wholeDataSet, int size){
    vector<HPatch> subSet;
    int randomSample;
    for(int i = 0; i < size; i++){
        randomSample = rand()% (wholeDataSet.size()/40);
        //cout << " chosing image " << randomSample << endl;
        for(int j = 0; j < 40; j++){
            //cout << "patch number " << randomSample*40+j << endl;
            subSet.push_back(wholeDataSet[randomSample*40+j]);
        }
    }
    return subSet;
}