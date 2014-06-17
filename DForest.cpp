//
//  DForest.cpp
//  fypFirstDraft
//
//  Created by callen on 31/3/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "DForest.h"
#include <ctime>
#define SIZE 800

string convertIntToString(int number)
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

void DForest::growForest(vector<HPatch> wholeDataSet, Mat* depthIntegral){
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
        treeFilename = "2dtreebiggest"+ convertIntToString(i) + ".txt";
        trees[i].write_tree(treeFilename);
    }
}

void DForest::loadTree(){
    cout << " inside load tree " << endl;
    string treeFilename = "";
    for(int i = 0; i < noTrees; i++){
        DTree tree(10);
        treeFilename = "2dtreesnew"+ convertIntToString(i) + ".txt";
        tree.read_tree(treeFilename);
        trees.push_back(tree);
        cout << "tree" << i << " load" << endl;
    }
}

void DForest::regressionEstimation(Mat test3D,boundingBox testBbox,vector<float> testGt,Mat img3D){
    PatchSet testPS(1000);
    testPS.sampleTestPatches(testBbox, img3D);
    //cout << sizeof(testPS.pSetArray) <<endl;
    //cout << testPS.pSet.size() <<endl;
    for(int i = 0; i < noTrees; i++){
        //cout << "doing tree " << i <<endl;
        trees[i].regressionEstimation(test3D, testBbox, testGt, estimatedMean, img3D,votes,testPS);
        //cout << "size mean after tree" << i+1 << " " << estimatedMean.size() << endl;
    }
    //cout << "votes : " << votes.size() << endl;
    //this->meanShift(1.f, 6.f, 5, 400);
}

void DForest::regressionEstimation2d(Mat test2D,boundingBox testBbox,vector<vector<float>> testGt){
    //clock_t t = clock();
    PatchSet testPS(1000);
    testPS.sampleTestPatches2d(testBbox, test2D);
    //t = clock() - t;
    //cout << "time " << (float)t/CLOCKS_PER_SEC << endl;
    for(int i = 0; i < noTrees; i++){
        //clock_t t = clock();
        trees[i].regressionEstimation2d(test2D, testBbox, testGt, estimatedMean,votes,testPS);
       // t = clock() - t;
       // cout << "time " << (float)t/CLOCKS_PER_SEC << endl;
        //cout << "size mean after tree" << i+1 << " " << estimatedMean.size() << endl;
    }
    //cout << "votes : " << votes.size() << endl;
    //this->meanShift(1.f, 6.f, 5, 400);
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