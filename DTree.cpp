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

void Node::setPatchSetBeforeSplit(const vector<HPatch>& PS){
    cout << " size of the data set " << PS.size() << endl;
    for (int j = 0; j < PS.size(); j++){
        //cout << " Patch " << j << " of image number : " << i << endl;
        cout << " px, py : " << PS[j].p_x << " , " << PS[j].p_y << endl;
        //cout << " f2 : " << PS[j].rectangles[1].x << " , " << PS[j].rectangles[1].y << endl;
        cout << " f1 : " << PS[j].rectangles[0].x << " , " << PS[j].rectangles[0].y << endl;
        cout << " f2 : " << PS[j].rectangles[1].x << " , " << PS[j].rectangles[1].y << endl;
    }
}

void Node::findBestT(){
    
}

void Node::splitPatchSet(int bin_test){
    
}



//load tree
bool DTree::read_tree(const string& fname){
    return 0;
}

//store tree
bool DTree::write_tree(const string& fname){
    return 0;
}

//compute infomation gain
float DTree::infoGain(vector<HPatch> patchSet, const vector<float>& randomThreshold){
    return 0.1;
}

//generate a set of random Thresholds for binary test
void DTree::generateRandomThreshold(int n){
    
    
    int rt = 0;
    for (int i = 0; i < n; i++){
        rt = -800 + rand() % 1600;
        rT.push_back(rt);
    }

}






//build the tree
void DTree::growTree(){
    
}