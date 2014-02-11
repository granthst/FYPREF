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

//load tree
bool DTree::read_tree(const string& fname){
    return 0;
}

//store tree
bool DTree::write_tree(const string& fname){
    return 0;
}

//compute infomation gain
float DTree::infoGain(vector<HPatch>, const vector<float>& randomThreshold){
    return 0.1;
}

//generate a set of random Thresholds for binary test
vector<int> DTree::generateRandomThreshold(){
    
    vector<int> rT;
    return rT;
}

//build the tree
void DTree::growTree(){
    
}