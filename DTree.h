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
#include "HPatch.h"
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
        test_set = false;
        leaf_set = false;
        left_child = 0;
        right_child = 0;
        depth = 0;
    }
    bool test_set, leaf_set;
    int depth;
    Node *left_child, *right_child;
    bool read(ifstream& fInp);
    void write(ofstream& fInp);
    
};


class DTree{

public:
    DTree(){m_root = 0;};
    DTree(int d){m_root = 0; max_depth = d;};
    ~DTree() { delete m_root;};
    
    Node *m_root;
    
    bool read_tree(const string& fname);
    
    bool write_tree(const string& fname);
    
    float infoGain(vector<HPatch>, const vector<float>& randomThreshold);
    
    vector<int> generateRandomThreshold();
    
    void growTree();
    
    int max_depth;

};