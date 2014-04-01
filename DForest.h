//
//  DForest.h
//  fypFirstDraft
//
//  Created by callen on 31/3/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#ifndef __fypFirstDraft__DForest__
#define __fypFirstDraft__DForest__

#include <iostream>
#include "DTree.h"
#endif /* defined(__fypFirstDraft__DForest__) */

class DForest{

public:
    DForest(){}
    DForest(int i){
        noTrees = i;
    }
    int noTrees;
    vector<DTree> trees;
    vector<vector<float>> estimatedMean;
    void growForest(vector<HPatch> wholeDataSet, vector<Mat> depthIntegral);
    void writeForest();
    void loadTree();
    void regressionEstimation(Mat test3D,boundingBox testBbox,vector<float> testGt,Mat img3D);
};