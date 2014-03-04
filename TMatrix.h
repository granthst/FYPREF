//
//  Matrix.h
//  fypFirstDraft
//
//  Created by callen on 13/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#ifndef __fypFirstDraft__Matrix__
#define __fypFirstDraft__Matrix__

#include <iostream>
#define n 3
#endif /* defined(__fypFirstDraft__Matrix__) */
class TMatrix{
    
public:
    TMatrix(){}
    float at(int r, int c);
    void put(float v, int r, int c);
    void swapRows(int r1,int r2);
    int searchRow(int row,int col);
    void UpperTMatrix();
    float determinant();
    float m[9];
    bool swaped = false;
};