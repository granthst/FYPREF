//
//  Matrix.cpp
//  fypFirstDraft
//
//  Created by callen on 13/2/14.
//  Copyright (c) 2014 callen. All rights reserved.
//

#include "TMatrix.h"

float TMatrix::at(int r, int c)
{
	return m[r*n+c];
}

void TMatrix::put(float v, int r, int c)
{
	m[r*n+c]=v;
}

void TMatrix::swapRows(int r1,int r2)
{
    
	float temp[n];
	int i,j,k;
	for( i=0;i<n;i++)
	{
		temp[i]=m[r1*n+i];
	}
	for( j=0;j<n;j++)
	{
        
		m[r1*n+j]=m[r2*n+j];
	}
	for( k=0;k<n;k++)
	{
        
        
		m[r2*n+k]=temp[k];
	}
    
}

int TMatrix::searchRow(int row,int col)
{
    int a;
    for( a=row;a<n;a++)
    {
        if(m[a*n+col]!=0)
            return a;
    }
    
    return row;
    
}



void TMatrix::UpperTMatrix()
{
	
	int col=0;
	int row=0;
	int row1=0;
	float factor=0;
	float inter=0;
	
	
	for(col=0;col<n;col++)
	{
        if(at(col,col)==0&&searchRow(col,col)!=row1)
        {
            swapRows(col,searchRow(col,col));
            swaped=true;
			
			for(row=0;row<n-1;row++)
            {
                if(row+1>col)
                {
                    if(at(row+1,col)!=0)
                    {
                        factor=at(row+1,col)/at(col,col);
                        int i;
                        for(i=col;i<n;i++)
                        {
                            inter=at(row+1,i)-factor*at(col,i);
                            put(inter,row+1,i);
                        }
                    }
                }
            }
			
        }
        
        if(at(col,col)!=0)
        {
			for(row=0;row<n-1;row++)
            {
                if(row+1>col)
                {
                    if(at(row+1,col)!=0)
                    {
                        factor=at(row+1,col)/at(col,col);
                        int i;
                        for( i=col;i<n;i++)
                        {
                            inter=at(row+1,i)-factor*at(col,i);
                            put(inter,row+1,i);
                        }
                    }
                }
            }	
			
            
        }
        
        
        
        
	}	
	
}



float TMatrix::determinant()
{	
    
    if(at(0,0)==0&&searchRow(0,0)==0)
        return 0.0;
    else
    {
        
		int i=0;
		UpperTMatrix();
		float d=1;
		while(i<(n*n))
		{	
			//cout<<m[i]<<endl;
			d=d*m[i]; 
			i=i+(n+1);
		}
		if(swaped)
            return -d;
		return d;
    }
    
	
}