//
//  Initial_values.cpp
//  clMDM
//
//  Created by Johan Van de Koppel on 10/09/2015.
//  Copyright (c) 2015 Johan Van de Koppel. All rights reserved.
//

#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
// Host initialization procedures
////////////////////////////////////////////////////////////////////////////////

void initialize_Cells(float *h_C,float *h_R,int hite,int width)
{
    int i,j;
    for(i=0;i<hite;i++)
    {
        for(j=0;j<width;j++)
        {
            h_C[i*width+j]=(float) (1+((h_R[i*width+j]<Frac)?1:0)-2*((h_R[i*width+j]<DFrac)?1:0));
        }
    }
}

void initialize_RandomCellMatrix(float* h_R,int hight,int width)
{
    int i,j;
    for(i=0;i<hight;i++)
    {
        for(j=0;j<width;j++)
        {
            h_R[i*width+j]=(float)rand() / (float)(RAND_MAX);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Prints the model name and additional info
////////////////////////////////////////////////////////////////////////////////

void Print_Label()
{
    //system("clear");
    printf("\n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n");
    printf(" * Mussel disturbance model                              * \n");
    printf(" * OpenCL implementation : Johan van de Koppel, 2014     * \n");
    printf(" * Following a model by Guichard et al 2003              * \n");
    printf(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n\n");
    
} // Print_Label


