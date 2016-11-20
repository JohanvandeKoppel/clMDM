//
//  MDM_Kernel.cl
//  clMDM
//
//  Created by Johan Van de Koppel on 11-09-14.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

//#include <stdio.h>
#include "Settings_and_Parameters.h"

//Disturbed NeighborHood Matrix
constant int	DNH[9]	=	{ 1,1,1,
                              1,0,1,
                              1,1,1};
//Growth NeighborHood Matrix
constant int	GNH[9]	=	{ 1,1,1,
                              1,0,1,
                              1,1,1};

////////////////////////////////////////////////////////////////////////////////
// GPU Device procedures, used within the simulation kernel
////////////////////////////////////////////////////////////////////////////////

// This function checks if the cell within matrix "C", labelled "current", has "value" condition
void CheckCellValue(__global float *C, __global float *Cgen, int condition, int current)
{
    if(C[current]==condition)
        Cgen[current]=1;
    else
        Cgen[current]=0;
}

// This function counts the number of neighbours of the cell labelled "current"
float Neighbourhood_filter(constant int* transfer, __global float* matrix, int current)
{
	return (transfer[0]*matrix[current-1-Grid_Width]	+	transfer[1]*matrix[current-Grid_Width]	+	transfer[2]*matrix[current-Grid_Width+1]	+
			transfer[3]*matrix[current-1]           +	transfer[4]*matrix[current]          	+	transfer[5]*matrix[current+1]	        +
			transfer[6]*matrix[current-1+Grid_Width]	+	transfer[7]*matrix[current+Grid_Width]	+	transfer[8]*matrix[current+1+Grid_Width]);
}


////////////////////////////////////////////////////////////////////////////////
// GPU Simulation kernel
////////////////////////////////////////////////////////////////////////////////

__kernel void
MusselDisturbKernel (__global float* d_Cell, __global float* d_Rdom,
                     __global float* d_CDist, __global float* d_COccu)
{
	int Nu0;
	float N;
    
	float Colonization;
    float Erosion;
	float BecomeEmpty;
    
    size_t current = get_global_id(0);
	
	int row		= floor((float)current/(float)Grid_Width);
	int column	= current%Grid_Width;

    
    // First, two matrices called d_CDist and d_COccu are created
    // if the current location is occupied (equals to 2), d_COccu becomes 1 at this location
    CheckCellValue(d_Cell, d_COccu,2,current);
    
    // if the current location is disturbed (equals to 0), d_CDist becomes 1 at this location
	CheckCellValue(d_Cell, d_CDist,0,current);
    
    // The above two lines are done for the entire matrices before the kernel continues,
    // ensured by this command:
    
	barrier(CLK_GLOBAL_MEM_FENCE);
    
	// Excluding the edges when processing the kernel, as the neighborhood is undefined
    if(row > 0 && row < Grid_Width-1 && column > 0 && column < Grid_Width-1)
    {
        // Here the calculations for state transfers are made
        
        //calculating whether there are any disturbed neighbors
        if (Neighbourhood_filter(DNH, d_CDist,current)	>	0.0f)
            Nu0=1;
        else
            Nu0=0;
        
        //calculating the fraction of occupied neighbors
        N=(float)(Neighbourhood_filter(GNH, d_COccu,current))	/	8.0f;
		
        // Here, a possible transition from occupied to disturbed is calculated
        if((d_COccu[current] == 1 ) && (d_Rdom[current]<=(float)(Alpha0*Nu0+Lambda0)))
            Erosion=-2;      // Change induced by erosion (from 2 to 0)
        else
            Erosion=0;       // No erosion
        
        // Here, a possible transition from empty to occupied is calculated
        if( !d_COccu[current] && !d_CDist[current] && (d_Rdom[current]<=(float)(N*Alpha2)) )
			Colonization=+1; // Change induced by colonization (from 1 to 2)
        else
			Colonization=0;  // No colonization
        
        // Here, the transition from disturbed to empty is calculated
        if( d_CDist[current] == 1)
			BecomeEmpty=1;  // Disturbed cells (0) become bare (1)
        else
			BecomeEmpty=0;  // No change
        
        // Here, the Cell values are updated
        d_Cell[current]	=	d_Cell[current] + Colonization + Erosion + BecomeEmpty;
    }
}

