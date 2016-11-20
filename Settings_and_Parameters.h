//
//  Settings_and_Parameter.cpp
//
//  Created by Johan Van de Koppel on 03-09-14.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

// The code now works only on an NVidia card, so choose that card here
// To find you computer's setup, switch Print_All_Devices to ON in
// Device_Utilities.h
#define Device_No         1    // 0: CPU; 1: Intel 4000; 2: Nvidia
#define ProgressBarWidth 45

// Thread block size
#define Block_Size_X	128               // 16
#define Block_Size_Y	128               // 16

// Number of Blocks
/* I define the Block_Number_ensions of the matrix as product of two numbers
Makes it easier to keep them a multiple of something (16, 32). For the Mersenne-
Twister code, it is essential to keep to GridWith and GridHeight a power of 2,
e.g. 512, 1024, 2048, etc */

#define Block_Number_X	16               // 32
#define Block_Number_Y	16               // 32

// Matrix Block_Number_ensions
// (chosen as multiples of the thread block size for simplicity)
#define Grid_Width  (Block_Size_X * Block_Number_X)			// Matrix A width
#define Grid_Height (Block_Size_Y * Block_Number_Y)			// Matrix A height
#define Grid_Size (Grid_Width*Grid_Height)                  // Grid Size

// Label definitions
#define HORIZONTAL	201
#define VERTICAL	202
#define CELLS		101
#define RAND		102
#define	DIS         0
#define INB         1
#define OCC         2

// Initial conditions
#define Frac		0.1     // The fraction of the domain this is covered with mussels
#define DFrac		0.05    // The fraction of the domain that is disturbed

// Model parameters
#define Alpha0		0.9	    // Disturbance probability if one of neighbours is disturbed
#define Alpha2		0.6	    // Recovery probability coefficient
#define Lambda0		0.002   // Probability that an occupied cell becoming disturbed by chance

// Simulation settings
#define EndTime		1200	// Number of timesteps in the whole simulation
#define Time		0
#define StoreNow	1
#define NumFrames	25
#define MAX_STORE	(NumFrames+1)

#define dT          1










