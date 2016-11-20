//
//  Main_program.cpp
//
//  Created by Johan Van de Koppel on 03-09-14.
//  Copyright (c) 2014 Johan Van de Koppel. All rights reserved.
//

#include <sys/time.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "Settings_and_Parameters.h"
#include "clMersenneTwister.h"
#include "clMersenneTwister_functions.h"
#include "Device_Utilities.h"
#include "Initial_values.h"

#define MAX_SOURCE_SIZE (0x100000)

// Forward definitions from functions at the end of this code file
void initialize_RandomCellMatrix(float*,int,int);
void initialize_Cells(float*,float*,int,int);

////////////////////////////////////////////////////////////////////////////////
// Main program code for Aridlands
////////////////////////////////////////////////////////////////////////////////

int main()
{
    
    /*----------Constant and variable definition------------------------------*/
    
    // The location of the kernel code is obtain from the __FILE__ macro
    const std::string SourcePath (__FILE__);
    const std::string PathName = SourcePath.substr (0,SourcePath.find_last_of("/")+1);
    
	unsigned int Grid_Memory = sizeof(float) * Grid_Size;
	unsigned int size_storegrid = Grid_Width * Grid_Height * MAX_STORE;
	unsigned int mem_size_storegrid = sizeof(float) * size_storegrid;
    
    /*----------Mersenne Twister specific settings and variables--------------*/
    
    const int seed = 66e6;                          // Seeding number for Random Number generation (not used)
    const int nPerRng = Grid_Size/MT_RNG_COUNT;     // # of recurrence steps, must be even if do Box-Muller transformation
    const int nRand = MT_RNG_COUNT * nPerRng;       // Output size
    
    const std::string MTKernelFile = PathName + "clMersenneTwister.cl";
    const std::string MTDataName = PathName + "MT_data/MersenneTwister.dat";
    const std::string MTRawName = PathName + "MT_data/MersenneTwister.raw";
    
    /*----------Defining and allocating memeory on host-----------------------*/
    
    // Defining and allocating the memory blocks
	float* h_Cells      = (float*) malloc(Grid_Memory);
    float* h_Rand       = (float*) malloc(Grid_Memory);
	float* h_storegrid	= (float*) malloc(mem_size_storegrid);
    
    /*----------Initializing the host arrays----------------------------------*/
    
    srand(50);
	initialize_RandomCellMatrix(h_Rand,Grid_Height,Grid_Width);
    initialize_Cells(h_Cells, h_Rand, Grid_Height,Grid_Width);
    
    /*----------Printing info to the screen ----------------------------------*/

    Print_Label();
    
	printf(" Current grid dimensions: %d x %d cells\n\n",
           Grid_Width, Grid_Height);
    
    /*----------Setting up the device and the kernel--------------------------*/
    
    cl_device_id* devices;
    cl_int err;
    
    cl_context context = CreateGPUcontext(devices);
    
    // Print the name of the device that is used
    printf(" Implementing PDE on device %d: ", Device_No);
    print_device_info(devices, (int)Device_No);
    printf("\n");
    
    // Create a command queue on the device
    cl_command_queue command_queue = clCreateCommandQueue(context, devices[Device_No], 0, &err);
    
    /* Create Buffer Objects */
	cl_mem d_Cells = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
	cl_mem d_CellsDisturbd = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
    cl_mem d_CellsOccupied = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
	cl_mem d_Rand = clCreateBuffer(context, CL_MEM_READ_WRITE, Grid_Memory, NULL, &err);
    
    for (int i=0;i<Grid_Size;i++) {h_Rand[i]=0;};
    
	/* Copy input data to the memory buffer */
	err = clEnqueueWriteBuffer(command_queue, d_Cells, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_Cells, 0, NULL, NULL);
	err = clEnqueueWriteBuffer(command_queue, d_Rand, CL_TRUE, 0, Grid_Width*Grid_Height*sizeof(float), h_Rand, 0, NULL, NULL);
    
    /*----------Mussel disturbance Kernel-------------------------------------*/
    
    cl_program MDM_program = BuildKernelFile("MDM_Kernel.cl", context, &devices[Device_No], &err);
    if (err!=0)  printf(" > Compile Program Error number: %d \n\n", err);
    
	/* Link kernel */
	cl_kernel MDM_kernel = clCreateKernel(MDM_program, "MusselDisturbKernel", &err);
    if (err!=0) printf(" > Create Kernel Error number: %d \n\n", err);
   
    /* Set OpenCL kernel arguments */
	//MusselDisturbKernel (float* d_Cell, float* d_Rdom, float* d_CDist, float* d_COccu)
    err = clSetKernelArg(MDM_kernel, 0, sizeof(cl_mem), (void *)&d_Cells);
    err = clSetKernelArg(MDM_kernel, 1, sizeof(cl_mem), (void *)&d_Rand);
    err = clSetKernelArg(MDM_kernel, 2, sizeof(cl_mem), (void *)&d_CellsDisturbd);
	err = clSetKernelArg(MDM_kernel, 3, sizeof(cl_mem), (void *)&d_CellsOccupied);
    
    /* ---------Setting up the Mersenne Twister kernel------------------------*/
    
    mt_struct_stripped *h_MT = (mt_struct_stripped*)malloc(sizeof(mt_struct_stripped)*MT_RNG_COUNT);
    
    #ifdef DO_BOXMULLER
        cl_kernel ckBoxMuller = NULL;               // BoxMuller kernel
    #endif
    
    loadMTGPU(MTDataName.c_str(), seed, h_MT, MT_RNG_COUNT);
    initMTRef(MTRawName.c_str());
    
    // Setting up the MT algorithm data on the GPU
    cl_mem d_MT = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(mt_struct_stripped)*MT_RNG_COUNT, NULL, &err);
    err |= clEnqueueWriteBuffer(command_queue, d_MT, CL_TRUE, 0,
                                   sizeof(mt_struct_stripped)*MT_RNG_COUNT, h_MT, 0, NULL, NULL);
    
    /*----------Loading and compiling of the Mersenne Twister kernel code-----*/
    
    cl_program MT_program = BuildKernelFile("clMersenneTwister.cl", context,
                                 &devices[Device_No], &err);
    
    // Linking to the specific kernel function
    cl_kernel clMersenneTwister = clCreateKernel(MT_program, "MersenneTwister", &err);
    
    // Setting the parameters of the kernel function
    err |= clSetKernelArg(clMersenneTwister, 0, sizeof(cl_mem), (void*)&d_Rand);
    err |= clSetKernelArg(clMersenneTwister, 1, sizeof(cl_mem), (void*)&d_MT);
    err |= clSetKernelArg(clMersenneTwister, 2, sizeof(int),    (void*)&nPerRng);
    
    #ifdef DO_BOXMULLER
        // Linking to the specific kernel function
        ckBoxMuller = clCreateKernel(cpProgram, "BoxMuller", &err);
        
        // Setting the parameters of the kernel function
        err |= clSetKernelArg(ckBoxMuller, 0, sizeof(cl_mem), (void*)&d_Rand);
        err |= clSetKernelArg(ckBoxMuller, 1, sizeof(int),    (void*)&nPerRng);
    #endif
    
    // Setting up MT Kernel execution parameters
    size_t globalWorkSize= MT_RNG_COUNT;    // 1D var for Total # of work items
    size_t localWorkSize = nPerRng;         // 1D var for # of work items in the work group
    
    //RandomInfo(nPerRng, Grid_Height, Grid_Width);
    
    /*----------End of Mersenne Twister Kernel--------------------------------*/
    
    /*----------Set op timer and progress bar---------------------------------*/
    
    // create and start timer
    struct timeval Time_Measured;
    gettimeofday(&Time_Measured, NULL);
    double Time_Begin=Time_Measured.tv_sec+(Time_Measured.tv_usec/1000000.0);

    /* Progress bar initiation */
    int RealBarWidth=std::min((int)NumFrames,(int)ProgressBarWidth);
    int BarCounter=0;
    double BarThresholds[RealBarWidth];
    for (int i=1;i<RealBarWidth;i++)
        { BarThresholds[i] = (double)(i+1)/(double)RealBarWidth*(double)NumFrames; };
    
    /* Print the reference bar */
    printf(" Progress: [");
    for (int i=0;i<RealBarWidth;i++) { printf("-"); }
    printf("]\n");
    fprintf(stderr, "           >");
    
    /*----------Kernel parameterization---------------------------------------*/

    size_t global_item_size= Grid_Size;     // 1D var for total # of work items
    size_t local_item_size = Block_Number_X*Block_Number_Y; // 1D var for # of work items in the work group

    /*----------Calculation loop----------------------------------------------*/
	for (int Counter=0;Counter<NumFrames;Counter++)
    {
        for (int Runtime=0;Runtime<floor((float)EndTime/NumFrames/dT);Runtime++)
        {
            
            // Running the kernel for 1 time
            err |= clEnqueueNDRangeKernel(command_queue, clMersenneTwister, 1, NULL,
                                             &globalWorkSize, &localWorkSize, 0, NULL, NULL);
            if (err!=0) { printf(" > MT Kernel Error number: %d \n\n", err); exit(-10);}

            
            // Execute MDM kernel
            err = clEnqueueNDRangeKernel(command_queue, MDM_kernel, 1, NULL,
                                         &global_item_size, &local_item_size, 0, NULL, NULL);
            if (err!=0) { printf(" > MDM Kernel Error number: %d \n\n", err); exit(-10);}
            
        }
        
        // Transfer result to host
        
        err |= clEnqueueReadBuffer(command_queue, d_Cells, CL_TRUE, 0,
                                      sizeof(cl_float) * nRand, h_Cells, 0, NULL, NULL);
        if (err!=0) { printf("Read Buffer Error: %d\n\n", err); }
        
        //Store values at this frame.
        memcpy(h_storegrid+(Counter*Grid_Size),h_Cells,Grid_Memory);
        
        // Progress the progress bar if time
        if ((float)(Counter+1)>=BarThresholds[BarCounter]) {
            fprintf(stderr,"*");
            BarCounter = BarCounter+1;}
            
    }
    
    printf("<\n\n");
    
    /*---------------------Report on time spending----------------------------*/
    
    gettimeofday(&Time_Measured, NULL);
    double Time_End=Time_Measured.tv_sec+(Time_Measured.tv_usec/1000000.0);
	printf(" Processing time: %4.5f (s) \n\n", Time_End-Time_Begin);
    
    /*---------------------Write to file now----------------------------------*/
    const std::string DataPath = PathName + "MDM_Output.dat";
    FILE * fp=fopen(DataPath.c_str(),"wb");

    int width_matrix = Grid_Width;
    int height_matrix = Grid_Height;
    int NumStored = NumFrames;
    int EndTimeVal = EndTime;

	// Storing parameters
	fwrite(&width_matrix,sizeof(int),1,fp);
	fwrite(&height_matrix,sizeof(int),1,fp);
	fwrite(&NumStored,sizeof(int),1,fp);
	fwrite(&EndTimeVal,sizeof(int),1,fp);
	
	for(int store_i=0;store_i<NumFrames;store_i++)
    {
		fwrite(&h_storegrid[store_i*Grid_Size],sizeof(float),Grid_Size,fp);
    }
    
	fclose(fp);
    
	/*---------------------Clean up memory------------------------------------*/
	
    // Freeing host space
    free(h_Cells);
	free(h_Rand);
	free(h_storegrid);
 
	// Freeing kernel and block space
	err = clFlush(command_queue);
	err = clFinish(command_queue);
	err = clReleaseKernel(MDM_kernel);
	err = clReleaseProgram(MDM_program);
	err = clReleaseMemObject(d_Cells);
	err = clReleaseMemObject(d_CellsDisturbd);
    err = clReleaseMemObject(d_CellsOccupied);
    err = clReleaseMemObject(d_Rand);
	err = clReleaseCommandQueue(command_queue);
	err = clReleaseContext(context);
    
    // Kernel memory release
    clReleaseKernel(clMersenneTwister);
    
    #ifdef DO_BOXMULLER
        clReleaseKernel(ckBoxMuller);
    #endif
    
    clReleaseProgram(MT_program);
    clReleaseMemObject(d_MT);
    
    free(h_MT);
    
    free(devices);
    
    #if defined(__APPLE__) && defined(__MACH__)
        system("say Simulation finished");
    #endif

	return 0;
}


