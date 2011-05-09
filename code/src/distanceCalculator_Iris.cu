


// Module responsible for loading data

//==============================================
//== Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "globals.cuh"
#include "errors.cuh"
#include "distanceCalculator_Iris.cuh"
#include "dataLoader.cuh"
#include <cutil_inline.h>
#include <cutil_math.h>
#include <cuda.h>
//#include <shrQATest.h>

//==============================================
//== Types

//==============================================
//== Constants and Globals
static float* dDistancesVector = 0;
static float* hDistancesVector = 0;

//==============================================
//== Declarations
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;

__global__ void calculateDistances(float* vector, uint numEntries, uint blockSize, uint gridSize);
__device__ float calculateEntries(dataEntry* first, dataEntry* second);
__device__ uint vectorIdx(uint x, uint y);
__device__ float sqr(float a);

#define BLOCK_SIZE 16
static const char* kIrisDistancesPath = "./data/iris_distances.data";

//==============================================
//== Functions
ErrorCode StartCalculatingDistances() {

	ErrorCode err = LoadData();
	if (err != errOk) {
		return err;
	}

	dataStore* data = GetCurrDataStore();
	if (data == 0) {
		return GetLastErrorCode();
	}
//--------------------------------------

	// Allocate CUDA array in device memory
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );

	cudaArray* cuArray;
	cudaMallocArray( &cuArray, &channelDesc, data->info.numEntries, 4 );
	
	// Copy to device memory some data located at address h_data 
	// in host memory 
	uint size = data->info.numEntries * sizeof(dataEntry);
	cudaMemcpyToArray( cuArray, 0, 0, data->dataVector, size, cudaMemcpyHostToDevice );

	// Set texture parameters
	texRef.addressMode[ 0] = cudaAddressModeClamp;
	texRef.addressMode[ 1] = cudaAddressModeClamp;
	texRef.filterMode = cudaFilterModeLinear;
	texRef.normalized = true;

	// Bind the array to the texture reference
	cudaBindTextureToArray(texRef, cuArray);

	// Allocate result of transformation in device memory
	uint outputSize = data->info.numEntries * (data->info.numEntries - 1) / 2;
	cudaMalloc( &dDistancesVector, outputSize * sizeof(float) );

	uint hGridSize = data->info.numEntries / BLOCK_SIZE;
	while (( hGridSize * BLOCK_SIZE ) < data->info.numEntries ) {
		hGridSize++;
	}

	dim3 dimBlock( BLOCK_SIZE, BLOCK_SIZE ); // thread per block
	dim3 dimGrid( hGridSize, hGridSize ); // blocks per grid

	calculateDistances<<<dimGrid, dimBlock>>>( dDistancesVector, data->info.numEntries, BLOCK_SIZE, hGridSize );

	hDistancesVector = (float*)malloc( outputSize * sizeof(float) );

	if ( hDistancesVector == 0 ) {
		SetError( errNoMemory );

		cudaFreeArray( cuArray );
		cudaFree( dDistancesVector );
		return errNoMemory;
	}

	cudaMemcpy( hDistancesVector, dDistancesVector, outputSize * sizeof(float), cudaMemcpyDeviceToHost );

	//Save results to file
	FILE * file = fopen( kIrisDistancesPath, "w" );
	size_t res = 0;
	if ( file ) {
		res = fwrite( &( data->info.numEntries ), sizeof(unsigned int), 1, file );
		if (res == 1) {
			res = fwrite( hDistancesVector, sizeof(float), outputSize, file );
		}
		if ( res!=1 && res!=outputSize ) {
			err = errFileWrite;
		}

		fclose( file );
	} else {
		err = errFileWrite;
	}	

	// Free device memory
	cudaFreeArray( cuArray );
	cudaFree( dDistancesVector );	

	return err;
}
//==============================================

__global__ void calculateDistances( float* vector, uint numEntries, uint blockSize, uint gridSize ) {
	// Check if thread of external block

	// global position of "first" thread
	uint firstCol = blockIdx.x * BLOCK_SIZE;
	uint col = firstCol + threadIdx.x;
	uint firstRow = blockIdx.y * BLOCK_SIZE;
	uint row = firstRow + threadIdx.y;
	bool boundryBlock = false;

	__shared__ dataEntry rowData[ BLOCK_SIZE];
	__shared__ dataEntry colData[ BLOCK_SIZE];

	// Check if this isn't external block
	if ( row >= col && col < numEntries && row < numEntries ) {
		// Check if we should care for loading colums here
		if ( threadIdx.y == 0 ) {
			// load columns
			float u = (float)col / (float)numEntries;			
			colData[ threadIdx.x].a = tex2D( texRef, u, 0.0f / (float)numEntries );
			colData[ threadIdx.x].a = tex2D( texRef, u, 1.0f / (float)numEntries );
			colData[ threadIdx.x].a = tex2D( texRef, u, 2.0f / (float)numEntries );
			colData[ threadIdx.x].a = tex2D( texRef, u, 3.0f / (float)numEntries );
		}
		if ( row = col ) {
			boundryBlock = true;
			// don't load rows here
		} 

		// Check if we should care for loading rows
		if ( threadIdx.x == 0 && !boundryBlock ) {
			// load rows as wel
			float u = (float)row / (float)numEntries;
			rowData[ threadIdx.y].a = tex2D( texRef, u, 0.0f / (float)numEntries );
			rowData[ threadIdx.y].a = tex2D( texRef, u, 1.0f / (float)numEntries );
			rowData[ threadIdx.y].a = tex2D( texRef, u, 2.0f / (float)numEntries );
			rowData[ threadIdx.y].a = tex2D( texRef, u, 3.0f / (float)numEntries );
		}

		// Sync up threads
		__syncthreads();
		// And do some calculations

		// something to do for us here ?
		if ( row > col ) {
			float distance = 0;
			if (boundryBlock) {
				distance = calculateEntries(&colData[threadIdx.x], &colData[threadIdx.y]);
			} else {
				distance = calculateEntries(&colData[threadIdx.x], &rowData[threadIdx.y]);
			}
			vector[vectorIdx(col, row)] = distance;
		}

	}
}
//==============================================

__device__ float calculateEntries(dataEntry* first, dataEntry* second) {
	float result = 0;

	result = sqr(first->a - second->a);
	result += sqr(first->b - second->b);
	result += sqr(first->c - second->c);
	result += sqr(first->d - second->d);
	result = sqrt(result);

	return result;
}
//==============================================

__device__ uint vectorIdx(uint x, uint y) {
	// Assuming that always y > x
	// Is this safe ?
	return y * (y - 1) / 2;
}
//==============================================

__device__ float sqr(float a) {
	return a * a;
}
//==============================================

float* GetDistances() {
	return hDistancesVector;
}
//==============================================

ErrorCode ReleaseDistances() {
	if ( hDistancesVector ) {
		free ( hDistancesVector );
	}

	return errOk;
}
//==============================================
