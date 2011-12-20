


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
static uint* dNeighbours = 0;
static uint* hNeighbours = 0;
#define BLOCK_SIZE 50

//==============================================
//== Declarations
texture<float, cudaTextureType1D, cudaReadModeElementType> texRef;

__global__ void calculateDistances(float* vector, uint numEntries, uint blockSize, uint gridSize);
__device__ float calculateEntries(dataEntry* first, dataEntry* second);
__device__ uint vectorIdx(uint x, uint y);
__device__ float sqr(float a);
__global__ void findNeighbours( uint numEntries, uint * output );


\
static const char* kIrisDistancesPath = "./data/iris_distances.data";

//==============================================
//== Functions
ErrorCode startCalculatingDistances() {

	ErrorCode err = LoadData();
	if (err != errOk) {
		return err;
	}

	dataStore* data = GetCurrDataStore();
	if (data == 0) {
		return GetLastErrorCode();
	}
//--------------------------------------

	// Create chanel descriptor for texture bindings
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );

	float* dData;
	// Allocate and copy data store to device memory
	uint size = data->info.numEntries * 4 * sizeof(float);
	cudaMalloc( &dData, size );
	cudaMemcpy( dData, data->dataVector, size, cudaMemcpyHostToDevice );

	// Set texture parameters
	texRef.addressMode[ 0] = cudaAddressModeWrap;
	texRef.addressMode[ 1] = cudaAddressModeWrap;
	texRef.filterMode = cudaFilterModeLinear;
	texRef.normalized = true;

	// Bind the array to the texture reference
	//cudaBindTextureToArray( texRef, cuArray );
	uint offset = 0;
	cudaBindTexture( &offset, &texRef, dData, &channelDesc, size );

	// Allocate result of transformation in device memory
	uint outputSize = data->info.numEntries * ( data->info.numEntries - 1 ) / 2;
	cudaMalloc( &dDistancesVector, outputSize * sizeof(float) );

	uint hGridSize = data->info.numEntries / BLOCK_SIZE;
	while (( hGridSize * BLOCK_SIZE ) < data->info.numEntries ) {
		hGridSize++;
	}

	dim3 dimBlock( BLOCK_SIZE, BLOCK_SIZE ); // thread per block
	dim3 dimGrid( hGridSize, hGridSize ); // blocks per grid

	calculateDistances<<<dimGrid, dimBlock>>>( dDistancesVector, data->info.numEntries, BLOCK_SIZE, hGridSize );

	cutilDeviceSynchronize();

	hDistancesVector = (float*)malloc( outputSize * sizeof(float) );

	if ( hDistancesVector == 0 ) {
		SetError( errNoMemory );

//		cudaFreeArray( cuArray );
		cudaFree( dData );
		cudaFree( dDistancesVector );
		return errNoMemory;
	}

	cudaMemcpy( hDistancesVector, dDistancesVector, outputSize * sizeof(float), cudaMemcpyDeviceToHost );	

	// no need for raw data - free it
	cudaFree( dData );

	// now bind distances to texture, so we could use it for neighbours
	cudaBindTexture( &offset, &texRef, dDistancesVector, &channelDesc, outputSize * sizeof(float) );	

	cudaMalloc( &dNeighbours, data->info.numEntries * MAX_NEIGHBOURS * sizeof(uint) );
	findNeighbours<<< BLOCK_SIZE, hGridSize >>>( data->info.numEntries, dNeighbours );

	cutilDeviceSynchronize();

	cudaFree( dDistancesVector );
	hNeighbours = (uint*)malloc( data->info.numEntries * MAX_NEIGHBOURS * sizeof(uint) );

	if ( hNeighbours == 0 ) {
		SetError( errNoMemory );

		cudaFree( dNeighbours );
		return errNoMemory;
	}

	cudaMemcpy( hNeighbours, dNeighbours, data->info.numEntries * MAX_NEIGHBOURS * sizeof(uint), cudaMemcpyDeviceToHost );

	//Save results to file
	saveDistanceData();

	cudaFree( dNeighbours );

	releaseDataStore();

	return err;
}
//==============================================
__global__ void calculateDistances( float* vector, uint numEntries, uint blockSize, uint gridSize ) {
	// global position of "first" thread
	uint firstCol = blockIdx.x * BLOCK_SIZE;
	uint col = firstCol + threadIdx.x;
	uint firstRow = blockIdx.y * BLOCK_SIZE;
	uint row = firstRow + threadIdx.y;
	bool boundryBlock = false;

	__shared__ dataEntry rowData[ BLOCK_SIZE];
	__shared__ dataEntry colData[ BLOCK_SIZE];

	// Check if this isn't external block
	if ( firstRow >= firstCol && col < numEntries && row < numEntries ) {
		// Check if we should care for loading colums here
		if ( threadIdx.y == 0 ) {
			// load columns
			colData[ threadIdx.x].a = tex1Dfetch( texRef, col*4+0 );
			colData[ threadIdx.x].b = tex1Dfetch( texRef, col*4+1 );
			colData[ threadIdx.x].c = tex1Dfetch( texRef, col*4+2 );
			colData[ threadIdx.x].d = tex1Dfetch( texRef, col*4+3 );
		}
		if ( row == col ) {
			boundryBlock = true;
			// don't load rows here
		} 

		// Check if we should care for loading rows
		if ( threadIdx.x == 0 && !boundryBlock ) {
			// load rows as wel
			rowData[ threadIdx.y].a = tex1Dfetch( texRef, row*4+0 );
			rowData[ threadIdx.y].b = tex1Dfetch( texRef, row*4+1 );
			rowData[ threadIdx.y].c = tex1Dfetch( texRef, row*4+2 );
			rowData[ threadIdx.y].d = tex1Dfetch( texRef, row*4+3 );
		}

		// Sync up threads
		__syncthreads();
		// And do some calculations

		// something to do for us here ?
		if ( row > col /*&& row == 1 && col == 0*/) {
			float distance = 0;
			if (boundryBlock) {
				distance = calculateEntries(&colData[threadIdx.x], &colData[threadIdx.y]);
			} else {
				distance = calculateEntries(&colData[threadIdx.x], &rowData[threadIdx.y]);
			}
			vector[vectorIdx(col, row)] = distance;
		} else {
			vector[vectorIdx( row, col )] = 2.0f;
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
	if ( y > x ) {
		return y * (y - 1) / 2 + x;
	} else {		
		return x * (x - 1) / 2 + y;
	}
	
}
//==============================================

__device__ float sqr(float a) {
	return a * a;
}
//==============================================

const float* getDistances() {
	if ( hDistancesVector == 0 ) {
		loadDistanceData();
	}
	FILE * dump;
	dump = fopen( "distances.txt", "w" );
	if ( dump ) {
		for ( int i = 0; i < numEntries(); i++ ) {
			for ( int j = 0; j <= i; j++ ) {
				if ( j == i ) {
					fprintf( dump, "\n");
				} else {
					fprintf( dump, " %f,", hDistancesVector[ i * (i - 1) / 2 + j] );
				}
			}
			
		}
		fclose( dump );
	}
	return hDistancesVector;
}
//==============================================

ErrorCode releaseDistances() {
	if ( hDistancesVector ) {
		free ( hDistancesVector );
	}

	return errOk;
}
//==============================================

__global__ void findNeighbours( uint numEntries, uint * output ) {
	uint neighbours[ MAX_NEIGHBOURS];
	float neighboursDistances[ MAX_NEIGHBOURS];
	float distance = 0;;
	int i;

	uint record = blockIdx.x * BLOCK_SIZE + threadIdx.x;

	if ( record >= numEntries ) {
		return;
	}

	for ( i = 0; i < MAX_NEIGHBOURS; i++ ) {
		neighbours[ i] = 0;
		neighboursDistances[ i] = 0;
	}

	// for each record in the data set
	for ( i = 0; i < numEntries; i++ ) {
		// if it's not the same		
		if ( record != i ) {
			// fetch distance
			distance = tex1Dfetch( texRef, vectorIdx( record, i ));			
			
			uint a = i;
			// for each neighbour already stored
			for ( int j = 0; j < MAX_NEIGHBOURS; j++ ) {
				// did we found proper one ?
				if ( neighboursDistances[ j] == 0 || distance < neighboursDistances[ j]) {					
					if ( neighboursDistances[ j] == 0 ) {
						// found empty entry - just save it here and break
						neighbours[ j] = a;
						neighboursDistances[ j] = distance;
						break;
					} else {
						// replace it and continue search with the one thrown out
						float d = neighboursDistances[ j];
						uint r = neighbours[ j];
						neighboursDistances[ j] = distance;
						neighbours[ j] = a;
						distance = d;
						a = r;
					}
				}
			} // for each neighbour
		}
	} // for each data entry

	// save results
	for ( i =0; i < MAX_NEIGHBOURS; i++ ) {
		output[ record * MAX_NEIGHBOURS + i] = neighbours[ i];
	}
}
//==============================================

const unsigned int* getNeighbours() {
	if ( hNeighbours == 0 ) {
		loadDistanceData();
	}
	return hNeighbours;
}
//==============================================

ErrorCode releaseNeighbours() {
	if ( hNeighbours != 0 ) {
		free( hNeighbours );
	}
	return errOk;
}
//==============================================

ErrorCode loadDistanceData() {
	unsigned int entries = 0;
	unsigned int inputSize = 0;

	ErrorCode ret = errOk;

	FILE * file = fopen( kIrisDistancesPath, "r" );
	size_t res = 0;
	if ( file ) {
		// read numOfEntries
		res = fread( &entries, sizeof(unsigned int), 1, file );
		if ( res == 1 ) {
			inputSize = entries * ( entries - 1 ) / 2;
			setNumEntries( entries );
			if ( hDistancesVector == 0 ) {
				hDistancesVector = (float*)malloc( inputSize * sizeof(float) );
			}
			// read distances
			res = fread( hDistancesVector, sizeof(float), inputSize, file );
			if ( res != inputSize ) {
				ret = errFileCorupted;
			} else {

				if ( hNeighbours == 0 ) {
					hNeighbours = (unsigned int*)malloc( entries * MAX_NEIGHBOURS * sizeof(unsigned int) );
				}
				// read neighbours
				res = fread( hNeighbours, MAX_NEIGHBOURS * sizeof(unsigned int), entries, file );
				if ( res != entries ) {
					res = errFileCorupted;
				}
			}
		}
		fclose( file );
		if ( ret != errOk ) {
			if ( hNeighbours != 0 ) free( hNeighbours );
			if ( hDistancesVector != 0 ) free( hDistancesVector );
			hNeighbours = 0;
			hDistancesVector = 0;
			setNumEntries( 0 );
		}
	}

	return ret;
}
//==============================================

ErrorCode saveDistanceData() {
	// check if we have something worh to save
	if ( hDistancesVector == 0 ) {
		return errNoData;
	}

	unsigned int entries = numEntries();
	unsigned int outputSize = entries * ( entries - 1 ) / 2;

	FILE * file = fopen( kIrisDistancesPath, "w" );
	size_t res = 0;
	ErrorCode err = errOk;
	if ( file ) {
		res = fwrite( &( entries ), sizeof(unsigned int), 1, file );
		if (res == 1) {
			res = fwrite( hDistancesVector, sizeof(float), outputSize, file );
		}

		if ( res == outputSize ) {
			if ( hNeighbours != 0 ) {
				res = 0;
				// read neighbours
				res = fwrite( hNeighbours, (MAX_NEIGHBOURS) * sizeof(unsigned int), entries, file );
			}
		}

		if ( res!=1 && res!=outputSize && res!=entries ) {
			err = errFileWrite;
		}

		fclose( file );
	} else {
		err = errFileWrite;
	}

	return err;
}
//==============================================