// Module responsible for loading data

//==============================================
//== Includes
#include "globals.cuh"
#include "errors.cuh"
#include <math.h>
#include "dataLoader.cuh"
#include "distanceCalculator.cuh"
#include <cutil_inline.h>
#include <cutil_math.h>
#include <cuda.h>

//==============================================
//== Globals

// num of entries per block
#define kBlockSize 16
// maximum dimensiosns ( used for some shared arrays )
#define kDimSize 14

static float*	dDistancesVector = 0;
static uint*	dNeighbours = 0;

//==============================================
//== Declarations
texture<float, cudaTextureType1D, cudaReadModeElementType> texRef;

__global__ void CalculateDistancesKernel( float* vector, uint numEntries, uint blockSize, uint gridSize, uint dimSize );
__global__ void CalculateNeighboursKernel( uint numEntries, uint * neighbours );
/**
 * Calculates distance betwen two data entries.
 * @param x			- [in] first entry
 * @param y			- [in] second entry
 * @param vector	- [in] vector with entries
 *
 * @return	calculated distance
 */
__device__ float CalculateEntries( float * x, float * y, uint dimSize );

__device__ unsigned int DistanceVIdx( unsigned int a, unsigned int b );

__device__ float sqr(float a);

//==============================================
//== Functions

/**
 * Calculates distances on given data.
 * @param dataSore - [in,out] 
 */
ErrorCode CalculateDistances( DataStore * dataStore );

/*
 * Checks for the distances already being saved.
 * @param	- [in] index of the loader to be used for data.
 *
 * @return	- true if file with saved distances for that loader exists.
 */
char CheckHaveSavedDistances( DataStore * dataStore );

/**
 * Saves Calculated results into file.
 * @param dataStore	- [in] pointer to structure with data to be saved.
 *
 * @return ErrorCode if any
 */
ErrorCode SaveCalculatedDistances( DataStore * dataStore );

/**
 * Loads data from file.
 * @param dataStore - [in,out] pointer to structure with info about file to be oppened.
 *
 * @return error code if any
 */
ErrorCode LoadCalculatedDistances( DataStore * dataStore );

/**
 * Calculates list of closests neighbours for each entry.
 * @param dataStore	- [in,out] holds required data and used to store results
 *
 * @return Error code if any.
 */
ErrorCode CalculateNeighbours( DataStore * dataStore );

//============================================================================


ErrorCode GetCalculatedDistances( unsigned int num, DataStore * dataStore ) {
	ErrorCode err;

	err = GetDataFromLoader( num, dataStore );

	if ( err != errOk ) {
		return err;
	}

	if ( CheckHaveSavedDistances( dataStore ) ) {
		err =  LoadCalculatedDistances( dataStore );
	} else {
		// We need to calculate them
		err = CalculateDistances( dataStore );

		if ( err == errOk ) {
			err = SaveCalculatedDistances( dataStore );
		}
	}

	return err;
}
//----------------------------------------------
//<<<( hGridSize, hGridSize ), ( BLOCK_SIZE, BLOCK_SIZE )>>>
__global__ void CalculateDistancesKernel( float* vector, uint numEntries, uint blockSize, uint gridSize, uint dimSize ) {
	// global position of "first" thread
	unsigned int firstCol = blockIdx.x * kBlockSize;
	unsigned int col = firstCol + threadIdx.x;
	unsigned int firstRow = blockIdx.y * kBlockSize;
	unsigned int row = firstRow + threadIdx.y;
	bool boundryBlock = false;

	__shared__ float rowData[ kBlockSize * kDimSize];
	__shared__ float colData[ kBlockSize * kDimSize];

	//--
	if ( col >= row || col >= numEntries || row >= numEntries ) {
		// External block, no calculations here
		return;
	}

	// Check if we should care for loading colums here
	if ( threadIdx.y == 0 ) {
		// load columns
		for ( int i = 0; i < dimSize; i++ ) {
			colData[ threadIdx.x * kDimSize + i] = tex1Dfetch( texRef, col * kDimSize + i );
		}
	}
	if ( row == col ) {
		boundryBlock = true;
		// don't load rows here
	} 

	// Check if we should care for loading rows
	if ( threadIdx.x == 0 && !boundryBlock ) {
		// load rows as wel
		for ( int i = 0; i < dimSize; i++ ) {
			rowData[ threadIdx.x * kDimSize + i] = tex1Dfetch( texRef, row * kDimSize + i );
		}
	}

	// Sync up threads
	__syncthreads();


	float distance = 0;
	if ( boundryBlock ) {
		distance = CalculateEntries( &colData[ threadIdx.x * kDimSize], &colData[ threadIdx.y * kDimSize], dimSize );
	} else {
		distance = CalculateEntries( &colData[ threadIdx.x * kDimSize], &rowData[ threadIdx.y * kDimSize], dimSize );
	}
	vector[ DistanceVIdx( col, row )] = distance;
}
//----------------------------------------------

char CheckHaveSavedDistances( DataStore * dataStore ) {
	unsigned int nameLen = 0;
	char * fileName = NULL;
	FILE * file = NULL;

	nameLen = strlen( dataStore->info.name );
	nameLen += strlen( "_distances.data" );
	nameLen += 1; // for null
	fileName = (char*)malloc( nameLen * sizeof(char) );
	sprintf( fileName, "%s_distances.data", dataStore->info.name );

	if ( fileName == NULL ) {
		reportError( errFailProcessData, "Failed to generate fileName. Got NULL%s", "" );
		return errFailProcessData;
	}

	file = fopen( fileName, "rb" );
	if ( file == NULL ) {
		return 0;
	}

	fclose( file );
	return 1;
}
//----------------------------------------------

ErrorCode CalculateDistances( DataStore * dataStore ) {
	ErrorCode err = errOk;

	//Check Params:
	//1) Wrong pointer
	if ( dataStore == NULL ) {
		reportError( errBadParam, "NULL given%s", "" );
		return SetLastErrorCode( errBadParam );
	}
	//2) got no data
	if ( dataStore->dataVector ==  NULL || dataStore->info.numEntries == 0 ) {
		reportError( errNoData, "vector=%x entries=%u", (unsigned int)dataStore->dataVector, dataStore->info.numEntries );
		return SetLastErrorCode( errNoData );
	}
	//=-

	// Create chanel descriptor for texture bindings
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );

	// Load data into device
	float* dData;
	uint dataSize = dataStore->info.numEntries * dataStore->info.dimSize * sizeof(float);
	cudaMalloc( &dData, dataSize );
	cudaMemcpy( dData, dataStore->dataVector, dataSize, cudaMemcpyHostToDevice );

	// Set texture parameters
	texRef.addressMode[ 0] = cudaAddressModeWrap;
	texRef.addressMode[ 1] = cudaAddressModeWrap;
	texRef.filterMode = cudaFilterModeLinear;
	texRef.normalized = true;

	// Bind the array to the texture reference
	uint offset = 0;
	cudaBindTexture( &offset, &texRef, dData, &channelDesc, dataSize );

	// Allocate result of transformation in device memory
	dataStore->info.distancesSize = dataStore->info.numEntries * ( dataStore->info.numEntries - 1 ) / 2;
	cudaMalloc( &dDistancesVector, dataStore->info.distancesSize * sizeof(float) );

	uint hGridSize = dataStore->info.numEntries / kBlockSize;
	while (( hGridSize * kBlockSize ) < dataStore->info.numEntries ) {
		hGridSize++;
	}

	dim3 dimBlock( kBlockSize, kBlockSize ); // thread per block
	dim3 dimGrid( hGridSize, hGridSize ); // blocks per grid

	CalculateDistancesKernel<<<dimGrid, dimBlock>>>( dDistancesVector, dataStore->info.numEntries, kBlockSize, hGridSize, dataStore->info.dimSize );

	cutilDeviceSynchronize();

	dataStore->distances = (float*)malloc( dataStore->info.distancesSize * sizeof(float) );

	cudaMemcpy( dataStore->distances, dDistancesVector, dataStore->info.distancesSize * sizeof(float), cudaMemcpyDeviceToHost );	
	// no need for raw data - free it
	cudaFree( dData );

	// now bind distances to texture, so we could use it for neighbours
	cudaBindTexture( &offset, &texRef, dDistancesVector, &channelDesc, dataStore->info.distancesSize * sizeof(float) );

	err = CalculateNeighbours( dataStore );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}

	return err;
}
//----------------------------------------------

__device__ unsigned int DistanceVIdx( unsigned int a, unsigned int b ) {
	if ( a > b ) {
		return a * (a - 1) / 2 + b;
	} else {		
		return b * (b - 1) / 2 + a;
	}
}
//----------------------------------------------

__device__ float CalculateEntries( float * x, float * y, uint dimSize ) {
	float result = 0;
	int i = 0;
	for (i = 0; i < dimSize; i++ ) {
		result += sqr( x[ i] - y[ i] );
	}
	result = sqrt( result );
	return result;
}
//----------------------------------------------

__device__ float sqr(float a) {
	return a * a;
}
//----------------------------------------------

ErrorCode LoadCalculatedDistances( DataStore * dataStore ) {
	// NAME_distances.data
	unsigned long nameLen = 0;
	char * fileName = NULL;
	FILE * file = NULL;
	unsigned long read = 0;

	if ( dataStore == NULL ||
		dataStore->info.name == NULL ) {
			reportError( errWrongParameter, "Got wrong parameters dataStore:%x, name:%x", (unsigned int)dataStore, (unsigned int)dataStore->info.name );
			return errWrongParameter;
	}

	nameLen = strlen( dataStore->info.name );
	nameLen += strlen( "_distances.data" );
	nameLen += 1; // for null
	fileName = (char*)malloc( nameLen * sizeof(char) );
	sprintf( fileName, "%s_distances.data", dataStore->info.name );

	if ( fileName == NULL ) {
		reportError( errFailProcessData, "Failed to generate fileName. Got NULL%s", "" );
		return errFailProcessData;
	}

	file = fopen( fileName, "rb" );
	if ( file == NULL ) {
		reportError( errFileRead, "Failed to open file(%s) for reading.", fileName );
		return errFileRead;
	}

	dataStore->info.distancesSize = dataStore->info.numEntries * ( dataStore->info.numEntries - 1 ) / 2;
	dataStore->distances = (float*)malloc( dataStore->info.distancesSize * sizeof(float) );
	checkAlloc( dataStore->distances )
		return GetLastErrorCode();
	}
	read = fread( dataStore->distances, sizeof(float), dataStore->info.distancesSize, file );
	if ( read != dataStore->info.distancesSize ) {
		reportError( errFileRead, "Failed to read distances data properly.%s", "" );
		return errFileRead;
	}

	dataStore->neighbours = (unsigned int*)malloc( dataStore->info.numEntries * kMaxNeighbours * sizeof(unsigned int) );
	checkAlloc( dataStore->neighbours )
		return GetLastErrorCode();
	}
	read = fread( dataStore->neighbours, kMaxNeighbours * sizeof(unsigned int), dataStore->info.numEntries, file );
	if ( read != dataStore->info.numEntries ) {
		reportError( errFileRead, "Failed to read neighbours data properly.%s", "" );
		return errFileRead;
	}

	fclose( file );

	return errOk;
}
//----------------------------------------------

ErrorCode SaveCalculatedDistances( DataStore * dataStore ) {
	// NAME_distances.data
	unsigned long nameLen = 0;
	char * fileName = NULL;
	FILE * file = NULL;
	unsigned long written = 0;

	if ( dataStore == NULL ||
		dataStore->info.name == NULL ||
		dataStore->distances == NULL ||
		dataStore->neighbours == NULL ) {
			reportError( errWrongParameter, "Got wrong parameters dataStore:%x, name:%x", (unsigned int)dataStore, (unsigned int)dataStore->info.name );
			return errWrongParameter;
	}

	nameLen = strlen( dataStore->info.name );
	nameLen += strlen( "_distances.data" );
	nameLen += 1; // for null
	fileName = (char*)malloc( nameLen * sizeof(char) );
	sprintf( fileName, "%s_distances.data", dataStore->info.name );

	if ( fileName == NULL ) {
		reportError( errFailProcessData, "Failed to generate fileName. Got NULL%s", "" );
		return errFailProcessData;
	}

	file = fopen( fileName, "wb" );
	if ( file == NULL ) {
		reportError( errFileWrite, "Failed to open file(%s) for writting.", fileName );
		return errFileWrite;
	}

	// write distances
	written = fwrite( dataStore->distances, sizeof(float), dataStore->info.distancesSize, file );
	if ( written != dataStore->info.distancesSize ) {
		reportError( errFileWrite, "Failed to write distances data properly.%s","" );
		return errFileWrite;
	}
	// write neighbours
	written = fwrite( dataStore->neighbours, kMaxNeighbours * sizeof(unsigned int), dataStore->info.numEntries, file );
	if ( written != dataStore->info.numEntries ) {
		reportError( errFileWrite, "Failed to write neighbours data properly.%s", "" );
		return errFileWrite;
	}

	fclose( file );

	return errOk;
}
//----------------------------------------------

ErrorCode CalculateNeighbours( DataStore *dataStore ) {
	ErrorCode err = errOk;
	unsigned int gridSize = 0;

	//Check Params:
	//1) Wrong pointer
	if ( dataStore == NULL ) {
		reportError( errBadParam, "NULL given%s", "" );
		return SetLastErrorCode( errBadParam );
	}
	//2) got no data
	if ( dataStore->distances ==  NULL || dataStore->info.distancesSize == 0 ) {
		reportError( errNoData, "vector=%x entries=%u", (unsigned int)dataStore->dataVector, dataStore->info.numEntries );
		return SetLastErrorCode( errNoData );
	}
	//=-

	gridSize = dataStore->info.numEntries / kBlockSize;
	while (( gridSize * kBlockSize ) < dataStore->info.numEntries ) {
		gridSize++;
	}

	cudaMalloc( &dNeighbours, dataStore->info.numEntries * kMaxNeighbours * sizeof(uint) );
	CalculateNeighboursKernel<<< gridSize, kBlockSize >>>( dataStore->info.numEntries, dNeighbours );

	cutilDeviceSynchronize();

	cudaFree( dDistancesVector );
	dataStore->neighbours = (unsigned int*)malloc( dataStore->info.numEntries * kMaxNeighbours * sizeof( unsigned int ) );

	checkAlloc( dataStore->neighbours )
		cudaFree( dNeighbours );
		return GetLastErrorCode();
	}

	cudaMemcpy( dataStore->neighbours, dNeighbours, dataStore->info.numEntries * kMaxNeighbours * sizeof(uint), cudaMemcpyDeviceToHost );

	//Save results to file
	//saveDistanceData(); // TODO:

	cudaFree( dNeighbours );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}

    return err;
}
//----------------------------------------------

//<<< hGridSize, kBlockSize >>>
__global__ void CalculateNeighboursKernel( uint numEntries, uint * output ) {
	unsigned int record = 0;
	unsigned int i = 0;
	unsigned int j = 0;
	float distance = 0;
	uint neighbours[ kMaxNeighbours];
	float neighboursDistances[ kMaxNeighbours];
	unsigned int candidate = 0;
	unsigned int a = 0;
	float b = 0;

	record = blockIdx.x * kBlockSize + threadIdx.x;

	if ( record >= numEntries ) {
		return;
	}

	// fill with first records from the list as candidates
	for ( j = 0; j < kMaxNeighbours; j++ ) {
		neighboursDistances[ j] = 0;
	}
	
	// for each record in the data set
	for ( i = 0; i < numEntries; i++ ) {
		if ( record == i ) {
			// skip this one
			continue;
		}

		distance = tex1Dfetch( texRef, DistanceVIdx( record, i ));	
		candidate = i;

		for ( j = 0; j < kMaxNeighbours; j++ ) {
			if ( neighbours[ record * kMaxNeighbours + j] == i ) {
				// it's already here - break;
				break;
			}

			if ( neighboursDistances[ j] == 0 ) {
				neighbours[ record * kMaxNeighbours + j] = candidate;
				neighboursDistances[ j] = distance;
				break;
			}

			if ( distance < neighboursDistances[ j] ) {
				// continue with a little buble sorting
				a = neighbours[ record * kMaxNeighbours + j];
				b = neighboursDistances[ j];
				neighbours[ record * kMaxNeighbours + j] = candidate;
				neighboursDistances[ j] = distance;
				candidate = a;
				distance = b;
			}
		} // for each neighbour
	} // for each data entry

	// save results
	for ( i = 0; i < kMaxNeighbours; i++ ) {
		output[ record * kMaxNeighbours + i] = neighbours[ i];
	}
}