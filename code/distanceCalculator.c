/**
 8/11/2011
 Jaroslaw Wojtasik

 noCuda

 distanceCalculator.c
 **/

#include "distanceCalculator.h"
#include "dataLoader.h"
#include "errors.h"
#include "loops.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
//==============================================
//== Globals

#define kBlockSize 16
#define kMaxNeighbours 10

//==============================================
//== Functions

/**
 * Kernel function to run inside the loop
 * @param loop - [in] structure that holds all context of the loop.
 */
void CalculateDistancesKernel( LoopContext loop );

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
char CheckHaveSavedDistances( unsigned int loaderIndex );

/**
 * Calculates distance betwen two data entries.
 * @param x			- [in] first entry
 * @param y			- [in] second entry
 * @param dataStore - [in] data store where entries are stored
 *
 * @return	calculated distance
 */
float CalculateEntries( unsigned int x, unsigned int y, DataStore * dataStore );

/**
 * Saves Calculated results into file.
 * @param dataStore	- [in] pointer to structure with data to be saved.
 *
 * @return ErrorCode if any
 */
ErrorCode SaveCalculatedDistances( DataStore * dataStore );

/**
 * Calculates list of closests neighbours for each entry.
 * @param dataStore	- [in,out] holds required data and used to store results
 *
 * @return Error code if any.
 */
ErrorCode CalculateNeighbours( DataStore * dataStore );

/**
 * Kernel function to run inside the loop
 * @param loop - [in] structure that holds all context of the loop.
 */
void CalculateNeighboursKernel( LoopContext loop );
//============================================================================


ErrorCode GetCalculatedDistances( unsigned int num, DataStore * dataStore ) {
	ErrorCode err;

	err = GetDataFromLoader( num, dataStore );

	if ( err != errOk ) {
		return err;
	}

	if ( CheckHaveSavedDistances( num ) ) {
		// TODO: load distances from the file
	} else {
		// We need to calculate them
		err = CalculateDistances( dataStore );

		if ( err == errOk ) {
			err = CalculateNeighbours( dataStore );
		}

		if ( err == errOk ) {
			err = SaveCalculatedDistances( dataStore );
		}
	}

	return err;
}
//----------------------------------------------

void CalculateDistancesKernel( LoopContext loop ) {
	// global position of "first" thread
	unsigned int firstCol = loop.blockIdx.x * kBlockSize;
	unsigned int col = firstCol + loop.threadIdx.x;
	unsigned int firstRow = loop.blockIdx.y * kBlockSize;
	unsigned int row = firstRow + loop.threadIdx.y;
	DistancesParams * params = (DistancesParams*)loop.params;

	//--
	if ( col >= row || col >= params->dataStore->info.numEntries || row >= params->dataStore->info.numEntries ) {
		// External block, no calculations here
		return;
	}

	params->dataStore->distances[ DistanceVIdx( col, row )] = CalculateEntries( col, row, params->dataStore );
}
//----------------------------------------------

char CheckHaveSavedDistances( unsigned int loaderIndex ) {
	return 0;
}
//----------------------------------------------

ErrorCode CalculateDistances( DataStore * dataStore ) {
	LoopDefinition distanceLoop;
	DistancesParams params;
	ErrorCode err = errOk;

	//Check Params:
	//1) Wrong pointer
	if ( dataStore == NULL ) {
		reportError( errBadParam, "NULL given" );
		return SetLastErrorCode( errBadParam );
	}
	//2) got no data
	if ( dataStore->dataVector ==  NULL || dataStore->info.numEntries == 0 ) {
		reportError( errNoData, "vector=%x entries=%u", dataStore->dataVector, dataStore->info.numEntries );
		return SetLastErrorCode( errNoData );
	}

	dataStore->info.distancesSize = dataStore->info.numEntries * ( dataStore->info.numEntries - 1 ) / 2;
	dataStore->distances = (float*)malloc( dataStore->info.distancesSize * sizeof(float) );

	checkAlloc( dataStore->distances )
		return GetLastErrorCode();
	}

	params.gridSize = dataStore->info.numEntries / kBlockSize;
	while (( params.gridSize * kBlockSize ) < dataStore->info.numEntries ) {
		params.gridSize++;
	}

	params.blockSize = kBlockSize;
	params.dataStore = dataStore;
	
	distanceLoop.blockSize.x = distanceLoop.blockSize.y = kBlockSize;
	distanceLoop.gridSize.x = distanceLoop.gridSize.y = params.gridSize;
	distanceLoop.blockSize.z = distanceLoop.gridSize.z = 1;
	distanceLoop.kernel = CalculateDistancesKernel;
	distanceLoop.params = (void*)&params;

	err = RunLoop( distanceLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error" );
	}

	return err;
}
//----------------------------------------------

unsigned int DistanceVIdx( unsigned int a, unsigned int b ) {
	if ( a > b ) {
		return a * (a - 1) / 2 + b;
	} else {		
		return b * (b - 1) / 2 + b;
	}
}
//----------------------------------------------

float CalculateEntries( unsigned int x, unsigned int y, DataStore * dataStore) {
	float result = 0;
	int i = 0;
	for (i = 0; i < dataStore->info.dimSize; i++ ) {
		result += pow( dataStore->dataVector[ x * dataStore->info.dimSize + i] -
			dataStore->dataVector[ y * dataStore->info.dimSize + i], 2.0 );
	}
	result = sqrt( result );
	return result;
}
//----------------------------------------------

ErrorCode SaveCalculatedDistances( DataStore * dataStore ) {
	// save the results
	return errOk;
}
//----------------------------------------------

ErrorCode CalculateNeighbours( DataStore *dataStore ) {
	LoopDefinition neighbourLoop;
	DistancesParams params;
	ErrorCode err = errOk;

	//Check Params:
	//1) Wrong pointer
	if ( dataStore == NULL ) {
		reportError( errBadParam, "NULL given" );
		return SetLastErrorCode( errBadParam );
	}
	//2) got no data
	if ( dataStore->distances ==  NULL || dataStore->info.distancesSize == 0 ) {
		reportError( errNoData, "vector=%x entries=%u", dataStore->dataVector, dataStore->info.numEntries );
		return SetLastErrorCode( errNoData );
	}

	dataStore->neighbours = (unsigned int*)malloc( dataStore->info.numEntries * kMaxNeighbours * sizeof( unsigned int ) );

	checkAlloc( dataStore->neighbours )
		return GetLastErrorCode();
	}

	params.gridSize = dataStore->info.numEntries / kBlockSize;
	while (( params.gridSize * kBlockSize ) < dataStore->info.numEntries ) {
		params.gridSize++;
	}

	params.blockSize = kBlockSize;
	params.dataStore = dataStore;
	
	neighbourLoop.blockSize.x = kBlockSize;
	neighbourLoop.gridSize.x = params.gridSize;
	neighbourLoop.blockSize.y = neighbourLoop.blockSize.z = neighbourLoop.gridSize.y = neighbourLoop.gridSize.z = 1;
	neighbourLoop.kernel = CalculateNeighboursKernel;
	neighbourLoop.params = (void*)&params;

	err = RunLoop( neighbourLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error" );
	}
}
//----------------------------------------------

void CalculateNeighboursKernel( LoopContext loop ) {
	unsigned int record = 0;
	int i = 0;
	int j = 0;
	float distance = 0;
	float distancesArray[ kMaxNeighbours];
	char shift = 0;
	unsigned int candidate = 0;
	unsigned int a = 0;
	float b = 0;

	DistancesParams *params = (DistancesParams*)loop.params;
	record = loop.blockIdx.x * kBlockSize + loop.threadIdx.x;

	// fill with first records from the list as candidates
	for ( j = 0; j < kMaxNeighbours; j++ ) {
		if ( j == record ) shift = 1;
		params->dataStore->neighbours[ record * kMaxNeighbours + j + shift] = j + shift;
		distancesArray[ j] = params->dataStore->distances[ DistanceVIdx( record, j + shift )];
	}
	
	// for each record in the data set
	for ( i = kMaxNeighbours; i < params->dataStore->info.numEntries; i++ ) {
		if ( record == i ) {
			// skip this one
			continue;
		}

		distance = params->dataStore->distances[ DistanceVIdx( record, i )];
		candidate = i;

		for ( j = 0; j < kMaxNeighbours; j++ ) {
			if ( params->dataStore->neighbours[ record * kMaxNeighbours + j] == i ) {
				// it's already here - break;
				break;
			}

			if ( distance < distancesArray[ j] ) {
				// continue with a little buble sorting
				a = params->dataStore->neighbours[ record * kMaxNeighbours + j];
				b = distancesArray[ j];
				params->dataStore->neighbours[ record * kMaxNeighbours + j] = candidate;
				distancesArray[ j] = distance;
				candidate = a;
				distance = b;
			}
		} // for each neighbour
	} // for each data entry
}