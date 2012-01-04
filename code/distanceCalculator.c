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

const unsigned char kBlockSize = 16;
static const unsigned char kMaxNeighbours = MAX_NEIGHBOURS;

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
char CheckHaveSavedDistances( DataStore * dataStore );

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

/**
 * Kernel function to run inside the loop
 * @param loop - [in] structure that holds all context of the loop.
 */
void CalculateNeighboursKernel( LoopContext loop );
//============================================================================


ErrorCode GetCalculatedDistances( unsigned int num, DataStore * dataStore ) {
	ErrorCode err;
	char preCalculated = 0;

	err = GetDataFromLoader( num, dataStore );

	if ( err != errOk ) {
		return err;
	}

	preCalculated = CheckHaveSavedDistances( dataStore );
	if ( preCalculated ) {
		err = LoadCalculatedDistances( dataStore );
	}

	if ( !preCalculated || err != errOk )  {
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

char CheckHaveSavedDistances( DataStore * dataStore ) {
	unsigned int nameLen = 0;
	char * fileName = NULL;
	FILE * file = NULL;

	nameLen = (unsigned int)strlen( dataStore->info.name );
	nameLen += (unsigned int)strlen( "_distances.data" );
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
	LoopDefinition distanceLoop;
	DistancesParams params;
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
		reportError( err, "Run loop returned with error%s", "" );
	}

	return err;
}
//----------------------------------------------

unsigned int DistanceVIdx( unsigned int a, unsigned int b ) {
	if ( a > b ) {
		return a * (a - 1) / 2 + b;
	} else {		
		return b * (b - 1) / 2 + a;
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

ErrorCode LoadCalculatedDistances( DataStore * dataStore ) {
	// NAME_distances.data
	unsigned int nameLen = 0;
	char * fileName = NULL;
	FILE * file = NULL;
	unsigned int read = 0;
	unsigned int i = 0;

	if ( dataStore == NULL ||
		dataStore->info.name == NULL ) {
			reportError( errWrongParameter, "Got wrong parameters dataStore or name are NULLs %s", "" );
			return errWrongParameter;
	}

	nameLen = (unsigned int)strlen( dataStore->info.name );
	nameLen += (unsigned int)strlen( "_distances.data" );
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
	read = (unsigned int)fread( dataStore->distances, sizeof(float), dataStore->info.distancesSize, file );
	if ( read != dataStore->info.distancesSize ) {
		reportError( errFileRead, "Failed to read distances data properly.%s", "" );
		return errFileRead;
	}

	dataStore->neighbours = (unsigned int*)malloc( dataStore->info.numEntries * kMaxNeighbours * sizeof(unsigned int) );
	checkAlloc( dataStore->neighbours )
		return GetLastErrorCode();
	}
	read = (unsigned int)fread( dataStore->neighbours, kMaxNeighbours * sizeof(unsigned int), dataStore->info.numEntries, file );
	if ( read != dataStore->info.numEntries ) {
		reportError( errFileRead, "Failed to read neighbours data properly.%s", "" );
		return errFileRead;
	}

	fclose( file );

	file = fopen( "distances.txt", "w" );
	for ( i = 0; i < dataStore->info.distancesSize; i++ ) {
		fprintf( file, " %f\n", dataStore->distances[ i] );
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
			reportError( errWrongParameter, "Got wrong parameters dataStore %s", "" );
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
	LoopDefinition neighbourLoop;
	DistancesParams params;
	ErrorCode err = errOk;

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
		reportError( err, "Run loop returned with error%s", "" );
	}

    return err;
}
//----------------------------------------------

void CalculateNeighboursKernel( LoopContext loop ) {
	unsigned int record = 0;
	unsigned int i = 0;
	unsigned int j = 0;
	float distance = 0;
	float distancesArray[ MAX_NEIGHBOURS];
	unsigned int candidate = 0;
	unsigned int a = 0;
	float b = 0;

	DistancesParams *params = (DistancesParams*)loop.params;
	record = loop.blockIdx.x * kBlockSize + loop.threadIdx.x;

	if ( record >= params->dataStore->info.numEntries ) {
		return;
	}

	// fill with first records from the list as candidates
	for ( j = 0; j < kMaxNeighbours; j++ ) {
		distancesArray[ j] = 0;
	}
	
	// for each record in the data set
	for ( i = 0; i < params->dataStore->info.numEntries; i++ ) {
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

			if ( distancesArray[ j] == 0 ) {
				params->dataStore->neighbours[ record * kMaxNeighbours + j] = candidate;
				distancesArray[ j] = distance;
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