/**
 18/10/2011
 Jaroslaw Wojtasik

 noCuda

 dataLoader_Iris.c
 **/

#include "dataLoader_Iris.h"
#include "errors.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//==============================================
//== Globals

static const char * kDataName = "iris";
static const char* kDataFilePath = "./data/iris.data";
static const unsigned char kIrisDimensions = 4;
static const unsigned int kIrisEntries = 150;

//==============================================
//== Functions

/**
 * Loads data into the data store.
 * @param irisStore	- [in] pointer to structure where data should be stored.
 *
 * @return 
 */
static ErrorCode LoadData( DataStore * irisStore );

ErrorCode IrisLoaderFunc( DataStore * irisStore, char loadAll );

ErrorCode IrisLoaderFunc( DataStore * irisStore, char loadAll ) {
	ErrorCode err = errOk;

	irisStore->info.dimSize = kIrisDimensions;
	irisStore->info.numEntries = kIrisEntries;
	irisStore->dataVector = NULL;
	irisStore->distances = NULL;
	irisStore->neighbours = NULL;
	irisStore->info.name = strdup( kDataName );

	if ( !loadAll ) {
		return err;
	}

	// - Read records from file
	err = LoadData( irisStore );

	return err;
}
//----------------------------------------------

void SetupIrisLoader( void ) {
	AddFunction( kDataName, IrisLoaderFunc );
}

ErrorCode LoadData( DataStore * irisStore ) {
	FILE * dataFile = NULL;
	int index = 0;
	char str[20];
	char read = 0;

	if ( irisStore == NULL ) {
		reportError ( errWrongParameter, "Should be not NULL.%s", "" );
		return SetLastErrorCode( errWrongParameter );
	}

	irisStore->dataVector = (float*)malloc( kIrisEntries * kIrisDimensions * sizeof(float) );

	checkAlloc( irisStore->dataVector )
		return errNoMemory;
	}

	irisStore->classes = (unsigned int*)malloc( kIrisEntries * sizeof(unsigned int) );

	checkAlloc( irisStore->classes )
		return errNoMemory;
	}
	
	// open file
	dataFile = fopen( kDataFilePath, "r" );
	if ( dataFile == 0 ) {
		free( irisStore->dataVector );
		reportError( errFileNotFound, "file:%s", kDataFilePath );
		return SetLastErrorCode( errFileNotFound );
	}

	// load data into memory
	while( !feof( dataFile ) && ( index < kIrisEntries )) {
		read = fscanf( dataFile, "%3f,%3f,%3f,%3f,%s",
			irisStore->dataVector + index * kIrisDimensions,
			irisStore->dataVector + index * kIrisDimensions + 1,
			irisStore->dataVector + index * kIrisDimensions + 2,
			irisStore->dataVector + index * kIrisDimensions + 3, str );
		index++;
		irisStore->classes[ index] = ( index % 50 ) + 1;
		// check if we read 5 elements, 4 values and the rest
		if ( read != 5 ) {
			break;
		}
	}

	if ( read != 5 ) {
		// we failed to read the file
		reportError( errFileCorupted, "we've failed to read the data at index: %d.", index );
		free( irisStore->dataVector );
		fclose( dataFile );

		return SetLastErrorCode( errFileCorupted );
	}

	fclose( dataFile );

	return errOk;
}