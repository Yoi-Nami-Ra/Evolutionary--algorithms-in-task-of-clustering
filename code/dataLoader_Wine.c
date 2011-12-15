/**
 15/12/2011
 Jaroslaw Wojtasik

 noCuda

 dataLoader_Wine.c
 **/

#include "dataLoader_Wine.h"
#include "errors.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//==============================================
//== Globals

static const char * kDataName = "wine";
static const char* kDataFilePath = "./data/wine.data";
static const unsigned char kWineDimensions = 13;
static const unsigned int kWineEntries = 178;

//==============================================
//== Functions

/**
 * Loads data into the data store.
 * @param wineStore	- [in] pointer to structure where data should be stored.
 *
 * @return 
 */
static ErrorCode LoadData( DataStore * wineStore );

ErrorCode WineLoaderFunc( DataStore * wineStore, char loadAll );

ErrorCode WineLoaderFunc( DataStore * wineStore, char loadAll ) {
	ErrorCode err = errOk;

	wineStore->info.dimSize = kWineDimensions;
	wineStore->info.numEntries = kWineEntries;
	wineStore->dataVector = NULL;
	wineStore->distances = NULL;
	wineStore->neighbours = NULL;
	wineStore->info.name = strdup( kDataName );

	if ( !loadAll ) {
		return err;
	}

	// - Read records from file
	err = LoadData( wineStore );

	return err;
}
//----------------------------------------------

void SetupWineLoader( void ) {
	AddFunction( kDataName, WineLoaderFunc );
}

ErrorCode LoadData( DataStore * wineStore ) {
	FILE * dataFile = NULL;
	int index = 0;
	char str[20];
	char read = 0;

	if ( wineStore == NULL ) {
		reportError ( errWrongParameter, "Should be not NULL.%s", "" );
		return SetLastErrorCode( errWrongParameter );
	}

	wineStore->dataVector = (float*)malloc( kWineEntries * kWineDimensions * sizeof(float) );

	checkAlloc( wineStore->dataVector )
		return errNoMemory;
	}

	wineStore->classes = (unsigned int*)malloc( kWineEntries * sizeof(unsigned int) );

	checkAlloc( wineStore->classes )
		return errNoMemory;
	}
	
	// open file
	dataFile = fopen( kDataFilePath, "r" );
	if ( dataFile == 0 ) {
		free( wineStore->dataVector );
		reportError( errFileNotFound, "file:%s", kDataFilePath );
		return SetLastErrorCode( errFileNotFound );
	}

	// load data into memory
	while( !feof( dataFile ) && ( index < kWineEntries )) {
		// 1,14.23,1.71,2.43,15.6,127,2.8,3.06,.28,2.29,5.64,1.04,3.92,1065
		read = fscanf( dataFile, "%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
			wineStore->classes + index,
			wineStore->dataVector + index * kWineDimensions,
			wineStore->dataVector + index * kWineDimensions + 1,
			wineStore->dataVector + index * kWineDimensions + 2,
			wineStore->dataVector + index * kWineDimensions + 3,
			wineStore->dataVector + index * kWineDimensions + 4,
			wineStore->dataVector + index * kWineDimensions + 5,
			wineStore->dataVector + index * kWineDimensions + 6,
			wineStore->dataVector + index * kWineDimensions + 7,
			wineStore->dataVector + index * kWineDimensions + 8,
			wineStore->dataVector + index * kWineDimensions + 9,
			wineStore->dataVector + index * kWineDimensions + 10,
			wineStore->dataVector + index * kWineDimensions + 11,
			wineStore->dataVector + index * kWineDimensions + 12
			);
		index++;
		// check if we read 14 elements, 1 class, and 13 dimensions
		if ( read != 14 ) {
			break;
		}
	}

	if ( read != 14 ) {
		// we failed to read the file
		reportError( errFileCorupted, "we've failed to read the data at index: %d.", index );
		free( wineStore->dataVector );
		fclose( dataFile );

		return SetLastErrorCode( errFileCorupted );
	}

	fclose( dataFile );

	return errOk;
}