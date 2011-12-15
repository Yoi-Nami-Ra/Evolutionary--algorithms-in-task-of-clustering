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
static const unsigned char kWineDimensions = 12; // acctualy 13, but first attribute holds class number
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

ErrorCode IrisLoaderFunc( DataStore * wineStore, char loadAll );

ErrorCode IrisLoaderFunc( DataStore * wineStore, char loadAll ) {
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

void SetupIrisLoader( void ) {
	AddFunction( kDataName, IrisLoaderFunc );
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

	wineStore->dataVector = (float*)malloc( kIrisEntries * kIrisDimensions * sizeof(float) );

	checkAlloc( wineStore->dataVector )
		return errNoMemory;
	}

	wineStore->classes = (unsigned int*)malloc( kIrisEntries * sizeof(unsigned int) );

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
	while( !feof( dataFile ) && ( index < kIrisEntries )) {
		read = fscanf( dataFile, "%3f,%3f,%3f,%3f,%s",
			wineStore->dataVector + index * kIrisDimensions,
			wineStore->dataVector + index * kIrisDimensions + 1,
			wineStore->dataVector + index * kIrisDimensions + 2,
			wineStore->dataVector + index * kIrisDimensions + 3, str );
		index++;
		wineStore->classes[ index] = ( index % 50 ) + 1;
		// check if we read 5 elements, 4 values and the rest
		if ( read != 5 ) {
			break;
		}
	}

	if ( read != 5 ) {
		// we failed to read the file
		reportError( errFileCorupted, "we've failed to read the data at index: %d.", index );
		free( wineStore->dataVector );
		fclose( dataFile );

		return SetLastErrorCode( errFileCorupted );
	}

	fclose( dataFile );

	return errOk;
}