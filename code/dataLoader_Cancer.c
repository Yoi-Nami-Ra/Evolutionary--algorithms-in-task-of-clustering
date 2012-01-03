//
//  dataLoader_Cancer.c
//  NoCuda
//
//  Created by Jarosław Wojtasik on 02.01.2012.
//

#include "dataLoader_Cancer.h"
#include "errors.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//==============================================
//== Globals

static const char * kDataName = "cancer";
static const char* kDataFilePath = "./data/wpbc.data";
static const unsigned char kCancerDimensions = 32; // first "ID" to be ignored, second holds class, the rest used for algorithms
static const unsigned int kCancerEntries = 198;

//==============================================
//== Functions

/**
 * Loads data into the data store.
 * @param cancerStore	- [in] pointer to structure where data should be stored.
 *
 * @return 
 */
static ErrorCode LoadData( DataStore * cancerStore );

ErrorCode CancerLoaderFunc( DataStore * cancerStore, char loadAll );

ErrorCode CancerLoaderFunc( DataStore * cancerStore, char loadAll ) {
	ErrorCode err = errOk;
    
	cancerStore->info.dimSize = kCancerDimensions;
	cancerStore->info.numEntries = kCancerEntries;
	cancerStore->dataVector = NULL;
	cancerStore->distances = NULL;
	cancerStore->neighbours = NULL;
	cancerStore->info.name = strdup( kDataName );
    
	if ( !loadAll ) {
		return err;
	}
    
	// - Read records from file
	err = LoadData( cancerStore );
    
	return err;
}
//----------------------------------------------

void SetupCancerLoader( void ) {
	AddFunction( kDataName, CancerLoaderFunc );
}

ErrorCode LoadData( DataStore * cancerStore ) {
	FILE * dataFile = NULL;
	int index = 0;
	char read = 0;
    unsigned int id;
    char class;
    
	if ( cancerStore == NULL ) {
		reportError ( errWrongParameter, "Should be not NULL.%s", "" );
		return SetLastErrorCode( errWrongParameter );
	}
    
	cancerStore->dataVector = (float*)malloc( kCancerEntries * kCancerDimensions * sizeof(float) );
    
	checkAlloc( cancerStore->dataVector )
        return errNoMemory;
    }

    cancerStore->classes = (unsigned int*)malloc( kCancerEntries * sizeof(unsigned int) );

    checkAlloc( cancerStore->classes )
        return errNoMemory;
    }

    // open file
    dataFile = fopen( kDataFilePath, "r" );
    if ( dataFile == 0 ) {
        free( cancerStore->dataVector );
        reportError( errFileNotFound, "file:%s", kDataFilePath );
        return SetLastErrorCode( errFileNotFound );
    }

    // load data into memory
    while( !feof( dataFile ) && ( index < kCancerEntries )) {
        // 1,14.23,1.71,2.43,15.6,127,2.8,3.06,.28,2.29,5.64,1.04,3.92,1065
        read = fscanf( dataFile, "%u,%c,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
                      &id,
                      &class,
                      cancerStore->dataVector + index * kCancerDimensions,
                      cancerStore->dataVector + index * kCancerDimensions + 1,
                      cancerStore->dataVector + index * kCancerDimensions + 2,
                      cancerStore->dataVector + index * kCancerDimensions + 3,
                      cancerStore->dataVector + index * kCancerDimensions + 4,
                      cancerStore->dataVector + index * kCancerDimensions + 5,
                      cancerStore->dataVector + index * kCancerDimensions + 6,
                      cancerStore->dataVector + index * kCancerDimensions + 7,
                      cancerStore->dataVector + index * kCancerDimensions + 8,
                      cancerStore->dataVector + index * kCancerDimensions + 9,
                      cancerStore->dataVector + index * kCancerDimensions + 10,
                      cancerStore->dataVector + index * kCancerDimensions + 11,
                      cancerStore->dataVector + index * kCancerDimensions + 12,
                      cancerStore->dataVector + index * kCancerDimensions + 13,
                      cancerStore->dataVector + index * kCancerDimensions + 14,
                      cancerStore->dataVector + index * kCancerDimensions + 15,
                      cancerStore->dataVector + index * kCancerDimensions + 16,
                      cancerStore->dataVector + index * kCancerDimensions + 17,
                      cancerStore->dataVector + index * kCancerDimensions + 18,
                      cancerStore->dataVector + index * kCancerDimensions + 19,
                      cancerStore->dataVector + index * kCancerDimensions + 20,
                      cancerStore->dataVector + index * kCancerDimensions + 21,
                      cancerStore->dataVector + index * kCancerDimensions + 22,
                      cancerStore->dataVector + index * kCancerDimensions + 23,
                      cancerStore->dataVector + index * kCancerDimensions + 24,
                      cancerStore->dataVector + index * kCancerDimensions + 25,
                      cancerStore->dataVector + index * kCancerDimensions + 26,
                      cancerStore->dataVector + index * kCancerDimensions + 27,
                      cancerStore->dataVector + index * kCancerDimensions + 28,
                      cancerStore->dataVector + index * kCancerDimensions + 29,
                      cancerStore->dataVector + index * kCancerDimensions + 30,
                      cancerStore->dataVector + index * kCancerDimensions + 31
                      );
        cancerStore->classes[ index] = class == 'R'?1:0; 
        index++;
        // check if we read 14 elements, 1 class, and 13 dimensions
        if ( read != 34 && read != 0 ) {
            break;
        }
    }

if ( read != 34 && read != 0 ) {
    // we failed to read the file
    reportError( errFileCorupted, "we've failed to read the data at index: %d.", index );
    free( cancerStore->dataVector );
    fclose( dataFile );
    
    return SetLastErrorCode( errFileCorupted );
}

fclose( dataFile );

return errOk;
}
