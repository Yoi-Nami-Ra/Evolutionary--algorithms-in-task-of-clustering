//
//  dataLoader_Cancer.cu
//  Cuda
//
//  Created by Jarosław Wojtasik on 02.01.2012.
//

#include "dataLoader_Cancer.cuh"
#include "errors.cuh"
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
    char tClass;
	char rest[ 20];

    
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
        // 119513,N,31,18.02,27.6,117.5,1013,0.09489,0.1036,0.1086,0.07055,0.1865,0.06333,0.6249,1.89,3.972,71.55,0.004433,0.01421,0.03233,0.009854,0.01694,0.003495,21.63,37.08,139.7,1436,0.1195,0.1926,0.314,0.117,0.2677,0.08113,5,5
        read = fscanf( dataFile, "%u,%c,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s",
                      &id,
                      &tClass,
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
                      cancerStore->dataVector + index * kCancerDimensions + 31,
					  rest
                      );
        cancerStore->classes[ index] = tClass == 'R'?1:0; 
        index++;
        // check if we read 14 elements, 1 class, and 13 dimensions
        if ( read != 35 && read != 34 ) {
            break;
        }
    }

if ( read != 35 && read != 34 ) {
    // we failed to read the file
    reportError( errFileCorupted, "we've failed to read the data at index: %d.", index );
    free( cancerStore->dataVector );
    fclose( dataFile );
    
    return SetLastErrorCode( errFileCorupted );
}

fclose( dataFile );

return errOk;
}
