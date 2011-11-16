/**
 16/11/2011
 Jaroslaw Wojtasik

 noCuda

 dataLoader_test.c
 **/

#include "errors.h"
#include "dataLoader_Test.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//==============================================
//== Globals

static const char * kDataName = "test";
static const unsigned char kTestDimensions = 2;
static const unsigned int kTestEntries = 6;
static float kData[] = { 3.0, 2.0,
	4.0, 5.0,
	7.0, 2.0,
	5.0, 9.0,
	7.0, 7.0,
	3.0, 7.0 };


//==============================================
//== Functions

/**
 * Loads data into the data store.
 * @param testStore	- [in] pointer to structure where data should be stored.
 *
 * @return 
 */
static ErrorCode LoadData( DataStore * testStore );

ErrorCode TestLoaderFunc( DataStore * testStore, char loadAll ) {
	ErrorCode err = errOk;

	testStore->info.dimSize = kTestDimensions;
	testStore->info.numEntries = kTestEntries;
	testStore->dataVector = NULL;
	testStore->distances = NULL;
	testStore->neighbours = NULL;
	testStore->info.name = strdup( kDataName );

	if ( !loadAll ) {
		return err;
	}

	// - Read records from file
	err = LoadData( testStore );

	return err;
}
//----------------------------------------------

void SetupTestLoader() {
	AddFunction( kDataName, TestLoaderFunc );
}

ErrorCode LoadData( DataStore * testStore ) {
	int i = 0;

	testStore->dataVector = kData;

	return errOk;
}