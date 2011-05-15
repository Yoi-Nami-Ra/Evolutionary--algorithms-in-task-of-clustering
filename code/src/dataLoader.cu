// Module responsible for loading data

//====================================================================
//== Includes
#include <stdio.h>
#include "globals.cuh"
#include "errors.cuh"
#include "dataLoader.cuh"

//====================================================================
//== Types

//====================================================================
//== Constants and Globals
static unsigned int dataStoreEntriesNumber;

//====================================================================
//== Functions

ErrorCode LoadData() {
	ErrorCode err = startLoadingData();
	if ( err == errOk ) {
		dataStoreEntriesNumber = GetCurrDataStore()->info.numEntries;
		if ( dataStoreEntriesNumber == 0 ) {
			printf("[E] LoadData(): 0 Entries after succesgul load\n");
			err = errGeneral;
		}
	}
	return err;
}
//====================================================================

unsigned int numEntries() {
	return dataStoreEntriesNumber;
}
//====================================================================