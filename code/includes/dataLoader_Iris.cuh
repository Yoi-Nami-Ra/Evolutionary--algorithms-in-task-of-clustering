// Header file for data loader of Iris database

#ifndef DATALOADER_IRIS_H
#define DATALOADER_IRIS_H

// Includes
#include <stdlib.h>
#include <stdio.h>
#include "errors.cuh"

// Data Types

typedef struct dataInfo {
	unsigned int numEntries; //< How Many etries in the data store
} dataInfo;

// specialization for Iris
typedef struct dataEntry {
	float a;
	float b;
	float c;
	float d;
} dataEntry;

typedef struct dataStore {
	dataInfo info;
	dataEntry* dataVector;
} dataStore;

// Functions

/*
 * Check if we didn't made data conversion already and saved it.
 * @ret true if the file already exists
 */
bool CheckAlreadyConverted();

/*
 * Start Loading Data. 
 */
 ErrorCode startLoadingData();

 /*
  * Returns pointer to current data store.
  */
 dataStore* GetCurrDataStore();

 /*
  * Cleans up the memory
  */
 ErrorCode releaseDataStore();


#endif // DATALOADER_IRIS_H