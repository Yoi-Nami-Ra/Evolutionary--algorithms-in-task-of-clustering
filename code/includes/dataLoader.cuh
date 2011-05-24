// Header file for dataLoader

#ifndef DATALOADER_H
#define DATALOADER_H


#include "globals.cuh"

// Here put includes with data loader specializations for all data source types
#ifdef DATABASE_IRIS
#include "dataLoader_Iris.cuh"
#endif // DATABASE_IRIS

/*
 * Loads data.
 */
ErrorCode LoadData();

/*
 * returns number of entries in data store.
 */
unsigned int numEntries();

/*
 * Set current number of entries
 */
void setNumEntries( unsigned int numEntries );


#endif //DATALOADER_H