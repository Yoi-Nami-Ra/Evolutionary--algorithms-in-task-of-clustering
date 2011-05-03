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

#endif //DATALOADER_H