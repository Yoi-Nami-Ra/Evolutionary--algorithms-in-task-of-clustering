/**
 29/09/2011
 Jaroslaw Wojtasik

 noCuda

 dataLoader.h
 **/

#ifndef DATALOADER_H
#define DATALOADER_H

#include "errors.h"

//==============================================
//== Types

typedef struct DataInfo {
	unsigned char	dimSize; //< How many dimensions
	unsigned int	numEntries; //< How Many etries in the data store
} DataInfo;

typedef struct DataStore {
	DataInfo	info; //< description of data held
	float*		dataVector; //< data
} DataStore;

typedef DataStore* (*LoaderFunc)();

typedef struct DataLoaderEntry {
	char		name[21]; //< name of the data source
	LoaderFunc	func; //< function to load data
} DataLoaderEntry;

//==============================================
//== Functions

ErrorCode AddFunction(const char * name, LoaderFunc * func);



#endif // DATALOADER_H