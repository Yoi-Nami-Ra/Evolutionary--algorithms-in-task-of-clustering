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

/**
 * Used to hold description of specific data store.
 */
typedef struct {
	unsigned char	dimSize;		///< How many dimensions
	unsigned int	numEntries;		///< How Many etries in the data store
	unsigned int	distancesSize;	///< Holds the size of array that holds calculated distances
	char			*name;			///< Name of the data store
} DataInfo;

/**
 * Represents specific data store and all data associated with it trough the program run.
 */
typedef struct {
	DataInfo		info;			///< description of data held
	float			*dataVector;	///< data
	float			*distances;		///< array to hold calculated distances
	unsigned int	*neighbours;	///< array of closest neighbours
	unsigned int	*classes;		///< storres assigment to predefined classes for RAND index
} DataStore;

/**
 * Function Prototype to generate Data store.
 * @param store Pointer to data store that should be returned
 * @param loadAll if true, description  an all data records are loaded.
 *
 * @return ErrorCode if any.
 */
typedef ErrorCode (*LoaderFunc)( DataStore*, char );

/**
 * Used to describe single Loader function.
 */
typedef struct {
	char		name[21]; //< name of the data source
	LoaderFunc loaderFunc;
} DataLoaderEntry;

//==============================================
//== Functions

/**
 * Adds Loader function to the list.
 */
ErrorCode AddFunction( const char * name, LoaderFunc func );

/**
 * Lists all loader functions
 */
ErrorCode ListAllFunctions( unsigned char *num, char ***list );

/**
 * Returns dataStore structure filled with details, but not with data.
 * @param num Number of the loader from the list.
 * @param store pointer to dataStore structure. Must be freed.
 * @return Error code if any
 */
ErrorCode LoaderDetails( unsigned char num, DataStore * store );

/**
 * Returns name of the given loader.
 * @param num[in] Number of the loader.
 * @param error[out] Pointer where error code could be stored.
 *
 * @return const char Pointer to string with name.
 */
const char * LoaderName( unsigned int num, ErrorCode * error );

/**
 * Returns data store struct filled with converted data.
 * @param num		- [in] index of the selected loader
 * @param dataStore	- [out] pointer where data would be stored
 *
 * @return error code if any
 */
ErrorCode GetDataFromLoader( unsigned int num, DataStore * dataStore );

#endif // DATALOADER_H