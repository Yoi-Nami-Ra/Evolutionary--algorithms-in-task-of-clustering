/**
 29/09/2011
 Jaroslaw Wojtasik

 noCuda

 dataLoader.c
 **/

#include "dataLoader.h"
#include "errors.h"
#include <stdlib.h>
#include <string.h>
//==============================================
//== Globals

static DataLoaderEntry ** gLoadersList = NULL; //< List of all available loaders
static unsigned char gLoadersListCount = 0; //< Number of loaders on the list
const unsigned char kStepSize = 4;

//==============================================
//== Functions

ErrorCode AddFunction( const char * name, LoaderFunc func ) {
	DataLoaderEntry* newEntry;
	int len;

	logDebug( "Adding Function (%s)", name );

	newEntry = (DataLoaderEntry*)malloc( sizeof(DataLoaderEntry) );
	checkAlloc( newEntry )
		return errNoMemory;
	}
	
	newEntry->loaderFunc = func;
	len = strlen( name );
	if ( len > 20 ) len = 20;
	memcpy( newEntry->name, name, len * sizeof(char) + 1 );
	newEntry->name[ len] = '\0';

	if ( gLoadersListCount % kStepSize == 0 ) {
		// resize
		gLoadersList = (DataLoaderEntry**)realloc( gLoadersList, sizeof(DataLoaderEntry*) * ( gLoadersListCount + kStepSize ) );
		checkAlloc( gLoadersList )
			free( newEntry );
			return errNoMemory;
		}
	}

	gLoadersList[ gLoadersListCount++] = newEntry;

	return errOk;
}
//----------------------------------------------

ErrorCode ListAllFunctions( unsigned char *num, char ***list ) {
	unsigned char i;

	*num = gLoadersListCount;
	list[ 0] = (char**)malloc( sizeof(char**) * gLoadersListCount );
	if ( *list == NULL ) {
		return SetLastErrorCode( errNoMemory );
	}
	for ( i = 0; i < (*num); i++ ) {
		list[ 0][ i] = gLoadersList[i]->name;
		if ( list[ 0][ i] == NULL ) {
			for (i--; i > 0; i--) {
				free( list[ 0][ i]);
			}
			return SetLastErrorCode( errNoMemory );
		}
	}
	return errOk;
}
//----------------------------------------------

ErrorCode LoaderDetails( unsigned char num, DataStore * store ) {
	LoaderFunc loaderFunction;

	if ( num >= gLoadersListCount ) {
		reportError( errOutOfRange, "with index:%d out of:%d",
			num, gLoadersListCount );
		return SetLastErrorCode( errOutOfRange );
	}

	loaderFunction = gLoadersList[ num]->loaderFunc;
	return loaderFunction( store, 0 );
}
//----------------------------------------------

const char * LoaderName( unsigned int num, ErrorCode * error ) {
	if ( num >= gLoadersListCount ) {
		reportError( errOutOfRange, "with index:%d out of:%d",
			num, gLoadersListCount );
		if ( error ) {
			*error = SetLastErrorCode( errOutOfRange );
		}
	}

	return gLoadersList[ num]->name;
}
//----------------------------------------------

ErrorCode GetDataFromLoader( unsigned int num, DataStore * dataStore ) {
	LoaderFunc loaderFunction;

	if ( num >= gLoadersListCount ) {
		reportError( errOutOfRange, "with index:%d out of:%d",
			num, gLoadersListCount );
		return SetLastErrorCode( errOutOfRange );
	}

	loaderFunction = gLoadersList[ num]->loaderFunc;
	return loaderFunction( dataStore, 1 );
}
//----------------------------------------------