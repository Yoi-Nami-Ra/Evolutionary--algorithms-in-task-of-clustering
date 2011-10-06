/**
 29/09/2011
 Jaroslaw Wojtasik

 noCuda

 dataLoader.cpp
 **/

#include "dataLoader.h"
#include "errors.h"
#include <stdlib.h>
#include <string.h>

//==============================================
//== Globals

static unsigned char gNumLoaders = 0; //< number of all loaders in the list

static DataLoaderEntry* gLoadersList = NULL; //< List of all available loaders

//==============================================
//== Functions

ErrorCode AddFunction( const char * name, LoaderFunc * func ) {
	logDebug( "Adding Function (%s)", name );

	DataLoaderEntry *newEntry = malloc( sizeof(DataLoaderEntry) );
	newEntry->func = func;
	int len = strlen( name );
	if ( len > 20 ) len = 20;
	newEntry->name = malloc( len * sizeof(char) );
	memcpy( newEntry->name, name, len * sizeof(char) );
	
	newEntry->name = 

	if ( gNumLoaders % 4 ) {
		// Some place left
		gLoadersList[ gNumLoaders ];

	} else {
	}
}
