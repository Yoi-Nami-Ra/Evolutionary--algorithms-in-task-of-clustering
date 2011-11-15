// Evo_vs2010_noCUDA.cpp : Defines the entry point for the console application.
//

#include <tchar.h>
#include "dataLoader.h"
#include "dataLoader_Iris.h"
#include "distanceCalculator.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * Displays list of available loaders.
 * @return (unsigned int) number of elements on the list
 */
unsigned char DisplayLoadersList() {
	ErrorCode err = errOk;
	unsigned char num = 0;
	char ** list = NULL;
	int i = 0;

	// 1) Get the list
	err = ListAllFunctions( &num, &list );
	if ( err != errOk ) {
		reportError( err );
	}
	// 2) Print
	for ( i = 0; i < num; i++ ) {
		printf( " %d : %s\n", i+1, list[ i] );
	}

	return num;
}

/**
 * Ask user for selection
 */
unsigned int PromptForSelection() {
	unsigned int result;

	printf( "\n >> ");
	scanf( "%i", &result );

	return result;
}

/**
 * Print Details
 */
void PrintDetails( unsigned int selected ) {
	DataStore dataStore;
	const char * loaderName;

	LoaderDetails( selected-1, &dataStore );
	loaderName = LoaderName( selected-1, NULL );

	printf( "Loader %i\n Name: %s\n Num Entries: %i\n Dimensions: %i\n",
		selected, loaderName, dataStore.info.numEntries, dataStore.info.dimSize );
}

/**
 * Confirm selection of this loader.
 * 
 * @return True if the user selected it.
 */
char ConfirmSelection() {
	int result;

	printf( "\n [y,n]>> ");
	do {
		result = getchar();
	} while ( result != 'y' && result != 'n' );

	return ( result == 'y' );
}

int _tmain(int argc, _TCHAR* argv[])
{
	unsigned int selectedLoader;
	unsigned char loadersCount;
	ErrorCode err = errOk;
	DataStore dataStore;

	// -- Setting up all loaders
	SetupIrisLoader();

	for (;;) {

		// -- Print list of available loaders
		loadersCount = DisplayLoadersList();

		// -- Ask the user for the list number
		selectedLoader = PromptForSelection();

		if ( selectedLoader > loadersCount || selectedLoader <= 0 ) {
			// -- Wrong loader selected
			printf(" Wrong loader selected \n Choose [1 - %i]\n\n", loadersCount);
			continue;
		}

		// -- Show selected loader details
		PrintDetails( selectedLoader );

		// -- Confirm
		if ( !ConfirmSelection() ) {
			continue;
		}

		// -- Load
		err = GetCalculatedDistances( selectedLoader - 1, &dataStore );

	}


	// -- 
	return 0;
}

