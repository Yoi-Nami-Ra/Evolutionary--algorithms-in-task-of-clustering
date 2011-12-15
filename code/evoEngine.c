/**
 17/11/2011
 Jaroslaw Wojtasik

 noCuda

 evoEngine.c
 **/

#include "evoEngine.h"
#include "dataLoader.h"
#include "dataLoader_Iris.h"
#include "dataLoader_Test.h"
#include "distanceCalculator.h"
#include "clustering.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * Displays list of available loaders.
 * @return (unsigned int) number of elements on the list
 */
unsigned char DisplayLoadersList( void );
unsigned char DisplayLoadersList( void ) {
	ErrorCode err = errOk;
	unsigned char num = 0;
	char ** list = NULL;
	int i = 0;

	// 1) Get the list
	err = ListAllFunctions( &num, &list );
	if ( err != errOk ) {
		reportError( err, "Error while fetching list of loaders%s.", "" );
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
unsigned int PromptForSelection( void );
unsigned int PromptForSelection( void ) {
	unsigned int result;

	printf( "\n >> ");
	scanf( "%i", &result );

	return result;
}

/**
 * Print Details
 */
void PrintDetails( unsigned int selected );
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
char ConfirmSelection( void );
char ConfirmSelection( void ) {
	int result;

	printf( "\n [y,n]>> ");
	do {
		result = getchar();
	} while ( result != 'y' && result != 'n' );

	return ( result == 'y' );
}

void runEvo( void ) {
	unsigned int selectedLoader;
	unsigned char loadersCount;
	int b = 0;
	ErrorCode err = errOk;
	DataStore dataStore;
	EvolutionProps props;

	// -- Setting up all loaders
	SetupIrisLoader();
	SetupTestLoader();

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
        
        if ( err != errOk ) {
            // error occured can't continue with the algorithms
            printf( " Error occured while preparing data for algorithms" );
        } else {
            // -- Run Algorithms
            props.blocksPerEntries = 0;
            props.crosFactor = 0;
            props.dataStore = &dataStore;
            props.dominanceCounts = NULL;
            props.dominanceMatrix = NULL;
            props.evoSteps = 3000; // 100 steps for the alg
            props.popSize = 4; // for now 256 members
            props.population = NULL;
            props.solutions = NULL;
            err = RunClustering( &props );
        }

		b++;
	}


	// -- 
	return;
}