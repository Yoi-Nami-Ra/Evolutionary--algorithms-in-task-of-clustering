/**
 17/11/2011
 Jaroslaw Wojtasik

 noCuda

 clustering.h
**/

#include "clustering.h"
#include "errors.h"
#include "loops.h"
#include <stdlib.h>
#include <time.h>

//==============================================
//== Globals

#define MAX_CLUSTER_SIZE 10
#define MEDOID_VECTOR_SIZE 20
#define CROS_FACTOR 3
#define OBJECTIVES 4

//==============================================
//== Functions

/**
 * Generates random population to start with.
 * @param props	- [in, out] For parameters and results of the operation.
 */
ErrorCode GenerateRandomPopulation( EvolutionProps * props );

/**
 * Main algorithms.
 */ 
ErrorCode RunAlgorithms( EvolutionProps * props );

//=============================================

/**
 * Just some defaults.
 */
ErrorCode GenerateDefaultProps( EvolutionProps * props ) {
//TODO: implement
	return errOk;
}
//----------------------------------------------

ErrorCode RunClustering( EvolutionProps * props ) {
	ErrorCode err = errOk;

	err = GenerateRandomPopulation( props );
	if ( err != errOk ) {
		reportError( err, "Failed to generate random population%s", "" );
		return err;
	}

	err = RunAlgorithms( props );

	if ( err != errOk ) {
		reportError( err, "Failed with algorithms%s", "" );
	}

	return err;
}
//----------------------------------------------

void GenerateRandomPopulationKernel( LoopContext loop ) {
}
//----------------------------------------------

ErrorCode GenerateRandomPopulation( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition distanceLoop;

	if ( props == NULL ) {
		reportError( errWrongParameter, "Evolution Propertis served with NULL pointer.%s", "" );
	}

	if ( props->population != NULL ) {
		free( props->population );
	}

	props->population = (PopMember*)malloc( props->popSize * sizeof(PopMember) );
	checkAlloc( props->population )
		return GetLastErrorCode();
	}

	srand( (unsigned int)time( NULL ) );

	distanceLoop.blockSize.x = props->popSize;
	distanceLoop.blockSize.y = distanceLoop.gridSize.x = 
		distanceLoop.gridSize.y = distanceLoop.blockSize.z =
		distanceLoop.gridSize.z = 1;
	distanceLoop.kernel = GenerateRandomPopulationKernel;
	distanceLoop.params = (void*)&props;

	err = RunLoop( distanceLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

ErrorCode RunAlgorithms( EvolutionProps * props ) {
//TODO: implement
	return errOk;
}
//----------------------------------------------
//----------------------------------------------