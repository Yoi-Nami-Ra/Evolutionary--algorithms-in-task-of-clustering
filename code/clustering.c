/**
 17/11/2011
 Jaroslaw Wojtasik

 noCuda

 clustering.h
**/

#include "clustering.h"
#include "errors.h"
#include "loops.h"
#include "distanceCalculator.h"
#include <stdlib.h>
#include <time.h>

//==============================================
//== Globals

#define MAX_CLUSTER_SIZE 10
#define MEDOID_VECTOR_SIZE 20
#define CROS_FACTOR 3
#define OBJECTIVES 4

const char threadsPerBlock = 50;

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

/**
 * Prepare enviroment to run algorithms.
 */
ErrorCode ConfigureAlgoritms( EvolutionProps * props );

ErrorCode MembershipAndDensity( EvolutionProps * props );

ErrorCode Connectivity( EvolutionProps * props );

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
	EvolutionProps * props = (EvolutionProps*)loop.params;
	char proposalOk;
	unsigned int proposal;
	unsigned int clustersSum = MEDOID_VECTOR_SIZE;
	int i,j;

	// local max size of a cluster
	props->population[ loop.threadIdx.x].attr.clusterMaxSize = rand() % MAX_CLUSTER_SIZE + 1;
	// local max list of neighbours
	props->population[ loop.threadIdx.x].attr.numNeighbours = rand() % MAX_NEIGHBOURS + 1;

	// for each medoid in vector
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		do {
			proposalOk = 1;
			// chose random data entry
			proposal = rand() % props->dataStore->info.numEntries;
			// and check if we didn't that already
			for ( j = 0; j < i; j++ ) {
				if ( proposal == props->population[ loop.threadIdx.x].medoids[ j] ) {
					proposalOk = 0;						
				}
			}
		} while ( !proposalOk );
		// save it
		props->population[ loop.threadIdx.x].medoids[ i] = proposal;

		// now clusters
		// if we didn't fill it up yet
		if ( clustersSum > 0 ) {
			proposal = rand() % props->population[ loop.threadIdx.x].attr.clusterMaxSize + 1;
			if ( clustersSum < proposal ) {
				proposal += clustersSum - proposal;
			}
			clustersSum -= proposal;
			props->population[ loop.threadIdx.x].clusters[ i] = proposal;
		} else {
			props->population[ loop.threadIdx.x].clusters[ i] = 0;
		}
	} // i
}
//----------------------------------------------

ErrorCode GenerateRandomPopulation( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition popGenLoop;

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

	popGenLoop.blockSize.x = props->popSize;
	popGenLoop.blockSize.y = popGenLoop.gridSize.x = 
		popGenLoop.gridSize.y = popGenLoop.blockSize.z =
		popGenLoop.gridSize.z = 1;
	popGenLoop.kernel = GenerateRandomPopulationKernel;
	popGenLoop.params = (void*)&props;

	err = RunLoop( popGenLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

ErrorCode RunAlgorithms( EvolutionProps * props ) {
	ErrorCode err = errOk;
	unsigned int i;

	err = ConfigureAlgoritms( props );
	if ( err != errOk ) {
		return err;
	}

	for ( i = 0; i < props->evoSteps; i++ ) {
		// mambership and density
		err = MembershipAndDensity( props );
		if ( err != errOk ) {
			break;
		}

		// connectivity

		// sum up results

		// disconnectivity

		// correctness

		// sorting

		// dominance count

		// setup fronts

		// group fronts

		// selection

		// 



	}

	return err;
}
//----------------------------------------------

ErrorCode ConfigureAlgoritms( EvolutionProps * props ) {
	unsigned int i;

	props->solutions = (Solution*)malloc( props->popSize * sizeof(Solution) );

	checkAlloc( props->solutions )
		return GetLastErrorCode();
	}
	// clear all densities
	for ( i = 0; i < props->popSize; i++ ) {
		props->solutions[ i].densities = 0;
	}

	props->blocksPerEntries = props->dataStore->info.numEntries / threadsPerBlock;
	while ( props->blocksPerEntries * threadsPerBlock <  props->dataStore->info.numEntries ) {
		props->blocksPerEntries++;
	}
}
//----------------------------------------------
void MembershipAndDensityKernel( LoopContext loop ) {
	// for this record in this specific solution (population member)
	unsigned int solution = loop.blockIdx.y;
	unsigned int record = loop.blockIdx.x + loop.threadIdx.x;
	unsigned int i;
	EvolutionProps * props = (EvolutionProps*)loop.params;
	float currentDistance = 0;
	float smallestDistance = props->dataStore->distances[ DistanceVIdx( record, 0 )];
	props->solutions[ solution].membership[ record] = 0;

	// find closest medoid to this record
	for ( i = 1; i < MEDOID_VECTOR_SIZE; i++ ) {
		 currentDistance = props->dataStore->distances[ DistanceVIdx( record, 0 )];
		 if ( currentDistance < smallestDistance ) {
			 smallestDistance = currentDistance;
			 props->solutions[ solution].membership[ record] = i;
		 }
	}

	
	props->solutions[ solution].densities += smallestDistance;
}
//----------------------------------------------

ErrorCode MembershipAndDensity( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition densityLoop;

	if ( props == NULL || props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	//<<( blocksPerEntires, populationSize ), threadsPerBlock>>
	densityLoop.gridSize.x = props->blocksPerEntries;
	densityLoop.gridSize.y = props->popSize;
	densityLoop.gridSize.z = 1;
	densityLoop.blockSize.x = threadsPerBlock;
	densityLoop.blockSize.y = densityLoop.blockSize.z = 1;
	densityLoop.kernel = MembershipAndDensityKernel;
	densityLoop.params = (void*)&props;

	err = RunLoop( densityLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

void ConnectivityKernel( LoopContext loop ) {
	unsigned int solution = loop.blockIdx.y;
	unsigned int record = loop.blockIdx.x + loop.threadIdx.x;
}
//----------------------------------------------
ErrorCode Connectivity( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition connectivityLoop;

	if ( props == NULL || props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	//<<( blocksPerEntires, populationSize ), threadsPerBlock>>
	connectivityLoop.gridSize.x = props->blocksPerEntries;
	connectivityLoop.gridSize.y = props->popSize;
	connectivityLoop.gridSize.z = 1;
	connectivityLoop.blockSize.x = threadsPerBlock;
	connectivityLoop.blockSize.y = connectivityLoop.blockSize.z = 1;
	connectivityLoop.kernel = ConnectivityKernel;
	connectivityLoop.params = (void*)&props;
}
//----------------------------------------------
//----------------------------------------------