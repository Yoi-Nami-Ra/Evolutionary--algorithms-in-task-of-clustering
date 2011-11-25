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
	unsigned int j;

	props->solutions = (Solution*)malloc( props->popSize * sizeof(Solution) );

	checkAlloc( props->solutions )
		return GetLastErrorCode();
	}
	// clear all densities
	for ( i = 0; i < props->popSize; i++ ) {
		props->solutions[ i].densities = 0;
		for ( j = 0; j < props->dataStore->info.numEntries; i++ ) {
			props->solutions[ i].recordMembership[j];
		}
	}

	props->blocksPerEntries = props->dataStore->info.numEntries / threadsPerBlock;
	while ( props->blocksPerEntries * threadsPerBlock <  props->dataStore->info.numEntries ) {
		props->blocksPerEntries++;
	}

	return errOk;
}
//----------------------------------------------
void MembershipAndDensityKernel( LoopContext loop ) {
	// for this record in this specific solution (population member)
	unsigned int solution = loop.blockIdx.y;
	unsigned int record = loop.blockIdx.x + loop.threadIdx.x;
	unsigned int i;
	EvolutionProps * props = (EvolutionProps*)loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int clusterPos = 0;
	unsigned int clusterSize = (*thisMember).clusters[ clusterPos];
	float currentDistance = 0;
	float smallestDistance = props->dataStore->distances[ DistanceVIdx( record, (*thisMember).medoids[ 0] )];
	(*thisSolution).recordMembership[ record] = clusterPos;
	clusterSize--;

	// find closest medoid to this record
	for ( i = 1; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( clusterSize <= 0 ) {
			clusterPos++;
			clusterSize = (*thisMember).clusters[ clusterPos];
		}

		currentDistance = props->dataStore->distances[ DistanceVIdx( record, (*thisMember).medoids[ i] )];
		 if ( currentDistance < smallestDistance ) {
			 smallestDistance = currentDistance;
			 (*thisSolution).clusterMembership[ record] = clusterPos;
		 }
	}
	
	(*thisSolution).densities += smallestDistance;
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
	EvolutionProps * props = (EvolutionProps*)loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int i;
	unsigned int memberOf = (*thisSolution).recordMembership[ record];

	(*thisSolution).connectivity = 0;
	for ( i = 0; i < (*thisMember).attr.numNeighbours; i++ ) {
		if ( memberOf == (*thisSolution).recordMembership[ props->dataStore->neighbours[ record * kMaxNeighbours + i]] ) {
			(*thisSolution).connectivity += 1.0 / (float)(*thisMember).attr.numNeighbours;
		}
	}
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

	err = RunLoop( connectivityLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

void DisconnectivityKernel( LoopContext loop ) {
	// For each medoid in solution
	// check distance to medoids from different clusters
	// for each such pair, find medoids with smaller distance
	// to any of the two, and count the ones found
	// -- more it founds, more the bigger "distance" is betwen clusters
	unsigned int solution = loop.blockIdx.x;
	unsigned int medoid = loop.threadIdx.x;
	EvolutionProps * props = (EvolutionProps*)loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int counts;
	unsigned int comparisions;
	float currDistance;
	unsigned int i, j;

	comparisions= 0;
	counts = 0;
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( i == medoid || (*thisSolution).clusterMembership[ i] == (*thisSolution).clusterMembership[ medoid] ) {
			// if same medoid or same cluster - skip
			continue;
		}

		comparisions++;
		currDistance = props->dataStore->distances[ DistanceVIdx( (*thisMember).medoids[ medoid], (*thisMember).medoids[ i] )];
		for ( j = 0; j < MEDOID_VECTOR_SIZE; j++ ) {
			if ( j == medoid || j == i ) {
				// if one of the pair - skip
				continue;
			}
			
			if ( props->dataStore->distances[ DistanceVIdx( (*thisMember).medoids[ medoid], (*thisMember).medoids[ i] )] <
				currDistance ) {
				counts++;
			}
		}
	}

	if ( medoid == 0 ) {
		// first medoid from the list - do the cleaning
		(*thisSolution).disconnectivity = 0;
	}

	(*thisSolution).disconnectivity += ((float)counts) / (float)comparisions;
}
//----------------------------------------------

ErrorCode Disconnectivity( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition disconnectivityLoop;

	if ( props == NULL || props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	// <<< populationSize, medoidsVectorSize >>>
	disconnectivityLoop.gridSize.x = props->popSize;
	disconnectivityLoop.gridSize.y = 1;
	disconnectivityLoop.gridSize.z = 1;
	disconnectivityLoop.blockSize.x = MEDOID_VECTOR_SIZE;
	disconnectivityLoop.blockSize.y = 1;
	disconnectivityLoop.blockSize.z = 1;
	disconnectivityLoop.kernel = DisconnectivityKernel;
	disconnectivityLoop.params = (void*)&props;

	err = RunLoop( disconnectivityLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

void CorrectnessKernel( LoopContext loop ) {
	// For this medoid check all other medoids from this solution.
	// +1 if there are some that repeat
	// +1 if they're both from different clusters
	unsigned int solution = loop.blockIdx.x;
	unsigned int medoid = loop.threadIdx.x;
	EvolutionProps * props = (EvolutionProps*)loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int i;

	if ( medoid == 0 ) {
		(*thisSolution).errors = 0;
	}

	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( i == medoid ) {
			continue;
		}

		if ( (*thisMember).medoids[ i] == (*thisMember).medoids[ medoid] ) {
			(*thisSolution).errors++;
		}

		if ( (*thisSolution).clusterMembership[ i] = (*thisSolution).clusterMembership[ medoid] ) {
			(*thisSolution).errors++;
		}
	}
}
//----------------------------------------------

ErrorCode Correctness( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition correctnessLoop;

	if ( props == NULL || props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	// <<< populationSize, medoidsVectorSize >>>
	correctnessLoop.gridSize.x = props->popSize;
	correctnessLoop.gridSize.y = 1;
	correctnessLoop.gridSize.z = 1;
	correctnessLoop.blockSize.x = MEDOID_VECTOR_SIZE;
	correctnessLoop.blockSize.y = 1;
	correctnessLoop.blockSize.z = 1;
	correctnessLoop.kernel = DisconnectivityKernel;
	correctnessLoop.params = (void*)&props;

	err = RunLoop( correctnessLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------
//----------------------------------------------