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
ErrorCode ConfigureAlgorithms( EvolutionProps * props );

ErrorCode MembershipAndDensity( EvolutionProps * props );

ErrorCode Connectivity( EvolutionProps * props );

ErrorCode Disconnectivity( EvolutionProps * props );

ErrorCode Correctness( EvolutionProps * props );

ErrorCode DominanceCount( EvolutionProps * props );

ErrorCode FrontDensity( FrontDensities * frontProps );

ErrorCode Crossing( BreedingTable * breedingProps );

float SolutionResult( Solution * solution, char objective );

//=============================================

/**
 * Just some defaults.
 */
ErrorCode GenerateDefaultProps( EvolutionProps * props ) {
//TODO: implement
	props;
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
			props->population[ loop.threadIdx.x].clusters[ i] = (char)proposal;
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
	unsigned int i, j, k;
	char * solutionsSelected = (char*)malloc( props->popSize * sizeof(char) );
	unsigned int solutionsLeft = 0;
	unsigned int currFront = 0;
	unsigned int currFrontSize = 0;
	unsigned int * solutionFronts = (unsigned int*)malloc( props->popSize * ( props->popSize + 1 ) * sizeof(unsigned int) );
	FrontDensities frontDensitiesProps;
	char * thisFrontSelection = (char*)malloc( props->popSize * sizeof(char) );
	unsigned int smallest = 0;
	BreedingTable breedingData;
	breedingData.table = NULL;
	//float * currFrontDensities = NULL;

	err = ConfigureAlgorithms( props );
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
		err = Connectivity( props );
		if ( err != errOk ) {
			break;
		}

		// disconnectivity
		err = Disconnectivity( props );
		if ( err != errOk ) {
			break;
		}

		// correctness
		err = Correctness( props );
		if ( err != errOk ) {
			break;
		}

		// sorting
		err = Correctness( props );
		if ( err != errOk ) {
			break;
		}

		// dominance count
		err = DominanceCount( props );
		if ( err != errOk ) {
			break;
		}

		// setup fronts
		solutionsLeft = props->popSize / 2; // no need to sort all
		currFront = 0;

		for ( j = 0; j < props->popSize; j++ ) {
			solutionsSelected[ j] = 0; // false
		}

		// group fronts
		while ( solutionsLeft > 0 ) {
			currFrontSize = 0;
			// select solutions for current front - where domination count is 0
			for ( j = 0; j < props->popSize && solutionsLeft > 0; j++ ) {				
				if ( !solutionsSelected[ j] && props->dominanceCounts[ j] == 0 ) {
					solutionFronts[ currFront * (props->popSize + 1) + (++currFrontSize)] = j;
					solutionsSelected[ j] = 1; // true
					solutionsLeft--;
				}
			}
			solutionFronts[ currFront * props->popSize + 0] = currFrontSize;
						
			if ( solutionsLeft > 0 ) {
				// for each solution dominated by solution from this front - reduce domination count
				for ( j = 0; j < currFrontSize; j++ ) {
					for ( k = 0; k < props->popSize; k++ ) {
						if ( props->dominanceMatrix[ solutionFronts[ currFront * ( props->popSize + 1 ) + j + 1] * props->popSize + k] ) {
							props->dominanceCounts[ k] -= 1;
						}
					}
				}
			}

			// now for next front
			currFront++;
		}

		// selection
		solutionsLeft = props->popSize / 2; // select half size of population
		for ( j = 0; j < props->popSize; j++ ) {
			solutionsSelected[ j] = 0; //false;
		}

		currFront = 0;
		while ( solutionsLeft > 0 ) {
			// if we need more than the current front can offer
			if ( solutionsLeft >= solutionFronts[ currFront * ( props->popSize + 1 ) + 0] ) {
				for ( j = 0; j < solutionFronts[ currFront * ( props->popSize + 1 ) + 0]; j++ ) {
					solutionsSelected[ solutionFronts[ currFront * ( props->popSize + 1 ) + j + 1]] = 1; //true;
					solutionsLeft--;
				}
			} else {
				// this front has more than we need
				currFrontSize = solutionFronts[ currFront * ( props->popSize + 1 ) + 0];

				// Calculate densities for solutions in this front
				frontDensitiesProps.props = props;
				frontDensitiesProps.front = &solutionFronts[ currFront * ( props->popSize + 1 ) + 1];
				frontDensitiesProps.frontSize = currFrontSize;
				frontDensitiesProps.densities = NULL;
				err = FrontDensity( &frontDensitiesProps );
				if ( err != errOk ) {
					break;
				}

				// Select first selectionLeft solutions and find the smallest one (bug density)
				smallest = 0;
				for ( j = 0; j < currFrontSize; j++ ) {
					thisFrontSelection [ j] = ( j < solutionsLeft );
					if ( thisFrontSelection[ j] ) {
						if ( frontDensitiesProps.densities[ j] != -1 && 
							( frontDensitiesProps.densities[ j] < frontDensitiesProps.densities[ smallest] || frontDensitiesProps.densities[ smallest] == -1 ) ) {
							smallest = j;
						}
					}
				} // for j

				// Now for each solution not selected at first, check if it's bigger than the smallest
				// If so, replece it with smallest
				if  ( frontDensitiesProps.densities[ smallest] != -1 ) {
					for ( j = 0; j < currFrontSize; j++ ) {
						if ( thisFrontSelection[ j]) {
							continue;
						}

						if ( frontDensitiesProps.densities[ j] == -1 || 
							frontDensitiesProps.densities[ j] > frontDensitiesProps.densities[ smallest] ) {
							thisFrontSelection[ smallest] = 0; //false;
							thisFrontSelection[ j] = 1; //true;
							smallest = j;
							for ( k = 0; k < j; k++ ) {
								if ( thisFrontSelection[ k] ) {
									if ( frontDensitiesProps.densities[ k] != -1 && 
										( frontDensitiesProps.densities[ k] < frontDensitiesProps.densities[ smallest] || 
										frontDensitiesProps.densities[ smallest] == -1 ) ) {
										smallest = k;
									}
								}
							} // for k
						}
					} // for j
				}

				// now mark solutions in main selection table
				for ( j = 0; j < currFrontSize; j++ ) {
					if ( thisFrontSelection[ j] ) {
						solutionsSelected[ solutionFronts[ currFront * ( props->popSize + 1 ) + j + 1]] = 1; //true;
						solutionsLeft--;
					}
				}// for j
			}

			currFront++;
		} // while

		// crossing		
		// breedingTable[ parent1, parent2, child, mutation probability]
		{
			unsigned int currParent1 = 0;
			unsigned int currParent2 = 0;
			unsigned int currChild = 0;

			if ( breedingData.table == NULL ) {
				breedingData.table = (BreedDescriptor*)malloc( props->popSize / 2 * sizeof(BreedDescriptor) );
			}

			srand( (unsigned int)time( 0 ) );
			// generate breeding Table
			for ( j = 0; j < props->popSize; j++ ) {
				if ( solutionsSelected[ j] ) {
					// place for parent
					if ( currParent1 <= currParent2 ) {
						// place taken by first parent
						breedingData.table[ currParent1++].parent1 = j;
						breedingData.table[ currParent1++].parent2 = j;
					} else {
						breedingData.table[ currParent2++].parent2 = j;
						breedingData.table[ currParent2++].parent1 = j;
					}
				} else {
					// place for child
					breedingData.table[ currChild].child = j;
					// mutation probability
					breedingData.table[ currChild++].factor = rand() % 100;
				}
			} // for j

			
		}

		// launch crossing
		breedingData.props = props;
		err = Crossing( &breedingData );
		if ( err != errOk ) {
			break;
		}
	}

	return err;
}
//----------------------------------------------

ErrorCode ConfigureAlgorithms( EvolutionProps * props ) {
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

	props->dominanceMatrix = (char*)malloc( props->popSize * props->popSize * sizeof(char) );
	checkAlloc( props->dominanceMatrix )
		return GetLastErrorCode();
	}

	props->dominanceCounts = (unsigned int*)malloc( props->popSize * sizeof(unsigned int) );

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
			 (*thisMember).clusterMembership[ record] = clusterPos;
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
	densityLoop.params = (void*)props;

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
			(*thisSolution).connectivity += (float)1.0 / (float)(*thisMember).attr.numNeighbours;
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
	connectivityLoop.params = (void*)props;

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
		if ( i == medoid || (*thisMember).clusterMembership[ i] == (*thisMember).clusterMembership[ medoid] ) {
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
	disconnectivityLoop.params = (void*)props;

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

		if ( (*thisMember).clusterMembership[ i] == (*thisMember).clusterMembership[ medoid] ) {
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
	correctnessLoop.kernel = CorrectnessKernel;
	correctnessLoop.params = (void*)props;

	err = RunLoop( correctnessLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------
void SortingKernel( LoopContext loop ) {
	// true if this solution (blockIdx.x) dominates the other one (threadIdx.x)
	char currDominance = 0;
	unsigned int i;
	char hasBetter = 0, hasWorse = 0;
	EvolutionProps * props = (EvolutionProps*)loop.params;
	unsigned int me = loop.blockIdx.x;
	unsigned int he = loop.threadIdx.x;

	for ( i = 0; i < OBJECTIVES; i++ ) {
		if ( hasWorse && hasBetter ) {
			// it's already equal so no need for more tests
			break;
		}
		switch ( i ) {
			case 0: {// Density
				// smaller better
				if ( props->solutions[ me].densities == props->solutions[ he].densities ) {
					continue;
				} else if ( props->solutions[ me].densities < props->solutions[ he].densities ) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
			case 3: {// Correctnes
				// smaller better
				if ( props->solutions[ me].errors == props->solutions[ he].errors ) {
					continue;
				} else if ( props->solutions[ me].errors < props->solutions[ he].errors ) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
			case 1: {// Connectivity
				// bigger better
				if ( props->solutions[ me].connectivity == props->solutions[ he].connectivity ) {
					continue;
				} else if ( props->solutions[ me].connectivity > props->solutions[ he].connectivity ) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
			case 2: { // Disconnectivity
				// bigger better
				if ( props->solutions[ me].disconnectivity == props->solutions[ he].disconnectivity ) {
					continue;
				} else if ( props->solutions[ me].disconnectivity > props->solutions[ he].disconnectivity ) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
		}
	}

	if ( hasBetter && !hasWorse ) {
		currDominance = 1;
	}

	props->dominanceMatrix[ me * props->popSize + he] = currDominance;
}
//----------------------------------------------
ErrorCode Sorting( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition sortingLoop;

	if ( props == NULL || props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	// <<< populationSize, populationSize >>>
	sortingLoop.gridSize.x = props->popSize;
	sortingLoop.gridSize.y = 1;
	sortingLoop.gridSize.z = 1;
	sortingLoop.blockSize.x = props->popSize;
	sortingLoop.blockSize.y = 1;
	sortingLoop.blockSize.z = 1;
	sortingLoop.kernel = SortingKernel;
	sortingLoop.params = (void*)props;

	err = RunLoop( sortingLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

void DominanceCountKernel( LoopContext loop ) {
	// Counts how many other solutions dominate over this one
	EvolutionProps * props = (EvolutionProps*)loop.params;
	unsigned int i;
	

	props->dominanceCounts[ loop.threadIdx.x] = 0;
	for ( i = 0; i < props->popSize; i++ ) {
		if ( props->dominanceMatrix[ i * props->popSize + loop.threadIdx.x] ) {
			// i dominates over threadIdx.x
			props->dominanceCounts[ loop.threadIdx.x]++ ;
		}
	}
}
//----------------------------------------------

ErrorCode DominanceCount( EvolutionProps * props ) {
	ErrorCode err = errOk;
	LoopDefinition dominanceLoop;

	if ( props == NULL || props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	// <<< 1, poulationSize >>>
	dominanceLoop.gridSize.x = 1;
	dominanceLoop.gridSize.y = 1;
	dominanceLoop.gridSize.z = 1;
	dominanceLoop.blockSize.x = props->popSize;
	dominanceLoop.blockSize.y = 1;
	dominanceLoop.blockSize.z = 1;
	dominanceLoop.kernel = DominanceCountKernel;
	dominanceLoop.params = (void*)props;

	err = RunLoop( dominanceLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

void FrontDensityKernel( LoopContext loop ) {
	FrontDensities * frontProps = (FrontDensities*)loop.params;

	char lesserFound = 0; //false;
	float lesserResult = 0;
	char biggerFound = 0; //false;
	float biggerResult = 0;

	float thisResult = SolutionResult( &frontProps->props->solutions[ frontProps->front[ loop.blockIdx.x]], (char)loop.threadIdx.x );
	float currResult;
	unsigned int i;

	// find the smallest score bigger than this
	// and biggest score smaller than this
	for ( i = 0; i < frontProps->frontSize; i++ ) {
		if ( loop.blockIdx.x == i ) {
			// skip if same
			continue;
		}

		currResult = SolutionResult( &frontProps->props->solutions[ frontProps->front[ i]], (char)loop.threadIdx.x );
		// check if lesser
		if ( thisResult > currResult ) {
			if ( !lesserFound ) {
				//lesser = i;
				lesserFound = 1; //true;
				lesserResult = currResult;
			} else {
				if ( lesserResult < currResult ) {
					//lesser = i;
					lesserFound = 1; //true;
					lesserResult = currResult;
				}
			}
		}

		// check if bigger
		if ( thisResult < currResult ) {
			if ( !biggerFound ) {
				//bigger = i;
				biggerFound = 1; //true;
				biggerResult = currResult;
			} else {
				if ( biggerResult > currResult ) {
					//bigger = i;
					biggerResult = currResult;
				}
			}
		}
	} // for each solution in this front

	// if first objective of this solution, clean
	if ( loop.threadIdx.x == 0 ) {
		frontProps->densities[ loop.blockIdx.x] = 0;
	}

	// is this an edge solution ?
	if ( !lesserFound || !biggerFound || frontProps->densities[ loop.blockIdx.x] == -1) {
		frontProps->densities[ loop.blockIdx.x] = -1;
	} else {
		frontProps->densities[ loop.blockIdx.x] += biggerResult - lesserResult;
	}
}
//----------------------------------------------

ErrorCode FrontDensity( FrontDensities * frontProps ) {
	ErrorCode err = errOk;
	LoopDefinition dominanceLoop;

	if ( frontProps == NULL || frontProps->props == NULL || frontProps->props->solutions == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	// <<< numSolutions, kryterions >>>
	dominanceLoop.gridSize.x = frontProps->props->popSize;
	dominanceLoop.gridSize.y = 1;
	dominanceLoop.gridSize.z = 1;
	dominanceLoop.blockSize.x = OBJECTIVES;
	dominanceLoop.blockSize.y = 1;
	dominanceLoop.blockSize.z = 1;
	dominanceLoop.kernel = FrontDensityKernel;
	dominanceLoop.params = (void*)frontProps;

	err = RunLoop( dominanceLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------

float SolutionResult( Solution * solution, char objective ) {
	switch ( objective ) {
		case 0: { // Density
			return solution->densities;
		} break;
		case 1: { // Connectivity
			return solution->connectivity;
		} break;
		case 2: { // Disconnectivity
			return solution->disconnectivity;
		} break;
		case 3: { // Correctness
			return solution->errors;
		} break;
		default: {
			return 0.0;
		}
	}
}
//----------------------------------------------

void CrossingKernel( LoopContext loop ) {
	BreedingTable * breedingProps = (BreedingTable*)loop.params;
	unsigned int i = 0, j = 0;
	unsigned int stepSize = MEDOID_VECTOR_SIZE / CROS_FACTOR;
	unsigned int thisParent1, thisParent2, thisChild;
	char mark = 0;
	unsigned char howMany;
	unsigned int currCluster;
	unsigned int clusterToWrite;
	unsigned int currMembership;
	unsigned int membersCount;

	thisParent1 = breedingProps->table[ loop.threadIdx.x].parent1;
	thisParent2 = breedingProps->table[ loop.threadIdx.x].parent2;
	thisChild = breedingProps->table[ loop.threadIdx.x].child;
	// generate cross template
	if ( loop.threadIdx.x == 0 ) {
		for ( i = 0, j = 0; i < MEDOID_VECTOR_SIZE; i++, j++ ) {
			if ( j >= stepSize ) {
				mark = ( mark ? 0 : 1 );
			}
			breedingProps->crossTemplate[ i] = mark;
		}
	}

	// start crossing
	/*
	 Exchange both medoids and their membership
	 */
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( breedingProps->crossTemplate[ i] ) {
			//clusterMembership
			breedingProps->props->population[ thisChild].medoids[ i] = breedingProps->props->population[ thisParent1].medoids[ i];
			breedingProps->props->population[ thisChild].clusterMembership[ i] =
				breedingProps->props->population[ thisParent1].clusterMembership[ i];
		} else {
			breedingProps->props->population[ thisChild].medoids[ i] = breedingProps->props->population[ thisParent2].medoids[ i];
			breedingProps->props->population[ thisChild].clusterMembership[ i] =
				breedingProps->props->population[ thisParent2].clusterMembership[ i];
		}
	}

	// copy attributes from the first parent
	breedingProps->props->population[ thisChild].attr = breedingProps->props->population[ thisParent1].attr;

	// mutation
	// 1) attributes
	if ( breedingProps->table[ loop.threadIdx.x].factor > 50 ) {
		if ( breedingProps->table[ loop.threadIdx.x].factor > 70 ) {
			if ( breedingProps->table[ loop.threadIdx.x].factor > 90 ) {
				// both
				breedingProps->props->population[ thisChild].attr.clusterMaxSize = rand() % MAX_CLUSTER_SIZE;
				breedingProps->props->population[ thisChild].attr.numNeighbours = rand() % MAX_NEIGHBOURS;
			} else {
				// neighbours
				breedingProps->props->population[ thisChild].attr.numNeighbours = rand() % MAX_NEIGHBOURS;
			}
		} else {
			// max cluster size
			breedingProps->props->population[ thisChild].attr.clusterMaxSize = rand() % MAX_CLUSTER_SIZE;
		}
	}

	// cleanup cluster membership and sizes
	currCluster = 0;
	clusterToWrite = 1;
	membersCount = 0;
	currMembership = breedingProps->props->population[ thisChild].clusterMembership[ 0];

	/*
	 Till membership doesn't change, and size of a cluster is in a norm count a medoid as a member of this cluster.
	 If membership changes or maximum size of a cluster is meet, make a new cluster+
	 */
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( breedingProps->props->population[ thisChild].clusterMembership[ i] == currMembership && 
			breedingProps->props->population[ thisChild].attr.clusterMaxSize ) {
			membersCount++;
		} else {
			breedingProps->props->population[ thisChild].clusters[ currCluster++] = membersCount;
			membersCount = 1;
			clusterToWrite++;
			breedingProps->props->population[ thisChild].clusterMembership[ i];
		}
		breedingProps->props->population[ thisChild].clusterMembership[ i] = clusterToWrite;
	}

	// mutation
	// 2) medoids
	howMany = 0;
	if ( breedingProps->table[ loop.threadIdx.x].factor > 20 ) {
		if ( breedingProps->table[ loop.threadIdx.x].factor > 50 ) {
			if ( breedingProps->table[ loop.threadIdx.x].factor > 80 ) {
				if ( breedingProps->table[ loop.threadIdx.x].factor > 95 ) {
					// 4
					howMany = 4;
				} else {
					// 3
					howMany = 3;
				}
			} else {
				// 2
				howMany = 2;
			}
		} else {
			// 1
			howMany = 1;
		}
	}
	// generate new random medoids
	for ( i = 0; i < howMany; i++ ) {
		breedingProps->props->population[ thisChild].medoids[ rand() % MEDOID_VECTOR_SIZE] =
			rand() % breedingProps->props->dataStore->info.numEntries;
	}


}
//----------------------------------------------

ErrorCode Crossing( BreedingTable * breedingProps ) {
	ErrorCode err = errOk;
	LoopDefinition crossingLoop;

	if ( breedingProps == NULL || breedingProps->props == NULL || breedingProps->table == NULL ) {
		return SetLastErrorCode( errWrongParameter );
	}

	// <<< 1, popSize/2 >>>
	crossingLoop.gridSize.x = 1;
	crossingLoop.gridSize.y = 1;
	crossingLoop.gridSize.z = 1;
	crossingLoop.blockSize.x = breedingProps->props->popSize/2;
	crossingLoop.blockSize.y = 1;
	crossingLoop.blockSize.z = 1;
	crossingLoop.kernel = CrossingKernel;
	crossingLoop.params = (void*)breedingProps;

	err = RunLoop( crossingLoop );

	if ( err != errOk ) {
		reportError( err, "Run loop returned with error%s", "" );
	}
	return err;
}
//----------------------------------------------
ErrorCode calculateBDI( void ) {
}
//----------------------------------------------

ErrorCode calculateDI( void ) {
}
//----------------------------------------------

ErrorCode calculateRand( void ) {
}
//----------------------------------------------
//----------------------------------------------