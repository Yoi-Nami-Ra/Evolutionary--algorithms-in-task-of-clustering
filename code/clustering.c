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

#pragma mark - Globals
//==============================================
//== Globals

#define CROS_FACTOR 3
#define OBJECTIVES 4

const char threadsPerBlock = 50;

#pragma mark - Function Prototypes
//==============================================
//== Functions

/**
 * Generates random population to start with.
 * @param props	- [in, out] For parameters and results of the operation.
 */
ErrorCode GenerateRandomPopulation(EvolutionProps * props);

/**
 * Main algorithms.
 */
ErrorCode RunAlgorithms(EvolutionProps * props);

/**
 * Prepare enviroment to run algorithms.
 */
ErrorCode ConfigureAlgorithms(EvolutionProps * props);

ErrorCode MembershipAndDensity(EvolutionProps * props);

ErrorCode Connectivity(EvolutionProps * props);

ErrorCode Disconnectivity(EvolutionProps * props);

ErrorCode Correctness(EvolutionProps * props);

ErrorCode DominanceCount(EvolutionProps * props);

ErrorCode Sorting(EvolutionProps * props);

ErrorCode FrontDensity(FrontDensities * frontProps);

ErrorCode Crossing(BreedingTable * breedingProps);

float SolutionResult(Solution * solution, char objective);

ErrorCode calculateBDI(EvolutionProps * props);

ErrorCode calculateDI(EvolutionProps * props);

ErrorCode calculateRand(EvolutionProps * props);

float Distance(EvolutionProps * props, unsigned int a, unsigned int b);

#pragma mark - -
//=============================================

/**
 * Just some defaults.
 */
ErrorCode DefaultProps(EvolutionProps * props, DataStore * dataStore) {
    props->blocksPerEntries = 0;
    props->crosFactor = 0;
    props->dataStore = dataStore;
    props->dominanceCounts = NULL;
    props->dominanceMatrix = NULL;
    props->evoSteps = 0;
    props->popSize = 0;
    props->medoidsVectorSize = 0;
    props->maxNeighbours = 0;
    props->maxClusterSize = 0;
    props->population = NULL;
    props->solutions = NULL;
    
    return ConfigureAlgorithms( props );
}
//----------------------------------------------

ErrorCode RunClustering(EvolutionProps * props) {
	ErrorCode err = errOk;

	err = GenerateRandomPopulation(props);
	if (err != errOk) {
		reportError( err, "Failed to generate random population%s", "");
		return err;
	}

	err = RunAlgorithms(props);

	if (err != errOk) {
		reportError( err, "Failed with algorithms%s", "");
	}

	return err;
}
//----------------------------------------------

float Distance(EvolutionProps * props, unsigned int a, unsigned int b) {
	if (a == b)
		return 0.0;
	return props->dataStore->distances[DistanceVIdx(a, b)];
}
//----------------------------------------------

void GenerateRandomPopulationKernel(LoopContext loop);
void GenerateRandomPopulationKernel(LoopContext loop) {
	EvolutionProps * props = (EvolutionProps*) loop.params;
	char proposalOk;
	unsigned int proposal;
	unsigned int clustersSum = props->medoidsVectorSize;
	int i, j;

	// local max size of a cluster
	props->population[loop.threadIdx.x].attr.clusterMaxSize = rand()
			% props->maxClusterSize + 1;
	// local max list of neighbours
	props->population[loop.threadIdx.x].attr.numNeighbours = rand()
			% props->maxNeighbours + 1;

	// for each medoid in vector
	for (i = 0; i < props->medoidsVectorSize; i++) {
		do {
			proposalOk = 1;
			// chose random data entry
			proposal = rand() % props->dataStore->info.numEntries;
			// and check if we didn't that already
			for (j = 0; j < i; j++) {
				if (proposal
						== props->population[loop.threadIdx.x].medoids[j]) {
					proposalOk = 0;
				}
			}
		} while (!proposalOk);
		// save it
		props->population[loop.threadIdx.x].medoids[i] = proposal;

		// now clusters
		// if we didn't fill it up yet
		if (clustersSum > 0) {
			proposal = rand()
					% props->population[loop.threadIdx.x].attr.clusterMaxSize
					+ 1;
			if (clustersSum < proposal) {
				proposal += clustersSum - proposal;
			}
			clustersSum -= proposal;
			props->population[loop.threadIdx.x].clusters[i] = (char) proposal;
		} else {
			props->population[loop.threadIdx.x].clusters[i] = 0;
		}
	} // i
}
//----------------------------------------------

ErrorCode GenerateRandomPopulation(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition popGenLoop;
    unsigned int i;

	logDebug( " GenerateRandomPopulation%s", "" );

	if (props == NULL) {
		reportError( errWrongParameter,
				"Evolution Propertis served with NULL pointer.%s", "" );
	}

    for ( i = 0; i < props->popSize; i++ ) {
        props->population[ i].medoids =
            props->population[ i].clusters = NULL;
    }

	srand((unsigned int) time(NULL));

	popGenLoop.blockSize.x = props->popSize;
	popGenLoop.blockSize.y = popGenLoop.gridSize.x = popGenLoop.gridSize.y =
			popGenLoop.blockSize.z = popGenLoop.gridSize.z = 1;
	popGenLoop.kernel = GenerateRandomPopulationKernel;
	popGenLoop.params = (void*) props;

	err = RunLoop(popGenLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

ErrorCode RunAlgorithms(EvolutionProps * props) {
	ErrorCode err = errOk;
	unsigned int i, j, k;
	char * solutionsSelected = (char*) malloc(props->popSize * sizeof(char));
	unsigned int solutionsLeft = 0;
	unsigned int currFront = 0;
	unsigned int currFrontSize = 0;
	unsigned int * solutionFronts = (unsigned int*) malloc(
			props->popSize * (props->popSize + 1) * sizeof(unsigned int));
	FrontDensities frontDensitiesProps;
	char * thisFrontSelection = (char*) malloc(props->popSize * sizeof(char));
	unsigned int smallest = 0;
	BreedingTable breedingData;
	float bestBDI = 100;
	float bestDI = 100;
	float bestRAND = 10;

	breedingData.table = NULL;

	err = ConfigureAlgorithms(props);
	if (err != errOk) {
		return err;
	}

    printf(" >>Evolution<<\n population:%u, steps:%u\n medoids:%u, maxClusterSize:%u, maxNeighbour:%u\n", props->popSize, props->evoSteps, props->medoidsVectorSize, props->maxClusterSize, props->maxNeighbours );
    // Starting evolution loop
	for (i = 0; i < props->evoSteps; i++) {
		logDebug( " Evolution Step:%u", i );
		// mambership and density
		err = MembershipAndDensity(props);
		if (err != errOk) {
			break;
		}

		// connectivity
		err = Connectivity(props);
		if (err != errOk) {
			break;
		}

		// disconnectivity
		err = Disconnectivity(props);
		if (err != errOk) {
			break;
		}

		// correctness
		err = Correctness(props);
		if (err != errOk) {
			break;
		}

		// sorting
		err = Sorting(props);
		if (err != errOk) {
			break;
		}

		// dominance count
		err = DominanceCount(props);
		if (err != errOk) {
			break;
		}

		// setup fronts
		solutionsLeft = props->popSize;
		currFront = 0;

		for (j = 0; j < props->popSize; j++) {
			solutionsSelected[j] = 0; // false
		}

		// group fronts
		while (solutionsLeft > 0) {
			currFrontSize = 0;
			// select solutions for current front - where domination count is 0
			for (j = 0; j < props->popSize && solutionsLeft > 0; j++) {
				if (!solutionsSelected[ j] && props->dominanceCounts[ j] == 0) {
					solutionFronts[ currFront * ( props->popSize + 1 )
							+ ( ++currFrontSize )] = j;
					solutionsSelected[ j] = 1; // true
					solutionsLeft--;
				}
			}
			solutionFronts[currFront * (props->popSize + 1) + 0] =
					currFrontSize;

			if (solutionsLeft > 0) {
				// for each solution dominated by solution from this front - reduce domination count
				for (j = 0; j < currFrontSize; j++) {
					for (k = 0; k < props->popSize; k++) {
						if ( props->dominanceMatrix[solutionFronts[currFront * (props->popSize + 1) + j + 1] * props->popSize + k] && ( props->dominanceCounts[ k] > 0 ) ) {
							props->dominanceCounts[ k] -= 1;
						}
					}
				}
			}

			// now for next front
			currFront++;
		} // while

		// selection
		solutionsLeft = props->popSize / 2; // select half size of population
		for (j = 0; j < props->popSize; j++) {
			solutionsSelected[j] = 0; //false;
		}

		currFront = 0;
		while (solutionsLeft > 0) {
			// if we need more than the current front can offer
			if (solutionsLeft
					>= solutionFronts[ currFront * (props->popSize + 1) + 0]) {
				for (j = 0;
						j < solutionFronts[ currFront * (props->popSize + 1) + 0];
						j++) {
					solutionsSelected[ solutionFronts[ currFront
							* (props->popSize + 1) + j + 1]] = 1; //true;
					solutionsLeft--;
				}
			} else {
				// this front has more than we need
				currFrontSize = solutionFronts[currFront * (props->popSize + 1)
						+ 0];

				// Calculate densities for solutions in this front
				frontDensitiesProps.props = props;
				frontDensitiesProps.front = &solutionFronts[currFront
						* (props->popSize + 1) + 1];
				frontDensitiesProps.frontSize = currFrontSize;
				frontDensitiesProps.densities = NULL;
				err = FrontDensity(&frontDensitiesProps);
				if (err != errOk) {
					break;
				}

				// Select first selectionLeft solutions and find the smallest one (bug density)
				smallest = 0;
				for (j = 0; j < currFrontSize; j++) {
					thisFrontSelection[j] = (j < solutionsLeft);
					if (thisFrontSelection[j]) {
						if (frontDensitiesProps.densities[j] != -1
								&& (frontDensitiesProps.densities[j]
										< frontDensitiesProps.densities[smallest]
										|| frontDensitiesProps.densities[smallest]
												== -1)) {
							smallest = j;
						}
					}
				} // for j

				// Now for each solution not selected at first, check if it's bigger than the smallest
				// If so, replece it with smallest
				if (frontDensitiesProps.densities[smallest] != -1) {
					for (j = 0; j < currFrontSize; j++) {
						if (thisFrontSelection[j]) {
							continue;
						}

						if (frontDensitiesProps.densities[j] == -1
								|| frontDensitiesProps.densities[j]
										> frontDensitiesProps.densities[smallest]) {
							thisFrontSelection[smallest] = 0; //false;
							thisFrontSelection[j] = 1; //true;
							smallest = j;
							for (k = 0; k < j; k++) {
								if (thisFrontSelection[k]) {
									if (frontDensitiesProps.densities[k] != -1
											&& (frontDensitiesProps.densities[k]
													< frontDensitiesProps.densities[smallest]
													|| frontDensitiesProps.densities[smallest]
															== -1)) {
										smallest = k;
									}
								}
							} // for k
						}
					} // for j
				}

				// now mark solutions in main selection table
				for (j = 0; j < currFrontSize; j++) {
					if (thisFrontSelection[j]) {
						solutionsSelected[solutionFronts[currFront
								* (props->popSize + 1) + j + 1]] = 1; //true;
						solutionsLeft--;
					}
				} // for j
			}

			currFront++;
		} // while

		// don't bother if it's the last pass
		if (props->evoSteps == i + 1) {
			break;
		}

		// crossing		
		// breedingTable[ parent1, parent2, child, mutation probability]
		{
			unsigned int currParent1 = 0;
			unsigned int currParent2 = 0;
			unsigned int currChild = 0;

			if (breedingData.table == NULL) {
				breedingData.table = (BreedDescriptor*) malloc(
						props->popSize / 2 * sizeof(BreedDescriptor));
			}

			srand((unsigned int) time(0));
			// generate breeding Table
			for (j = 0; j < props->popSize; j++) {
				if (solutionsSelected[ j]) {
					// place for parent
					if (currParent1 <= currParent2) {
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
		err = Crossing(&breedingData);
		if (err != errOk) {
			break;
		}

		// gather results
		calculateBDI(props);
		calculateDI(props);
		calculateRand(props);
		if (i == 0) {
			bestBDI = props->solutions[j].resultBDI;
			bestDI = props->solutions[j].resultDI;
			bestRAND = props->solutions[j].resultRand;
		}
        
        {
            char thisSolution = 0;
            for (j = 0; j < props->popSize; j++) {
                if (bestBDI > props->solutions[j].resultBDI) {
                    bestBDI = props->solutions[j].resultBDI;
                    printf("\n");
                    logMessage( " BDI: %f", bestBDI);
                    thisSolution = 1;
                }
                
                if (bestDI > props->solutions[j].resultDI) {
                    bestDI = props->solutions[j].resultDI;
                    printf("\n");
                    logMessage( " DI: %f", bestDI);
                    thisSolution = 1;
                }
                
                if (bestRAND < props->solutions[j].resultRand) {
                    bestRAND = props->solutions[j].resultRand;
                    printf("\n");
                    logMessage( " RAND: %f", bestRAND);
                    thisSolution = 1;
                }
                
                if (thisSolution) {
                    logMessage("  // density:        %f", props->solutions[j].densities );
                    logMessage("  // connectivity:   %f", props->solutions[j].connectivity );
                    logMessage("  // disconnectivity:%f", props->solutions[j].disconnectivity );
                    logMessage("  // errors:         %f", props->solutions[j].errors );
                    thisSolution = 0;
                }
            }
        }
		
        
        if (!(i % 100)) {
            printf("0\n");
        }
        if (i % 10) {
            printf(".");
        } else {
            printf("/");
        }
        
	} // evolution for

	// gather results
	calculateBDI(props);
	calculateDI(props);
	calculateRand(props);
	logMessage( "\n == Results == <%u>", solutionFronts[0]);
	for (i = 0; i < solutionFronts[0]; i++) {
		if (solutionsSelected[solutionFronts[i + 1]]) {
			logMessage(" = Solution[ %u]:", solutionFronts[ i+1]);
			logMessage("   BDI: %f",
					props->solutions[ solutionFronts[ i+1]].resultBDI);
			logMessage("   DI: %f",
					props->solutions[ solutionFronts[ i+1]].resultDI);
			logMessage("   Rand: %f",
					props->solutions[ solutionFronts[ i+1]].resultRand);
		}
	}
    
    FILE * resDump = fopen("rands.txt", "w");
    if ( resDump ) {
        for ( i = 0; i < props->popSize; i++ ) {
            fprintf(resDump, " %f %s %s\n", props->solutions[ i].resultRand, solutionsSelected[ i]?"<-":"", (props->solutions[ i].resultRand>80)?"=====":"" );
            fprintf(resDump, "  // density:        %f\n", props->solutions[ i].densities );
            fprintf(resDump, "  // connectivity:   %f\n", props->solutions[ i].connectivity );
            fprintf(resDump, "  // disconnectivity:%f\n", props->solutions[ i].disconnectivity );
            fprintf(resDump, "  // errors:         %f\n", props->solutions[ i].errors );
        }
        fclose( resDump );
    }

	return err;
}
//----------------------------------------------

ErrorCode ConfigureAlgorithms(EvolutionProps * props) {
	unsigned int i;
	unsigned int j;

	logDebug(" Configure Algorithms %s", "" );
    
    if ( props->population == NULL ) {
        props->population = (PopMember*)malloc( props->popSize * sizeof(PopMember) );
    } else {
        props->population = (PopMember*)realloc( props->population, props->popSize * sizeof(PopMember) );
    }
    
    if ( props->population[loop.threadIdx.x].medoids == NULL ) {
        props->population[loop.threadIdx.x].medoids = (unsigned int*)malloc( props->medoidsVectorSize * sizeof(unsigned int) );
        props->population[loop.threadIdx.x].clusters = (unsigned int*)malloc( props->medoidsVectorSize * sizeof(unsigned int) );
    } else {
        props->population[loop.threadIdx.x].medoids = (unsigned int*)realloc( props->population[loop.threadIdx.x].medoids, 
                                                                             props->medoidsVectorSize * sizeof(unsigned int) );
        props->population[loop.threadIdx.x].clusters = (unsigned int*)realloc( props->population[loop.threadIdx.x].medoids, 
                                                                              props->medoidsVectorSize * sizeof(unsigned int) );
    }
    
    
    
    if (frontProps->densities == NULL) {
        frontProps->densities = (float*)malloc(frontProps->frontSize * sizeof(float));
    } else {
        frontProps->densities = (float*)realloc(frontProps->densities, frontProps->frontSize * sizeof(float));
    }
    
    props->solutions[i].recordMembership = (unsigned int*) malloc(
                                                                  props->dataStore->info.numEntries * sizeof(unsigned int));

	props->solutions = (Solution*) malloc(props->popSize * sizeof(Solution));

	checkAlloc( props->solutions )
		return GetLastErrorCode();
	}

	// clear all densities
	for (i = 0; i < props->popSize; i++) {
		props->solutions[i].densities = 0;
		
		for (j = 0; j < props->dataStore->info.numEntries; j++) {
			props->solutions[i].recordMembership[j] = 0;
		}
	}

	props->blocksPerEntries = props->dataStore->info.numEntries
			/ threadsPerBlock;
	while (props->blocksPerEntries * threadsPerBlock
			< props->dataStore->info.numEntries) {
		props->blocksPerEntries++;
	}

	props->dominanceMatrix = (char*) malloc(
			props->popSize * props->popSize * sizeof(char));
	checkAlloc( props->dominanceMatrix )
		return GetLastErrorCode();
	}

	props->dominanceCounts = (unsigned int*) malloc(
			props->popSize * sizeof(unsigned int));

	return errOk;
}
//----------------------------------------------

void MembershipAndDensityKernel(LoopContext loop);
void MembershipAndDensityKernel(LoopContext loop) {
	// for this record in this specific solution (population member)
	unsigned int solution = loop.blockIdx.y;
	unsigned int record = loop.blockIdx.x * threadsPerBlock + loop.threadIdx.x;
	unsigned int i;
	EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int clusterPos = 0;
	unsigned int clusterSize = 0; //(*thisMember).clusters[ clusterPos];
	float currentDistance = 0;
	float smallestDistance = 0;
	unsigned int smallestDClusterPos = 0;

	if (record >= props->dataStore->info.numEntries) {
		return;
	}

	smallestDistance = Distance(props, record, (*thisMember).medoids[0]) + 0.1;
	(*thisSolution).recordMembership[record] = clusterPos;

	if (record == 0) {
		// if the first record
		// clean cluster densities
		for (i = 0; i < props->medoidsVectorSize; i++) {
			(*thisSolution).clusterDensities[i] = 0;
		}
	}

	// find closest medoid to this record
	for (i = 0; i < props->medoidsVectorSize; i++, clusterSize--) {
		if (clusterSize <= 0) {
			(i == 0) ? clusterPos = 0 : clusterPos++;
			clusterSize = (*thisMember).clusters[ clusterPos];
		}

		(*thisMember).clusterMembership[ i] = clusterPos + 1;

		currentDistance = Distance(props, record, (*thisMember).medoids[i]);
		if (currentDistance < smallestDistance) {
			smallestDClusterPos = clusterPos;
			smallestDistance = currentDistance;
			(*thisSolution).recordMembership[record] = clusterPos + 1;
		}
	}

	(*thisSolution).numOfClusters = clusterPos + 1;
	(*thisSolution).clusterDensities[smallestDClusterPos] += smallestDistance;
	(*thisSolution).densities += smallestDistance;
}
//----------------------------------------------

ErrorCode MembershipAndDensity(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition densityLoop;

	if (props == NULL || props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == MembershipAndDensity %s", "" );

	//<<( blocksPerEntires, populationSize ), threadsPerBlock>>
	densityLoop.gridSize.x = props->blocksPerEntries;
	densityLoop.gridSize.y = props->popSize;
	densityLoop.gridSize.z = 1;
	densityLoop.blockSize.x = threadsPerBlock;
	densityLoop.blockSize.y = densityLoop.blockSize.z = 1;
	densityLoop.kernel = MembershipAndDensityKernel;
	densityLoop.params = (void*) props;

	err = RunLoop(densityLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void ConnectivityKernel(LoopContext loop);
void ConnectivityKernel(LoopContext loop) {
	unsigned int solution = loop.blockIdx.y;
	unsigned int record = loop.blockIdx.x * threadsPerBlock + loop.threadIdx.x;
	EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int i;
	unsigned int memberOf;

	if (record >= props->dataStore->info.numEntries)
		return;

	memberOf = (*thisSolution).recordMembership[record];

	if (record == 0) {
		(*thisSolution).connectivity = 0;
	}

	// for each record - how many of its neighbours belong to the same cluster
	// 1.0 means all of them

	for (i = 0; i < (*thisMember).attr.numNeighbours; i++) {
		if (memberOf
				== (*thisSolution).recordMembership[props->dataStore->neighbours[record
						* kMaxNeighbours + i]]) {
			(*thisSolution).connectivity += (float) 1.0
					/ (float) (*thisMember).attr.numNeighbours;
		}
	}
}
//----------------------------------------------
ErrorCode Connectivity(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition connectivityLoop;

	if (props == NULL || props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == Connectivity %s", "" );
	//<<( blocksPerEntires, populationSize ), threadsPerBlock>>
	connectivityLoop.gridSize.x = props->blocksPerEntries;
	connectivityLoop.gridSize.y = props->popSize;
	connectivityLoop.gridSize.z = 1;
	connectivityLoop.blockSize.x = threadsPerBlock;
	connectivityLoop.blockSize.y = connectivityLoop.blockSize.z = 1;
	connectivityLoop.kernel = ConnectivityKernel;
	connectivityLoop.params = (void*) props;

	err = RunLoop(connectivityLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void DisconnectivityKernel(LoopContext loop);
void DisconnectivityKernel(LoopContext loop) {
	// For each medoid in solution
	// check distance to medoids from different clusters
	// for each such pair, find medoids with smaller distance
	// to any of the two, and count the ones found
	// -- more it founds, more the bigger "distance" is betwen clusters
	unsigned int solution = loop.blockIdx.x;
	unsigned int medoid = loop.threadIdx.x;
	EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int counts;
	unsigned int comparisions;
	float currDistance;
	unsigned int i, j;

	// for each two medoids not from the same cluster
	// find how far they are from each other using MND (mutual neighbour distance)
	// MND( A, B ) = NN(A,B) + NN(B,A)
	// NN(x,y) - how many objects there are around x in radius of |xy|

	comparisions = 1; // start with 1 as we compare those two medoids
	counts = 2; // start with 2 for both medoids
	for (i = 0; i < props->medoidsVectorSize; i++) {
		if (i == medoid
				|| (*thisMember).clusterMembership[i]
						== (*thisMember).clusterMembership[medoid]) {
			// if same medoid or same cluster - skip
			continue;
		}

		comparisions++;
		currDistance = Distance(props, (*thisMember).medoids[medoid],
				(*thisMember).medoids[i]);
		for (j = 0; j < props->medoidsVectorSize; j++) {
			if (j == medoid || j == i) {
				// if one of the pair - skip
				continue;
			}

			if (Distance(props, (*thisMember).medoids[medoid],
					(*thisMember).medoids[j]) < currDistance) {
				counts++;
			}

			if (Distance(props, (*thisMember).medoids[i],
					(*thisMember).medoids[j]) < currDistance) {
				counts++;
			}
		}
	}

	if (medoid == 0) {
		// first medoid from the list - do the cleaning
		(*thisSolution).disconnectivity = 0;
	}

	(*thisSolution).disconnectivity += ((float) counts) / (float) comparisions;

	// calculate densities for each cluster
	if (medoid < (*thisSolution).numOfClusters) {
		(*thisSolution).clusterDensities[medoid] /=
				(float) (*thisMember).clusters[medoid];
	}
}
//----------------------------------------------

ErrorCode Disconnectivity(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition disconnectivityLoop;

	if (props == NULL || props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == Disconnectivity %s", "" );
	// <<< populationSize, medoidsVectorSize >>>
	disconnectivityLoop.gridSize.x = props->popSize;
	disconnectivityLoop.gridSize.y = 1;
	disconnectivityLoop.gridSize.z = 1;
	disconnectivityLoop.blockSize.x = props->medoidsVectorSize;
	disconnectivityLoop.blockSize.y = 1;
	disconnectivityLoop.blockSize.z = 1;
	disconnectivityLoop.kernel = DisconnectivityKernel;
	disconnectivityLoop.params = (void*) props;

	err = RunLoop(disconnectivityLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void CorrectnessKernel(LoopContext loop);
void CorrectnessKernel(LoopContext loop) {
	// For this medoid check all other medoids from this solution.
	// +1 if there are some that repeat
	// +1 if they're both from different clusters
	unsigned int solution = loop.blockIdx.x;
	unsigned int medoid = loop.threadIdx.x;
	EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + solution;
	PopMember *thisMember = props->population + solution;
	unsigned int i;

	if (medoid == 0) {
		(*thisSolution).errors = 0;
	}

	for (i = 0; i < props->medoidsVectorSize; i++) {
		if (i == medoid) {
			continue;
		}

		if ((*thisMember).medoids[i] == (*thisMember).medoids[medoid]) {
			(*thisSolution).errors++;
		}
	}
}
//----------------------------------------------

ErrorCode Correctness(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition correctnessLoop;

	if (props == NULL || props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == Correctness %s", "" );
	// <<< populationSize, medoidsVectorSize >>>
	correctnessLoop.gridSize.x = props->popSize;
	correctnessLoop.gridSize.y = 1;
	correctnessLoop.gridSize.z = 1;
	correctnessLoop.blockSize.x = props->medoidsVectorSize;
	correctnessLoop.blockSize.y = 1;
	correctnessLoop.blockSize.z = 1;
	correctnessLoop.kernel = CorrectnessKernel;
	correctnessLoop.params = (void*) props;

	err = RunLoop(correctnessLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------
void SortingKernel(LoopContext loop);
void SortingKernel(LoopContext loop) {
	// true if this solution (blockIdx.x) dominates the other one (threadIdx.x)
	char currDominance = 0;
	unsigned int i;
	char hasBetter = 0, hasWorse = 0;
	EvolutionProps * props = (EvolutionProps*) loop.params;
	unsigned int me = loop.blockIdx.x;
	unsigned int he = loop.threadIdx.x;

	for (i = 0; i < OBJECTIVES; i++) {
		if ((hasWorse && hasBetter) || (me == he)) {
			// it's already equal so no need for more tests
			// and there is no need to compare them self
			break;
		}
		switch (i) {
			case 0: { // Density
					// smaller better
					if (props->solutions[me].densities
							== props->solutions[he].densities) {
						continue;
					} else if (props->solutions[me].densities
							< props->solutions[he].densities) {
						hasBetter = 1;
					} else {
						hasWorse = 1;
					}
				} break;
			case 3: { // Correctnes
				// smaller better
				if (props->solutions[me].errors == props->solutions[he].errors) {
					continue;
				} else if (props->solutions[me].errors
						< props->solutions[he].errors) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
			case 1: { // Connectivity
				// bigger better
				if (props->solutions[me].connectivity
						== props->solutions[he].connectivity) {
					continue;
				} else if (props->solutions[me].connectivity
						> props->solutions[he].connectivity) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
			case 2: { // Disconnectivity
				// bigger better
				if (props->solutions[me].disconnectivity
						== props->solutions[he].disconnectivity) {
					continue;
				} else if (props->solutions[me].disconnectivity
						> props->solutions[he].disconnectivity) {
					hasBetter = 1;
				} else {
					hasWorse = 1;
				}
			} break;
		}
	}

	if (hasBetter && !hasWorse) {
		currDominance = 1;
	}

	props->dominanceMatrix[me * props->popSize + he] = currDominance;
}
//----------------------------------------------

ErrorCode Sorting(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition sortingLoop;

	if (props == NULL || props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == Sorting %s", "" );
	// <<< populationSize, populationSize >>>
	sortingLoop.gridSize.x = props->popSize;
	sortingLoop.gridSize.y = 1;
	sortingLoop.gridSize.z = 1;
	sortingLoop.blockSize.x = props->popSize;
	sortingLoop.blockSize.y = 1;
	sortingLoop.blockSize.z = 1;
	sortingLoop.kernel = SortingKernel;
	sortingLoop.params = (void*) props;

	err = RunLoop(sortingLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void DominanceCountKernel(LoopContext loop);
void DominanceCountKernel(LoopContext loop) {
	// Counts how many other solutions dominate over this one
	EvolutionProps * props = (EvolutionProps*) loop.params;
	unsigned int i;

	props->dominanceCounts[loop.threadIdx.x] = 0;
	for (i = 0; i < props->popSize; i++) {
		if (props->dominanceMatrix[i * props->popSize + loop.threadIdx.x]) {
			// i dominates over threadIdx.x
			props->dominanceCounts[loop.threadIdx.x]++;
		}
	}
}
//----------------------------------------------

ErrorCode DominanceCount(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition dominanceLoop;

	if (props == NULL || props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == DominanceCount %s", "" );
	// <<< 1, poulationSize >>>
	dominanceLoop.gridSize.x = 1;
	dominanceLoop.gridSize.y = 1;
	dominanceLoop.gridSize.z = 1;
	dominanceLoop.blockSize.x = props->popSize;
	dominanceLoop.blockSize.y = 1;
	dominanceLoop.blockSize.z = 1;
	dominanceLoop.kernel = DominanceCountKernel;
	dominanceLoop.params = (void*) props;

	err = RunLoop(dominanceLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void FrontDensityKernel(LoopContext loop);
void FrontDensityKernel(LoopContext loop) {
	FrontDensities * frontProps = (FrontDensities*) loop.params;

	char lesserFound = 0; //false;
	float lesserResult = 0;
	char biggerFound = 0; //false;
	float biggerResult = 0;

	float thisResult = SolutionResult(
			&frontProps->props->solutions[frontProps->front[loop.blockIdx.x]],
			(char) loop.threadIdx.x);
	float currResult;
	unsigned int i;

	// find the smallest score bigger than this
	// and biggest score smaller than this
	for (i = 0; i < frontProps->frontSize; i++) {
		if (loop.blockIdx.x == i) {
			// skip if same
			continue;
		}

		currResult = SolutionResult(
				&frontProps->props->solutions[frontProps->front[i]],
				(char) loop.threadIdx.x);
		// check if lesser
		if (thisResult > currResult) {
			if (!lesserFound) {
				//lesser = i;
				lesserFound = 1; //true;
				lesserResult = currResult;
			} else {
				if (lesserResult < currResult) {
					//lesser = i;
					lesserFound = 1; //true;
					lesserResult = currResult;
				}
			}
		}

		// check if bigger
		if (thisResult < currResult) {
			if (!biggerFound) {
				//bigger = i;
				biggerFound = 1; //true;
				biggerResult = currResult;
			} else {
				if (biggerResult > currResult) {
					//bigger = i;
					biggerResult = currResult;
				}
			}
		}
	} // for each solution in this front

	// if first objective of this solution, clean
	if (loop.threadIdx.x == 0) {
		frontProps->densities[loop.blockIdx.x] = 0;
	}

	// is this an edge solution ?
	if (!lesserFound || !biggerFound
			|| frontProps->densities[loop.blockIdx.x] == -1) {
		frontProps->densities[loop.blockIdx.x] = -1;
	} else {
		frontProps->densities[loop.blockIdx.x] += biggerResult - lesserResult;
	}
}
//----------------------------------------------

ErrorCode FrontDensity(FrontDensities * frontProps) {
	ErrorCode err = errOk;
	LoopDefinition dominanceLoop;

	if (frontProps == NULL || frontProps->props == NULL
			|| frontProps->props->solutions == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == FrontDensity %s", "" );
	// <<< numSolutions, kryterions >>>
	dominanceLoop.gridSize.x = frontProps->props->popSize;
	dominanceLoop.gridSize.y = 1;
	dominanceLoop.gridSize.z = 1;
	dominanceLoop.blockSize.x = OBJECTIVES;
	dominanceLoop.blockSize.y = 1;
	dominanceLoop.blockSize.z = 1;
	dominanceLoop.kernel = FrontDensityKernel;
	dominanceLoop.params = (void*) frontProps;

	err = RunLoop(dominanceLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

float SolutionResult(Solution * solution, char objective) {
	switch (objective) {
	case 0: { // Density
		return solution->densities;
	}
		break;
	case 1: { // Connectivity
		return solution->connectivity;
	}
		break;
	case 2: { // Disconnectivity
		return solution->disconnectivity;
	}
		break;
	case 3: { // Correctness
		return solution->errors;
	}
		break;
	default: {
		return 0.0;
	}
	}
}
//----------------------------------------------

void CrossingKernel(LoopContext loop);
void CrossingKernel(LoopContext loop) {
	BreedingTable * breedingProps = (BreedingTable*) loop.params;
	unsigned int i = 0;
	unsigned int thisParent1, thisParent2, thisChild;
	unsigned char howMany;
	unsigned int currCluster;
	unsigned int clusterToWrite;
	unsigned int currMembership;
	unsigned int membersCount;

	thisParent1 = breedingProps->table[loop.threadIdx.x].parent1;
	thisParent2 = breedingProps->table[loop.threadIdx.x].parent2;
	thisChild = breedingProps->table[loop.threadIdx.x].child;

	// start crossing
	/*
	 Exchange both medoids and their membership
	 */
	howMany = 0;
	currCluster = 0;
	for (i = 0; i < breedingProps->props->medoidsVectorSize; i++) {
		/*
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
		 */
		breedingProps->props->population[thisChild].clusterMembership[i] = 0;
		breedingProps->props->population[thisChild].clusters[i] = 0;
		breedingProps->props->population[thisChild].medoids[i] =
				breedingProps->props->population[thisParent1].medoids[i];
		if (breedingProps->props->population[thisParent1].clusterMembership[i]
				% 2) {
			currCluster =
					breedingProps->props->population[thisChild].clusterMembership[i] =
							breedingProps->props->population[thisParent1].clusterMembership[i];
			howMany++;
		}

		if (!(breedingProps->props->population[thisParent2].clusterMembership[i]
				% 2)) {
			if (breedingProps->props->population[thisChild].clusterMembership[i]
					== 0) {
				breedingProps->props->population[thisChild].medoids[i] =
						breedingProps->props->population[thisParent2].medoids[i];
				currCluster =
						breedingProps->props->population[thisChild].clusterMembership[i] =
								breedingProps->props->population[thisParent2].clusterMembership[i];
				howMany++;
			} else {
				if (breedingProps->props->population[thisChild].clusterMembership[i]
						> breedingProps->props->population[thisParent2].clusterMembership[i]) {
					breedingProps->props->population[thisChild].medoids[i] =
							breedingProps->props->population[thisParent2].medoids[i];
					currCluster =
							breedingProps->props->population[thisChild].clusterMembership[i] =
									breedingProps->props->population[thisParent2].clusterMembership[i];
				}
			}
		}

		//if still 0
		if (breedingProps->props->population[thisChild].clusterMembership[i]
				== 0) {
			breedingProps->props->population[thisChild].clusterMembership[i] =
					currCluster;
			howMany++;
		}
	}

	if ( howMany < breedingProps->props->medoidsVectorSize ) {
		logMessage(" aa %s", "");
	}

	// copy attributes from the first parent
	breedingProps->props->population[thisChild].attr =
			breedingProps->props->population[thisParent1].attr;

	// mutation
	// 1) attributes
	if (breedingProps->table[loop.threadIdx.x].factor > 50) {
		if (breedingProps->table[loop.threadIdx.x].factor > 70) {
			if (breedingProps->table[loop.threadIdx.x].factor > 90) {
				// both
				breedingProps->props->population[thisChild].attr.clusterMaxSize =
						rand() % breedingProps->props->medoidsVectorSize;
				breedingProps->props->population[thisChild].attr.numNeighbours =
						rand() % breedingProps->props->medoidsVectorSize;
			} else {
				// neighbours
				breedingProps->props->population[thisChild].attr.numNeighbours =
						rand() % breedingProps->props->medoidsVectorSize;
			}
		} else {
			// max cluster size
			breedingProps->props->population[thisChild].attr.clusterMaxSize =
					rand() % breedingProps->props->medoidsVectorSize;
		}
	}

	// cleanup cluster membership and sizes
	currCluster = 0;
	clusterToWrite = 1;
	membersCount = 0;
	currMembership =
			breedingProps->props->population[thisChild].clusterMembership[0];

	/*
	 Till membership doesn't change, and size of a cluster is in a norm count a medoid as a member of this cluster.
	 If membership changes or maximum size of a cluster is meet, make a new cluster+
	 */
	for (i = 0; i < breedingProps->props->medoidsVectorSize; i++) {
		if (breedingProps->props->population[thisChild].clusterMembership[i]
				== currMembership
				&& breedingProps->props->population[thisChild].attr.clusterMaxSize
						>= membersCount) {
			membersCount++;
		} else {
			breedingProps->props->population[thisChild].clusters[currCluster++] =
					membersCount;
			currMembership =
					breedingProps->props->population[thisChild].clusterMembership[i];
			membersCount = 1;
			clusterToWrite++;
		}
		breedingProps->props->population[thisChild].clusterMembership[i] =
				clusterToWrite;
		breedingProps->props->population[thisChild].clusters[currCluster] =
				membersCount;
	}

	// mutation
	// 2) medoids
	howMany = 0;
	if (breedingProps->table[loop.threadIdx.x].factor > 20) {
		if (breedingProps->table[loop.threadIdx.x].factor > 50) {
			if (breedingProps->table[loop.threadIdx.x].factor > 80) {
				if (breedingProps->table[loop.threadIdx.x].factor > 95) {
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
	for (i = 0; i < howMany; i++) {
		breedingProps->props->population[thisChild].medoids[rand()
				% breedingProps->props->medoidsVectorSize] = rand()
				% breedingProps->props->dataStore->info.numEntries;
	}

}
//----------------------------------------------

ErrorCode Crossing(BreedingTable * breedingProps) {
	ErrorCode err = errOk;
	LoopDefinition crossingLoop;

	if (breedingProps == NULL || breedingProps->props == NULL
			|| breedingProps->table == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == Crossing %s", "" );
	// <<< 1, popSize/2 >>>
	crossingLoop.gridSize.x = 1;
	crossingLoop.gridSize.y = 1;
	crossingLoop.gridSize.z = 1;
	crossingLoop.blockSize.x = breedingProps->props->popSize / 2;
	crossingLoop.blockSize.y = 1;
	crossingLoop.blockSize.z = 1;
	crossingLoop.kernel = CrossingKernel;
	crossingLoop.params = (void*) breedingProps;

	err = RunLoop(crossingLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void BDIKernel(LoopContext loop);
void BDIKernel(LoopContext loop) {
	unsigned int i;
	unsigned int j;
	float prevDistance;
	float currDistance;
    float * numbers;
    EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + loop.threadIdx.x;
	PopMember *thisMember = props->population + loop.threadIdx.x;
	unsigned int clusterCountedIn = 0;
    
    numbers = (float*)malloc( props->medoidsVectorSize * sizeof(float) );
	// for each pair of medoids taht doesn't belong to the same cluster
	// find the lowest result of:
	// ( densities[ a] + densities[ b] ) / distance( a, b )

	// some cleaning first
	for (i = 0; i < props->medoidsVectorSize; i++) {
		numbers[i] = 0;
	}

	for (i = 0; i < props->medoidsVectorSize; i++) {
		prevDistance = 0; // clean it for new search
		for (j = 0; j < props->medoidsVectorSize; j++) {
            // if not the same record
            if ( (*thisMember).medoids[i] == (*thisMember).medoids[j] ) {
                continue;
            }
			// if not the same and from the same cluster
			if (j != i
					&& thisMember->clusterMembership[i]
							!= thisMember->clusterMembership[j]) {
				currDistance = ((thisSolution->clusterDensities[i]
						+ thisSolution->clusterDensities[j])
						/ Distance(props, (*thisMember).medoids[i],
								(*thisMember).medoids[j]));
				// find the biggest BDI distance
				if (currDistance > prevDistance || prevDistance == 0) {
					prevDistance = currDistance;
				}
			}
		} // for j
		  // save it
		if (numbers[thisMember->clusterMembership[i] - 1] < prevDistance
				|| numbers[thisMember->clusterMembership[i] - 1] == 0) {
			numbers[thisMember->clusterMembership[i] - 1] = prevDistance;
			clusterCountedIn++;
		}
	}

	thisSolution->resultBDI = 0;
	for (i = 0; i < thisSolution->numOfClusters; i++) {
		thisSolution->resultBDI += numbers[i];
	}

	thisSolution->resultBDI /= (float) thisSolution->numOfClusters;
    
    free( numbers );
}
//----------------------------------------------

ErrorCode calculateBDI(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition bdiLoop;

	if (props == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == calculateBDI %s", "" );
	// <<< 1, populationSize >>>
	bdiLoop.gridSize.x = 1;
	bdiLoop.gridSize.y = 1;
	bdiLoop.gridSize.z = 1;
	bdiLoop.blockSize.x = props->popSize;
	bdiLoop.blockSize.y = 1;
	bdiLoop.blockSize.z = 1;
	bdiLoop.kernel = BDIKernel;
	bdiLoop.params = (void*) props;

	err = RunLoop(bdiLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void DIKernel(LoopContext loop);
void DIKernel(LoopContext loop) {
	// Find smallest distance betwen two medoids from different clusters
	// Find biggest density
	// Divide one by another - there you go
	EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + loop.threadIdx.x;
	PopMember *thisMember = props->population + loop.threadIdx.x;

	float currVal = 0;
	float smallestDistance = 0;
	float biggestDensity = 0;
	unsigned int i, j;

	for (i = 0; i < thisSolution->numOfClusters; i++) {
		for (j = 0; j < thisSolution->numOfClusters; j++) {
            // if not the same record
            if ( (*thisMember).medoids[i] == (*thisMember).medoids[j] ) {
                continue;
            }
			// if not the same and from the same cluster
			if (j != i
					&& thisMember->clusterMembership[i]
							!= thisMember->clusterMembership[j]) {
				currVal = Distance(props, (*thisMember).medoids[i],
						(*thisMember).medoids[j]);
				if ( ( currVal != 0 ) && ( ( currVal < smallestDistance )
						|| ( smallestDistance == 0 ) ) ) {
					smallestDistance = currVal;
				}
			}
		}
	}

	for (i = 0; i < thisSolution->numOfClusters; i++) {
		if (thisSolution->clusterDensities[i] > biggestDensity
				|| biggestDensity == 0) {
			biggestDensity = thisSolution->clusterDensities[i];
		}
	}

	thisSolution->resultDI = smallestDistance / biggestDensity;
}
//----------------------------------------------

ErrorCode calculateDI(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition diLoop;

	if (props == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == calculateDI %s", "" );
	// <<< 1, populationSize >>>
	diLoop.gridSize.x = 1;
	diLoop.gridSize.y = 1;
	diLoop.gridSize.z = 1;
	diLoop.blockSize.x = props->popSize;
	diLoop.blockSize.y = 1;
	diLoop.blockSize.z = 1;
	diLoop.kernel = DIKernel;
	diLoop.params = (void*) props;

	err = RunLoop(diLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------

void RandKernel(LoopContext loop);
void RandKernel(LoopContext loop) {
	unsigned int i, j;
	unsigned int t = 0, // true positive and true negative
			f = 0; // false positive and false negative
	char we, they;
	EvolutionProps * props = (EvolutionProps*) loop.params;
	Solution *thisSolution = props->solutions + loop.threadIdx.x;

	for (i = 0; i < props->dataStore->info.numEntries; i++) {
		for (j = 0; j < props->dataStore->info.numEntries; j++) {
			if (thisSolution->recordMembership[i] == 0) {
				//logMessage( " Break %s", "");
			}
			we = (thisSolution->recordMembership[i]
					== thisSolution->recordMembership[j]);
			they =
					(props->dataStore->classes[i]
							== props->dataStore->classes[j]);

			if (we == they) {
				t++;
			} else {
				f++;
			}
		}
	}

	thisSolution->resultRand = (float) t / ((float) t + (float) f);
}
//----------------------------------------------

ErrorCode calculateRand(EvolutionProps * props) {
	ErrorCode err = errOk;
	LoopDefinition randLoop;

	if (props == NULL) {
		return SetLastErrorCode(errWrongParameter);
	}

	logDebug(" == calculateRand %s", "" );
	// <<< 1, populationSize >>>
	randLoop.gridSize.x = 1;
	randLoop.gridSize.y = 1;
	randLoop.gridSize.z = 1;
	randLoop.blockSize.x = props->popSize;
	randLoop.blockSize.y = 1;
	randLoop.blockSize.z = 1;
	randLoop.kernel = RandKernel;
	randLoop.params = (void*) props;

	err = RunLoop(randLoop);

	if (err != errOk) {
		reportError( err, "Run loop returned with error%s", "");
	}
	return err;
}
//----------------------------------------------
//----------------------------------------------
