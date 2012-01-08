


//====================================================================
//== Includes
#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
//#include <math.h>
#include "errors.cuh"
#include "clustering.cuh"
#include "distanceCalculator.cuh"
#include "dataLoader.cuh"
#include <cutil_inline.h>
#include <cutil_math.h>
#include <cuda.h>
#include <time.h>
#include <curand_kernel.h>
//====================================================================
//== Types

//====================================================================
//== Defines
ErrorCode translateCudaError( cudaError cuErr ) {
	ErrorCode err = errOk;

	switch( cuErr ) {
		case cudaErrorMemoryAllocation:
			err = errNoMemory; break;
		case cudaSuccess:
			err = errOk; break;
		default:
			err = errGeneral;
	};

	return err;
}

#define enableErrorControl ErrorCode myErr = errOk;\
	cudaError_t cuErr = cudaSuccess;\
	unsigned int lineNum

#define checkCuda(x) cuErr = x; if ( cuErr != cudaSuccess ) { lineNum = __LINE__; goto cudaErrorLabel; }

// catch and log errors
#define catchCudaError cudaErrorLabel:\
	if ( cuErr != cudaSuccess ) {\
		printf( "[E][cuda] (%d) \"%s\"\n", cuErr, cudaGetErrorString( cuErr ) );\
		return translateCudaError( cuErr );\
	}

//====================================================================
//== Constants and Globals
// Host
static float * dDistances;
static float * dNeighbours;

static unsigned int populationSize;
const char threadsPerBlock = 50;
unsigned int blocksPerEntires = 0;

unsigned int * hDominanceCounts = 0;
bool * hDominanceMatrix = 0;

unsigned int hNumEntries;

// Device
__constant__ unsigned int dNumEntries; //< number of records - database size.
__constant__ unsigned int dThreadsPerBlock; //< how many entries to be calculated in single block.
__constant__ unsigned int dBlocksPerSolution; //< numEntries / threadsPerBlock.
__constant__ unsigned int dPopulationSize; //< size of population for EA.


/*
 * array to hold EA fitnes results.
 * [ populationSize x numObjectives x blocksPerSolution]
 */
__constant__ float * dFitnesResults;
float * fitnesResults;

/*
 * array to hold current population.
 * [ populationSize]
 */
__device__ unit * dPopulationPool;

/*
 * array to hold membership of each record for each solution.
 * [ populationSize x numEntries]
 */
__constant__ char * dMembership;

/*
 * describes dominance relationszhip
 * [ populationSize x populationSize]
 * if [ x, y]==true then x dominates y
 */
bool * dDominanceMatrix;

/*
 * array to hold dominance counts ( how many other solutions are dominating this one ).
 * [ populationSize]
 */
unsigned int * dDominanceCounts;

/*
 * array to hold current Breeding plan for EA.
 * [ populationSize/2 x 4]
 * Where it describes:
 * first parent ( solution number )
 * second parent ( solution number )
 * child ( solution number where resulting child should be placed )
 */
__constant__ breedDescriptor * dBreedingTable;

/*
 * global state for al curand operations
 */
__device__ curandState randState;


//====================================================================
//== Declarations
texture<float, cudaTextureType1D, cudaReadModeElementType> texRefDistances;
texture<float, cudaTextureType1D, cudaReadModeElementType> texRefNeighbour;

// host
ErrorCode runAlgorithms( DataStore * dataStore, unsigned int steps, algResults * results );
void hostRandomPopulation( unsigned int popSize, unit * dPopulationPool );
// device
__global__ void randomPopulation();
__global__ void kernelMembershipAndDensity();
__device__ float distance( unsigned int a, unsigned int b );
__device__ uint distanceIdx(uint x, uint y);
__device__ unsigned int fitnesResultIndex( unsigned int solution, char objective, char block );
__global__ void kernelConnectivity();
__global__ void kernelDisconnectivity();
__global__ void kernelCorectness();
__global__ void kernelSorting( bool * dominanceMatrix );
__global__ void kernelDominanceCount( bool * dominanceMatrix, unsigned int * dominanceCount );
__global__ void kernelFrontDensity( unsigned int * front, unsigned int frontSize, float * frontDensities );
__global__ void kernelCrossing( breedDescriptor * breedingTable );
__global__ void kernelSumResults();

//====================================================================
//== Functions

ErrorCode generateRandomPopulation( unsigned int popSize ) {
	populationSize = popSize;
	cudaMemcpyToSymbol( dPopulationSize, &populationSize, sizeof(unsigned int) );
	unit * populationPool;
	cudaMalloc( &populationPool, popSize * sizeof(unit) );
	cudaMemcpyToSymbol( dPopulationPool, &populationPool, sizeof(unit*) );

	//randomPopulation<<< 1, popSize>>>();
	hostRandomPopulation( popSize, populationPool );
	//cutilDeviceSynchronize();

	cudaError_t cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf("[E][cuda] Error in random population" );
	}

	return errOk;
}
//====================================================================

ErrorCode runClustering( unsigned int popSize, unsigned int steps, DataStore * dataStore, algResults * results ) {
	// Bind textures
	//   Chanel descriptor
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );
	//   Texture properties	
	texRefDistances.addressMode[ 0] = cudaAddressModeWrap;
	texRefDistances.addressMode[ 1] = cudaAddressModeWrap;
	texRefDistances.filterMode = cudaFilterModeLinear;
	texRefDistances.normalized = true;
	texRefNeighbour.addressMode[ 0] = cudaAddressModeWrap;
	texRefNeighbour.addressMode[ 1] = cudaAddressModeWrap;
	texRefNeighbour.filterMode = cudaFilterModeLinear;
	texRefNeighbour.normalized = true;

	//   Allocate memory for distances
	unsigned int offset = 0;

	hNumEntries = dataStore->info.numEntries;
	cudaMemcpyToSymbol( dNumEntries, &hNumEntries, sizeof(unsigned int) );
	
	unsigned int distancesSize = hNumEntries * ( hNumEntries -1 ) /2 * sizeof(float);
	cudaMalloc( &dDistances, distancesSize );
	cudaMemcpy( dDistances, dataStore->distances, dataStore->info.distancesSize, cudaMemcpyHostToDevice );
	//   bind distances to texture
	cudaBindTexture( &offset, &texRefDistances, dDistances, &channelDesc, dataStore->info.distancesSize );

	//   Allocate memory for neighbours	
	unsigned int neighbourSize = hNumEntries * kMaxNeighbours * sizeof(unsigned int);
	cudaMalloc( &dNeighbours, neighbourSize );
	cudaMemcpy( dNeighbours, dataStore->neighbours, neighbourSize, cudaMemcpyHostToDevice );
	//   bind neighbours to texture
	cudaBindTexture( &offset, &texRefNeighbour, dNeighbours, &channelDesc, neighbourSize );

	cudaError_t cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf("[E][cuda:%u] Cuda Error\n" );
	}

	// generate startup population
	generateRandomPopulation( popSize );

	runAlgorithms( dataStore, steps, results );

	return errOk;
}
//====================================================================

void hostRandomPopulation( unsigned int popSize, unit * dPopulationPool ) {
	srand( time( 0 ));

	index2d threadIdx;

	unit * populationPool = (unit*)malloc( popSize * sizeof(unit) );

//	bool proposalOk = false;

	for ( threadIdx.x = 0;  threadIdx.x < populationSize; threadIdx.x++ ) {
		// attributes
		populationPool[ threadIdx.x].attr.clusterMaxSize = rand() % MAX_CLUSTER_SIZE + 1;
		populationPool[ threadIdx.x].attr.numNeighbours = rand() % kMaxNeighboursToUSe + 1;

		// medoids and clusters
		unsigned int clustersSum = MEDOID_VECTOR_SIZE;
		unsigned int proposal;
		bool proposalOk = false;
		for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
			do {
				proposalOk = true;
				proposal = rand() % hNumEntries;
				for ( int j = 0; j < i; j++ ) {
					if ( proposal == populationPool[ threadIdx.x].medoids[ j] ) {
						proposalOk = false;						
					}
				}
			} while ( !proposalOk );
			populationPool[ threadIdx.x].medoids[ i] = proposal;

			if ( clustersSum > 0 ) {
				proposal = rand() % populationPool[ threadIdx.x].attr.clusterMaxSize + 1;
				if ( clustersSum < proposal ) {
					proposal += clustersSum - proposal;
				}
				clustersSum -= proposal;
				populationPool[ threadIdx.x].clusters[ i] = proposal;
			} else {
				populationPool[ threadIdx.x].clusters[ i] = 0;
			}
		} // for each medoid in vector
	}

	cudaMemcpy( dPopulationPool, populationPool, popSize * sizeof(unit), cudaMemcpyHostToDevice );
}
//====================================================================


// <<< 1, populationSize >>>
__global__ void kernelRandomPopulation() {
	if ( threadIdx.x == 0 ) {
		curand_init( dPopulationSize, dNumEntries, 0, &randState );
	}

	// attributes
	dPopulationPool[ threadIdx.x].attr.clusterMaxSize = curand( &randState ) % MAX_CLUSTER_SIZE + 1;
	dPopulationPool[ threadIdx.x].attr.numNeighbours = curand( &randState ) % kMaxNeighboursToUSe + 1;

	// medoids and clusters
	unsigned int clustersSum = MEDOID_VECTOR_SIZE;
	unsigned int proposal;
	bool proposalOk = false;
	for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		do {
			proposalOk = true;
			proposal = curand( &randState ) % dNumEntries;
			for ( int j = 0; j < i; j++ ) {
				if ( proposal == dPopulationPool[ threadIdx.x].medoids[ j] ) {
					proposalOk = false;
					break;
				}
			}
		} while ( !proposalOk );
		dPopulationPool[ threadIdx.x].medoids[ i] = proposal;

		if ( clustersSum > 0 ) {
			proposal = curand( &randState ) % dPopulationPool[ threadIdx.x].attr.clusterMaxSize + 1;
			if ( clustersSum < proposal ) {
				proposal += clustersSum - proposal;
			}
			clustersSum -= proposal;
			dPopulationPool[ threadIdx.x].clusters[ i] = proposal;
		} else {
			dPopulationPool[ threadIdx.x].clusters[ i] = 0;
		}
	} // for each medoid in vector
}
//====================================================================

ErrorCode configureAlgoritms() {
	// Set things up
	//float * fitnesResults = 0;
	char * membership = 0;
//	cudaError err = cudaSuccess;

	blocksPerEntires = hNumEntries / threadsPerBlock;
	while ( blocksPerEntires * threadsPerBlock <  hNumEntries ) {
		blocksPerEntires++;
	}

	enableErrorControl;

	checkCuda( cudaMemcpyToSymbol( dBlocksPerSolution, &blocksPerEntires, sizeof(unsigned int) ));
	checkCuda( cudaMemcpyToSymbol( dThreadsPerBlock, &threadsPerBlock, sizeof(unsigned int) ));

	// populationSize x numObjectives x blocksPerSolution
	checkCuda( cudaMalloc( &fitnesResults, populationSize * OBJECTIVES * blocksPerEntires * sizeof(float)));
	checkCuda( cudaMemcpyToSymbol( dFitnesResults, &fitnesResults, sizeof(float*) ));
	checkCuda( cudaMalloc( &membership, populationSize * hNumEntries * sizeof(char) ));
	checkCuda( cudaMemcpyToSymbol( dMembership, &membership, sizeof(char*) ));

	// population dominations
	bool * dominanceMatrix = 0;
	unsigned int * dominanceCounts = 0;
	checkCuda( cudaMalloc( &dDominanceMatrix, populationSize * populationSize * sizeof(bool) ));	
	checkCuda( cudaMalloc( &dDominanceCounts, populationSize * sizeof(unsigned int) ));
	
	//=-

	hDominanceMatrix = (bool*)malloc( populationSize * populationSize * sizeof(bool));
	hDominanceCounts = (unsigned int*)malloc( populationSize * sizeof(unsigned int));
	
	return errOk;

	catchCudaError;
}

ErrorCode runAlgorithms( DataStore * dataStore, unsigned int steps, algResults * results ) {

	enableErrorControl;

	// Set ups
	myErr = configureAlgoritms();	

	if (myErr != errOk) {
		return myErr;
	}

	dim3 dimGrid( blocksPerEntires, populationSize );
	//dim3 dimBlock( threadsPerBlock );

	//dim3 dimGrid2( populationSize );
	//dim3 dimBlock2( MEDOID_VECTOR_SIZE );

	unsigned int solutionsLeft = 0;
	bool * solutionsSelected = (bool*)malloc( populationSize * sizeof(bool) );
	unsigned int * solutionFronts = (unsigned int*)malloc( populationSize * ( populationSize + 1 ) * sizeof(unsigned int));
	unsigned int currFront;
	unsigned int currFrontSize;
	
	float * dFrontDensities = 0;
	cudaMalloc( &dFrontDensities, populationSize * sizeof(float) );	
	float * hFrontDensities = (float*)malloc( populationSize * sizeof(float) );

	breedDescriptor * breedingTable;
	unsigned int halfPopulation = populationSize / 2;
	cudaMalloc( &breedingTable, halfPopulation * sizeof(breedDescriptor));
	cudaMemcpyToSymbol( dBreedingTable, &breedingTable, sizeof(breedDescriptor*) );

	for (int i = 0; i < steps; i++ ) {
		// membership and density phase
		kernelMembershipAndDensity<<<dimGrid, threadsPerBlock>>>();
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cuda] After kernelMembershipAndDensity - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// connectivity phase
		kernelConnectivity<<<dimGrid, threadsPerBlock>>>();
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cuda] After kernelConnectivity - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// sum up results		
		kernelSumResults<<< populationSize,  2 >>>();
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cuda] After kernelSumResults - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// disconnectivity phase
		kernelDisconnectivity<<<populationSize, MEDOID_VECTOR_SIZE>>>();
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cuda] After kernelDisconnectivity - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// correctness phase
		kernelCorectness<<<populationSize, MEDOID_VECTOR_SIZE>>>();
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cude] After kernelCorectness - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		float * hResults = (float*)malloc( populationSize * OBJECTIVES * blocksPerEntires * sizeof(float) );
		cudaMemcpy( hResults, fitnesResults, populationSize * OBJECTIVES * blocksPerEntires * sizeof(float), cudaMemcpyDeviceToHost );
		
		// sorting
		kernelSorting<<<populationSize, populationSize>>>( dDominanceMatrix );
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cude] After kernelSorting - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// dominance count
		kernelDominanceCount<<<1, populationSize>>>( dDominanceMatrix, dDominanceCounts );
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();

		if ( cuErr != cudaSuccess ) {
			printf( "[E][cude] After kernelDominanceCount - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		cudaMemcpy( hDominanceMatrix, dDominanceMatrix, populationSize * populationSize * sizeof(bool), cudaMemcpyDeviceToHost );
		cudaMemcpy( hDominanceCounts, dDominanceCounts, populationSize * sizeof( unsigned int ), cudaMemcpyDeviceToHost );

		cuErr = cudaGetLastError();
		if ( cuErr != cudaSuccess ) {
			printf( "[E] After Dominance memory copies - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// setup fronts
		solutionsLeft = populationSize;
		currFront = 0;
		
		int j;
		for ( j = 0; j < populationSize; j++ ) {
			solutionsSelected[ j] = false;
		}

		// front grouping phase
		while ( solutionsLeft > 0 ) {
			currFrontSize = 0;
			// select solutions for current front - where domination count is 0
			for ( j = 0; j < populationSize && solutionsLeft > 0; j++ ) {				
				if ( !solutionsSelected[ j] && hDominanceCounts[ j] == 0 ) {
					solutionFronts[ currFront * ( populationSize + 1 ) + (++currFrontSize)] = j;
					solutionsSelected[ j] = true;
					solutionsLeft--;
				}
			}
			solutionFronts[ currFront * ( populationSize + 1 ) + 0] = currFrontSize;
						
			if ( solutionsLeft > 0 ) {
				// for each solution dominated by solution from this front - reduce domination count
				for ( j = 0; j < currFrontSize; j++ ) {
					for ( int k = 0; k < populationSize; k++ ) {
						if ( hDominanceMatrix[ solutionFronts[ currFront * ( populationSize + 1 ) + j + 1] * populationSize + k] && hDominanceCounts[ k] > 0 ) {
							hDominanceCounts[ k] -= 1;
						}
					}
				}
			}

			if ( currFrontSize == 0 ) {
				for( int k = 0; k < populationSize; k++ ) {
					if ( hDominanceCounts[ k] > 0 ) {
						hDominanceCounts[ k] -= 1;
					}
				}
			} else {
				// now for next front
				currFront++;
			}
		}

		// selection
		solutionsLeft = populationSize / 2; // select half size of population
		for ( j = 0; j < populationSize; j++ ) {
			solutionsSelected[ j] = false;
		}

		currFront = 0;
		while ( solutionsLeft > 0 ) {
			// if we need more than the current front can offer
			if ( solutionsLeft >= solutionFronts[ currFront * ( populationSize + 1 ) + 0] ) {
				for ( j = 0; j < solutionFronts[ currFront * ( populationSize + 1 ) + 0]; j++ ) {
					solutionsSelected[ solutionFronts[ currFront * ( populationSize + 1 ) + j + 1]] = true;
					solutionsLeft--;
				}
			} else {
				// this front has more than we need
				currFrontSize = solutionFronts[ currFront * ( populationSize + 1 ) + 0];

				unsigned int * dFront;

				cudaMalloc( &dFront, currFrontSize * sizeof(unsigned int) );

				cudaMemcpy( dFront, &solutionFronts[ currFront * ( populationSize + 1 ) + 1], currFrontSize * sizeof(unsigned int), cudaMemcpyHostToDevice );

				// Calculate densities for solutions in this front
				kernelFrontDensity<<<currFrontSize, 4>>>( dFront, currFrontSize, dFrontDensities );

				cutilDeviceSynchronize();
				cudaFree( dFront );

				cuErr = cudaGetLastError();
				if ( cuErr != cudaSuccess ) {
					printf( "[E][cude] After kernelFrontDensity - %s\n", cudaGetErrorString( cuErr ));
					break;
				}

				// Export results to Host
				cudaMemcpy( hFrontDensities, dFrontDensities, currFrontSize * sizeof(float), cudaMemcpyDeviceToHost );
				cuErr = cudaGetLastError();
				if ( cuErr != cudaSuccess ) {
					printf( "[E] After copying front densities to host - %s\n", cudaGetErrorString( cuErr ));
					break;
				}
				
				bool * thisFrontSelection = (bool*)malloc( populationSize * sizeof(bool));
				unsigned int smallest = 0;

				// Select first selectionLeft solutions and find the smallest one (bug density)
				smallest = 0;
				for ( j = 0; j < currFrontSize; j++ ) {
					thisFrontSelection [ j] = ( j < solutionsLeft );
					if ( thisFrontSelection[ j] ) {
						if ( hFrontDensities[ j] != -1 && ( hFrontDensities[ j] < hFrontDensities[ smallest] || hFrontDensities[ smallest] == -1 ) ) {
							smallest = j;
						}
					}
				} // for j

				// Now for each solution not selected at first, check if it's bigger than the smallest
				// If so, find new smallest
				if  ( hFrontDensities[ smallest] != -1 ) {
					for ( j = 0; j < currFrontSize; j++ ) {
						if ( thisFrontSelection[ j] ) continue;

						if ( hFrontDensities[ j] == -1 || hFrontDensities[ j] > hFrontDensities[ smallest] ) {
							thisFrontSelection[ smallest] = false;
							thisFrontSelection[ j] = true;
							smallest = j;
							for ( int k = 0; k < j; k++ ) {
								if ( thisFrontSelection[ k] ) {
									if ( hFrontDensities[ k] != -1 && ( hFrontDensities[ k] < hFrontDensities[ smallest] || hFrontDensities[ smallest] == -1 ) ) {
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
						solutionsSelected[ solutionFronts[ currFront * ( populationSize + 1 ) + j + 1]] = true;
						solutionsLeft--;
					}
				}// for j
			}

			currFront++;
		} // while

		// crossing
		// breedingTable[ parent1, parent2, child, mutation probability]
		breedDescriptor * hBreedingTable = (breedDescriptor*)malloc( halfPopulation * sizeof(breedDescriptor) );
		unsigned int currParent1 = 0;
		unsigned int currParent2 = 0;
		unsigned int currChild = 0;

		srand( clock() );
		// generate breeding Table
		for ( j = 0; j < populationSize; j++ ) {
			if ( solutionsSelected[ j] ) {
				// place for parent
				if ( currParent1 <= currParent2 ) {
					// place taken by first parent
					hBreedingTable[ currParent1++].parent1 = j;
					hBreedingTable[ currParent1++].parent2 = j;
				} else {
					hBreedingTable[ currParent2++].parent2 = j;
					hBreedingTable[ currParent2++].parent1 = j;
				}
			} else {
				// place for child
				hBreedingTable[ currChild].child = j;
				// mutation probability
				hBreedingTable[ currChild++].factor = rand() % 100; 
			}
		} // for j

		cudaMemcpy( breedingTable, hBreedingTable, halfPopulation * sizeof(breedDescriptor), cudaMemcpyHostToDevice );
		cuErr = cudaGetLastError();
		if ( cuErr != cudaSuccess ) {
			printf( "[E] After copying breeding table to device - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

		// launch crossing
		kernelCrossing<<< 1, halfPopulation >>>( breedingTable );
		cutilDeviceSynchronize();
		cuErr = cudaGetLastError();
		if ( cuErr != cudaSuccess ) {
			printf( "[E][cude] After kernelCrossing - %s\n", cudaGetErrorString( cuErr ));
			break;
		}

	} // evolution

	
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] After evolution - %s\n", cudaGetErrorString( cuErr ));
	}

	calculateBDI( results->bdi, results->clusters );
	calculateDI( results->di );
	calculateRand( dataStore, results->rand );

	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] After evolution calculation - %s\n", cudaGetErrorString( cuErr ));
	}

	return errOk;
}
//====================================================================
// <<<( blocks, populationSize ), ( threadsPerBlock ) >>>
__global__ void kernelMembershipAndDensity() {
	unsigned int solution = blockIdx.y;
	unsigned int record = blockIdx.x * threadsPerBlock + threadIdx.x;
	unsigned int zeros = 0;

	__shared__ unit thisSolution;
	__shared__ float density[ 256]; // shared table to hold density results for futher calculation

	if ( threadIdx.x == 0 ) {
		// first thread should load solution for this block to shared memory
		thisSolution = dPopulationPool[ solution];
	}

	// Sync up threads
	__syncthreads();

	density[ threadIdx.x] = 0;

	float prevDistance = distance( record, thisSolution.medoids[ 0] ) + 0.1;
	float currDistance;
	unsigned int res = 0;
	
	for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		currDistance = distance ( record, thisSolution.medoids[ i] );
		if ( ( currDistance < prevDistance ) && ( currDistance != 0.0 ) ) {
			prevDistance = currDistance;
			res = i;
		}
		if ( currDistance == 0.0 ) {
			zeros++;
		}
	}

	density[ threadIdx.x] = prevDistance / (float)(dNumEntries - zeros);

	dMembership[ solution * dNumEntries + record] = res;

	// Sync up threads
	__syncthreads();

	if ( threadIdx.x == 0 ) {
		// sum solutions from all threads in this block
		currDistance = 0;
		for ( int i = 0; i < dThreadsPerBlock; i++ ) {
			currDistance += density[ 0];
		}

		// sum all solutions for this block
		dFitnesResults[ fitnesResultIndex( solution, 0, blockIdx.x )] = currDistance;		
	}	
}
//====================================================================

__device__ float distance( unsigned int a, unsigned int b ) {
	if ( a == b ) return 0.0f;
	return tex1Dfetch( texRefDistances, distanceIdx( a, b ));
}
//====================================================================

// [ populationSize x objectives x blocksPerSolution ]
__device__ unsigned int fitnesResultIndex( unsigned int solution, char objective, char block ) {
	return ( solution * OBJECTIVES * dBlocksPerSolution +
		objective * dBlocksPerSolution +
		block );
}
//====================================================================

__device__ uint distanceIdx(uint x, uint y) {
	if ( y > x ) {
		return y * (y - 1) / 2 + x;
	} else {		
		return x * (x - 1) / 2 + y;
	}
}
//====================================================================

__device__ unsigned int neighbour( unsigned int record, unsigned int num ) {
	return tex1Dfetch( texRefNeighbour, record * kMaxNeighbours + num );
}
//====================================================================

__global__ void kernelConnectivity() {
	
	unsigned int solution = blockIdx.y;
	unsigned int record = blockIdx.x * dThreadsPerBlock + threadIdx.x;
	unit thisSolution = dPopulationPool[ solution];
	unsigned int memberOf = dMembership[ solution * dNumEntries + record];
	unsigned int numOfNeighbours = thisSolution.attr.numNeighbours;

	__shared__ float blockResults[256];
	float result = 0;

	for ( int i = 0; i < numOfNeighbours; i++ ) {
		if ( memberOf == dMembership[ solution * dNumEntries + i] ) {
			result += 1.0f / (float)numOfNeighbours;
		}
	}

	blockResults[ threadIdx.x] = result;

	// Sync and sum solutions for this block
	__syncthreads();
	if ( threadIdx.x == 0 ) {
		// sum solutions from all threads in this block
		result = 0;
		for ( int i = 0; i < dThreadsPerBlock; i++ ) {
			result += blockResults[ i];
		}

		// sum all solutions for this block
		dFitnesResults[ fitnesResultIndex( solution, 1, blockIdx.x)] = result;
	}
}
//====================================================================

// <<< populationSize,  2 >>>
__global__ void kernelSumResults() {
	float result = 0;
	
	for ( int i = 0; i < dBlocksPerSolution; i++ ) {
		result += dFitnesResults[ fitnesResultIndex( blockIdx.x, threadIdx.x, i )];
	}

	dFitnesResults[ fitnesResultIndex( blockIdx.x, threadIdx.x, 0 )] = result;
}
//====================================================================

// <<< populationSize, medoidsVectorSize >>>
__global__ void kernelDisconnectivity() {

	__shared__ unsigned int clusters[ MEDOID_VECTOR_SIZE];
	__shared__ float disconnectivities[ MEDOID_VECTOR_SIZE];

	unsigned int comparisions;
	unsigned int counts;

	// label medoids to their clusters
	if ( threadIdx.x == 0 ) {
		for ( int i=0, j=0; i<MEDOID_VECTOR_SIZE && j<MEDOID_VECTOR_SIZE; i++ ) {
			for (int k=0; k<dPopulationPool[ blockIdx.x].clusters[i]; k++,j++) {
				clusters[ j] = i;
			}
		}
	}

	__syncthreads();

	float currDistance  = 0;
	
	comparisions = 1; // start with 1 as we compare those two medoids
	counts = 2; // start with 2 for both medoids
	// For each medoid in the vector
	for ( unsigned int i=0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( i == threadIdx.x || clusters[ threadIdx.x] == clusters[ i] ) {
			// if medoid the same or same cluster - skip
			continue;
		}
		comparisions++;
		currDistance = distance( dPopulationPool[ blockIdx.x].medoids[ threadIdx.x], dPopulationPool[ blockIdx.x].medoids[ i] );
		// now find MND for this medoid
		for ( unsigned int j = 0; j < MEDOID_VECTOR_SIZE; j++ ) {
			if ( j == i || j == threadIdx.x ) {
				continue;
			}
			// counts if someone else is closer to our medoid
			if ( distance( dPopulationPool[ blockIdx.x].medoids[ threadIdx.x],
				dPopulationPool[ blockIdx.x].medoids[ j] ) < currDistance )  {
				counts++;
			}
			// counts if someone else is closer to the other
			if ( distance( dPopulationPool[ blockIdx.x].medoids[ i],
				dPopulationPool[ blockIdx.x].medoids[ j] ) < currDistance ) {
				counts++;
			}
		}
	}

	disconnectivities[ threadIdx.x] += counts / comparisions;

	__syncthreads();

	if ( threadIdx.x == 0 ) {
		float disconnectivity = 0;
		for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
			// CHANGED 12-12-2011: 
			//-- compars += comparisions[ i];
			//-- count += counts[ i];
			disconnectivity += disconnectivities[ i];
		}
		//-- dFitnesResults[ fitnesResultIndex( blockIdx.x, 2, 0)] = count/compars;
		dFitnesResults[ fitnesResultIndex( blockIdx.x, 2, 0)] = disconnectivity;
	}
}
//====================================================================

// <<< populationSize, medoidsVectorSize >>>
__global__ void kernelCorectness() {

	__shared__ char checks[ MEDOID_VECTOR_SIZE];
	__shared__ unsigned int medoids[ MEDOID_VECTOR_SIZE];

	if ( threadIdx.x == 0 ) {
		memcpy( medoids, dPopulationPool[ blockIdx.x].medoids, MEDOID_VECTOR_SIZE * sizeof(unsigned int) );
	}

	__syncthreads();

	unsigned int thisMedoid = medoids[ threadIdx.x];
	char count = 0;
	
	for ( int i=0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if (i == threadIdx.x) {
			continue;
		}

		if ( thisMedoid == medoids[ i] ) {
			count += 2;
		}
	}

	checks[ threadIdx.x] = count;

	__syncthreads();

	count = 0;
	if ( threadIdx.x == 0 )  {
		for ( int i = 0; i < dThreadsPerBlock; i++ ) {
			count += checks[ i];
		}
		dFitnesResults[ fitnesResultIndex( blockIdx.x, 3, 0 )] = count;
	}
}
//====================================================================

// <<< populationSize, populationSize >>>
__global__ void kernelSorting( bool * dominanceMatrix ) {

	__shared__ float thisSolutionFitnesResults[ 4];

	// true if this solution (blockIdx.x) dominates the other one (threadIdx.x)
	bool currDominance = false;
	float heCurrResult;
	float meCurrResult;
	bool hasBetter = false;
	bool hasWorse = false;

	if ( threadIdx.x == 0 ) {		
		thisSolutionFitnesResults[ 0] = dFitnesResults[ fitnesResultIndex( blockIdx.x, 0, 0)];
		thisSolutionFitnesResults[ 1] = dFitnesResults[ fitnesResultIndex( blockIdx.x, 1, 0)];
		thisSolutionFitnesResults[ 2] = dFitnesResults[ fitnesResultIndex( blockIdx.x, 2, 0)];
		thisSolutionFitnesResults[ 3] = dFitnesResults[ fitnesResultIndex( blockIdx.x, 3, 0)];
	}

	__syncthreads();

	for ( int i = 0; i < 1; i++ ) {
		if ( hasWorse && hasBetter || threadIdx.x == blockIdx.x ) {
			// theyre already "equal" stop comparing
			break;
		}
		heCurrResult = dFitnesResults[ fitnesResultIndex( threadIdx.x, i, 0 )];
		meCurrResult = dFitnesResults[ fitnesResultIndex( blockIdx.x, i, 0 )];
		if ( heCurrResult != meCurrResult) {
			switch ( i ) {
				case 0: {// Density
					if ( heCurrResult < meCurrResult ) {
						hasWorse = true;						
					} else {
						hasBetter = true;
					}
				} break;
				case 3: {// Correctnes
					// smaller better
					if ( heCurrResult < meCurrResult ) {
						hasWorse = true;						
					} else {
						hasBetter = true;
					}
				} break;
				case 1: {// Connectivity
					if ( heCurrResult > meCurrResult ) {
						hasWorse = true;						
					} else {
						hasBetter = true;
					}
				} break;
				case 2: { // Disconnectivity
					// bigger better
					if ( heCurrResult > meCurrResult ) {
						hasWorse = true;						
					} else {
						hasBetter = true;
					}
				} break;
			}; // switch
		} // if
	} // for

	if ( hasBetter && !hasWorse ) {
		currDominance = true;
	}

	// if blockIdx.x dominates over threadIdx.x then true
	dominanceMatrix[ blockIdx.x * dPopulationSize + threadIdx.x] = currDominance;
}
//====================================================================

// <<< 1, poulationSize >>>
__global__ void kernelDominanceCount( bool * dominanceMatrix, unsigned int * dominanceCount ) {

	unsigned int count = 0;

	for ( int i = 0; i < dPopulationSize; i++ ) {
		if ( dominanceMatrix[ i * dPopulationSize + threadIdx.x] ) {
			// i dominates over threadIdx.x
			count ++ ;
		}
	}

	dominanceCount[ threadIdx.x] = count;
}
//====================================================================

// <<< numSolutions, kryterions >>>
__global__ void kernelFrontDensity( unsigned int * front, unsigned int frontSize, float * frontDensities ) {

	__shared__ float solutionDensities [ 4];

	bool lesserFound = false;
	float lesserResult;
	bool biggerFound = false;
	float biggerResult;

	float thisResult = dFitnesResults[ fitnesResultIndex( front[ blockIdx.x], threadIdx.x, 0 )];
	float currResult;

	

	for ( int i = 0; i < frontSize; i++ ) {
		if ( blockIdx.x == i ) {
			// skip if same
			continue;
		}

		currResult = dFitnesResults[ fitnesResultIndex( front[ i], threadIdx.x, 0)];
		// check if lesser
		if ( thisResult > currResult ) {
			if ( !lesserFound ) {
				//lesser = i;
				lesserFound = true;
				lesserResult = currResult;
			} else {
				if ( lesserResult < currResult ) {
					//lesser = i;
					lesserFound = true;
					lesserResult = currResult;
				}
			}
		}

		// check if bigger
		if ( thisResult < currResult ) {
			if ( !biggerFound ) {
				//bigger = i;
				biggerFound = true;
				biggerResult = currResult;
			} else {
				if ( biggerResult > currResult ) {
					//bigger = i;
					biggerResult = currResult;
				}
			}
		}
	} // for each solution in this front

	// is this edge solution ?
	if ( !lesserFound || !biggerFound ) {
		solutionDensities[ threadIdx.x] = -1;  
	} else {
		solutionDensities[ threadIdx.x] = biggerResult - lesserResult;
	}

	if ( threadIdx.x == 0 ) {
		for ( int i = 1; i < 4; i++ ) {
			if ( solutionDensities[ i] != -1 ) {
				solutionDensities[ 0] += solutionDensities[ i];
			} else {
				solutionDensities[ 0] = -1;
				break;
			}
		}
		frontDensities[ blockIdx.x] = solutionDensities[ 0];
	}
}
//====================================================================

__global__ void kernelCrossing( breedDescriptor * breedingTable ) {
	char parent1Clusters[ MEDOID_VECTOR_SIZE];
	char parent2Clusters[ MEDOID_VECTOR_SIZE];
	
	breedDescriptor descriptor = breedingTable[ threadIdx.x];	
	
	curandState randState;
	curand_init( dPopulationSize, dNumEntries, 0, &randState );	

	// calculate cluster groupings for parent 1
	unsigned int toSign = 0;
	unsigned int j = 0;
	for ( unsigned int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( toSign == 0 ) {
			toSign = dPopulationPool[ descriptor.parent1].clusters[ j++];			
		}
		parent1Clusters[ i] = j - 1;
		if ( toSign != 0 ) toSign--;
	}

	toSign = 0;
	j = 0;
	for ( unsigned int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( toSign == 0 ) {
			toSign = dPopulationPool[ descriptor.parent2].clusters[ j++];			
		}
		parent2Clusters[ i] = j - 1;
		if ( toSign != 0 ) toSign--;
	}
		
	char childCluster[ MEDOID_VECTOR_SIZE];
	unit childUnit;
	unsigned char howMany = 0;
	unsigned int currCluster = 0;

	// exchange data and put it into child
	for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		childCluster[ i] = parent1Clusters[ i];
		childUnit.medoids[ i] = dPopulationPool[ descriptor.parent1].medoids[ i];

		if ( parent1Clusters[ i] % 2 ) {
			currCluster = childCluster[ i] = parent1Clusters[ i];
			howMany++;
		}

		if ( !( parent2Clusters[ i] % 2 ) ) {
            if ( childCluster[ i] == 0 ) {
                childUnit.medoids[ i] = dPopulationPool[ descriptor.parent2].medoids[ i];
                currCluster = childCluster[ i] = parent2Clusters[ i];
                howMany++;
            } else {
                if ( childCluster[ i] > parent2Clusters[ i] ) {
                    childUnit.medoids[ i] = dPopulationPool[ descriptor.parent2].medoids[ i];
                    currCluster = childCluster[ i] = parent2Clusters[ i];
                }
            }
        }

		//if still 0
		if ( childCluster[ i] == 0 ) {
			childCluster[ i] = currCluster;
			howMany++;
		}
	} // for

	// copy attributes from first parent
	childUnit.attr.clusterMaxSize = dPopulationPool[ descriptor.parent1].attr.clusterMaxSize;
	childUnit.attr.numNeighbours = dPopulationPool[ descriptor.parent1].attr.numNeighbours;

	// mutation
	char mutationProb = descriptor.factor;	
	// 1) attributes
	if ( mutationProb > 50 ) {
		if ( mutationProb > 70 ) {
			if ( mutationProb > 90 ) {
				// both
				childUnit.attr.clusterMaxSize = curand( &randState ) % MAX_CLUSTER_SIZE;
				childUnit.attr.numNeighbours = curand( &randState ) % kMaxNeighbours;
			} else {
				// neighbours
				childUnit.attr.numNeighbours = curand( &randState ) % kMaxNeighbours;
			}
		} else {
			// max cluster size
			childUnit.attr.clusterMaxSize = curand( &randState ) % MAX_CLUSTER_SIZE;
		}
	}

	// 2) medoids
	if ( mutationProb > 20 ) {
		if ( mutationProb > 50 ) {
			if ( mutationProb > 80 ) {
				if ( mutationProb > 95 ) {
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
	for ( char i = 0; i < howMany; i++ ) {
		childUnit.medoids[ curand( &randState ) % MEDOID_VECTOR_SIZE] = curand( &randState ) % dNumEntries;
	}
	
	// check groups sizes, if needed regroup
	unsigned int index = 0;
	char cluster = childCluster[ 0];
	childUnit.clusters[ 0] = 0;
	for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( childCluster[ i] == cluster && childUnit.clusters[ index] < childUnit.attr.clusterMaxSize ) {
			childUnit.clusters[ index]++;
		} else {
			childUnit.clusters[ ++index] = 1;
			cluster = childCluster[ i];
		}
	}
	for ( int i = index+1; i < MEDOID_VECTOR_SIZE; i++ ) {
		childUnit.clusters[ i] = 0;
	}

	// now save the child into memory
	dPopulationPool[ descriptor.child] = childUnit;
}
//====================================================================

__device__ unsigned int devCalculateClustersAndDensities( unsigned char * clusters, float * densities ) {
	// which medoid to what cluster
	float numbers[ MEDOID_VECTOR_SIZE]; //< number of entries per cluster


	// calculate cluster groupings for p
	unsigned int toSign = 0;
	unsigned int j = 0;
	unsigned int i = 0;
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( toSign == 0 ) {
			toSign = dPopulationPool[ threadIdx.x].clusters[ j++];			
		}
		clusters[ i] = j - 1;
		if ( toSign != 0 ) toSign--;
	}
	
	// == 1 == Densities
	// first clear the array
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		densities[ i] = 0;
		numbers[ i] = 0;
	}

	// now calculate it
	for ( i = 0; i < dNumEntries; i++ ) {
		densities[ clusters[ dMembership[ threadIdx.x * dNumEntries + i]]] += distance( i, dMembership[threadIdx.x * dNumEntries + i] );
		numbers[ clusters[ dMembership[ threadIdx.x * dNumEntries + i]]] += 1;
	}

	// means
	for ( i = 0; i < j; i++ ) {
		densities[ i] /= (float)numbers[ i];
	}

	// return number of clusters in this solution
	return j;
}

// <<< 1, dPopulationSize >>>
__global__ void kernelBDI( float * indexes, unsigned char * clustersCount ) {

	// which medoid to what cluster
	unsigned char clusters[ MEDOID_VECTOR_SIZE]; //< what medoid to what cluster
	float densities[ MEDOID_VECTOR_SIZE]; //< densities of each cluster
	float numbers[ MEDOID_VECTOR_SIZE]; //< number of entries per cluster


	// calculate cluster groupings for p
	//unsigned int toSign = 0;
	unsigned int j = 0;
	unsigned int i = 0;
	
	j = devCalculateClustersAndDensities( clusters, densities );
	
	// == 2 == Distances betwen clusters
	unsigned int k = j;
	float prevDistance;
	float currDistance;

	// first clear the array once more
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {		
		numbers[ i] = 0;
	}

	// for each medoid
	for ( i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
		// compare it withe every else medoid
		prevDistance = 0;
		for ( j = 0 ; j < MEDOID_VECTOR_SIZE; j++ ) {
			// if not the same and from the same cluster
			if ( j != i && clusters[ i] != clusters[ j] &&  densities[ i] != 0 &&  densities[ j] != 0) {
				currDistance = (( densities[ i] + densities[ j] ) /
					distance( dPopulationPool[ threadIdx.x].medoids[ i],
					dPopulationPool[ threadIdx.x].medoids[ j] ));

				// find the biggest one
				if ( currDistance > prevDistance || prevDistance == 0 ) {
					prevDistance = currDistance;
				}
			}
		} // for j

		// save it
		if ( numbers[ clusters[ i]] < prevDistance || numbers[ clusters[ i]] == 0 ) {
			numbers[ clusters[ i]] = prevDistance;
		}
	} // for i
	
	// Now calculate the index
	currDistance = 0;
	for( i = 0; i < k; i++ ) {
		currDistance += numbers[ i];
	}
	currDistance /= (float)k;
	
	// return results
	indexes[ threadIdx.x] = currDistance;
	clustersCount[ threadIdx.x] = k;
}
//====================================================================
ErrorCode calculateBDI( preciseResult & topBDI, preciseResult & clusters ) {
	float * dResults;
	unsigned char * dClusters;
	enableErrorControl;

	cudaMalloc( &dResults, populationSize * sizeof(float) );
	cudaMalloc( &dClusters, populationSize * sizeof(unsigned int) );

	kernelBDI<<< 1, populationSize >>>( dResults, dClusters );
	cutilDeviceSynchronize();
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] After kernelBDI - %s\n", cudaGetErrorString( cuErr ));
	}

	float * hResults = (float*)malloc( populationSize * sizeof(float) );
	unsigned char * hClusters = (unsigned char*)malloc( populationSize * sizeof(unsigned int) );
	cudaMemcpy( hResults, dResults, populationSize * sizeof(float), cudaMemcpyDeviceToHost );
	cudaMemcpy( hClusters, dClusters, populationSize * sizeof(unsigned int), cudaMemcpyDeviceToHost );

	//printf( " \n BDI index: \n" );
	topBDI.min = topBDI.max = hResults[ 0];
	topBDI.sum = 0.0;
	clusters.min = clusters.max = hClusters[ 0];
	clusters.sum = 0.0;
	for ( int i = 1; i < populationSize; i++ ) {
		if ( topBDI.min == 0 || hResults[ i] < topBDI.min ) {
			topBDI.min = hResults[ i];
		}
		if ( topBDI.max == 0 || hResults[ i] > topBDI.max ) {
			topBDI.max = hResults[ i];
		}
		topBDI.sum += hResults[ i];

		if ( clusters.min == 0 || hClusters[ i] < clusters.min ) {
			clusters.min = hClusters[ i];
		}
		if ( clusters.max == 0 || hClusters[ i] > clusters.max ) {
			clusters.max = hClusters[ i];
		}
		clusters.sum += hClusters[ i];
	}

	return errOk;
}
//====================================================================

// <<< 1, dPopulationSize >>>
__global__ void kernelCalculateDI( float * indexes ) {
	unsigned char clusters[ MEDOID_VECTOR_SIZE];
	float densities[ MEDOID_VECTOR_SIZE];
	float lowestDensity = 0;

	// Find smallest distance betwen two medoids from different clusters
	// Find biggest density
	// Divide one by another - there you go

	unsigned int k = devCalculateClustersAndDensities( clusters, densities );
	unsigned int j = 0, i = 0;

	float prevDistance = 0;
	float currDistnace = 0;

	for ( i = 0; i < k; i++ ) {
		for ( j = 0; j < k; j++ ) {
			if ( i != j && clusters[ i] != clusters[ j] ) {
				currDistnace = distance( dPopulationPool[ threadIdx.x].medoids[ i],
					dPopulationPool[ threadIdx.x].medoids[ j] );

				if ( prevDistance > currDistnace || prevDistance == 0 ) {
					prevDistance = currDistnace;
				}
			}
		} // j
	}// i

	lowestDensity = prevDistance;
	prevDistance = 0;
	for ( i = 0; i < k; i++ ) {
		if ( prevDistance < densities[ i] || prevDistance == 0 ) {
			prevDistance = densities[ i];
		}
	}

	lowestDensity /= prevDistance;

	indexes[ threadIdx.x] = lowestDensity;
}

//====================================================================
ErrorCode calculateDI( preciseResult & topDi ) {
	enableErrorControl;
	float * dResults;
	cudaMalloc( &dResults, populationSize * sizeof(float) );

	kernelCalculateDI<<< 1, populationSize >>>( dResults );
	cutilDeviceSynchronize();
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] After kernelDI - %s\n", cudaGetErrorString( cuErr ));
	}

	float * hResults = (float*)malloc( populationSize * sizeof(float) );
	cudaMemcpy( hResults, dResults, populationSize * sizeof(float), cudaMemcpyDeviceToHost );
	unsigned int res = 0;

	//printf( " \n DI index: \n" );
	topDi.min = topDi.max = hResults[ 0];
	topDi.sum = 0.0;
	for ( int i = 1; i < populationSize; i++ ) {
		if ( topDi.min == 0 || hResults[ i] < topDi.min ) {
			topDi.min = hResults[ i];
		}
		if ( topDi.max == 0 || hResults[ i] > topDi.max ) {
			topDi.max = hResults[ i];
		}
		topDi.sum += hResults[ i];
	}

	return errOk;
}
//====================================================================

// <<< 1, dPopulationSize >>>
__global__ void kernelCalculateRand( unsigned int * preclasified, float * rand) {
	unsigned char clusters[ MEDOID_VECTOR_SIZE];
	float densities[ MEDOID_VECTOR_SIZE];
	
	unsigned int k = devCalculateClustersAndDensities( clusters, densities );

	unsigned int i,j;
	unsigned int t = 0, f = 0;
	bool we, they;

	for ( i = 0; i < dNumEntries; i++ ) {
		for ( j = 0; j < dNumEntries; j++ ) {
			if ( i != j ) {
				we = ( clusters[ dMembership[ threadIdx.x * dNumEntries + i]] ==
				clusters[ dMembership[ threadIdx.x * dNumEntries + j]] );

				they = ( preclasified[ i] == preclasified[ j] );

				if ( we == they ) {
					t++;
				} else {
					f++;
				}
			}
		} // for j
	} // for i

	rand[ threadIdx.x] = (float)( t ) / (float)( t + f );
}
//====================================================================

ErrorCode calculateRand( DataStore * dataStore, preciseResult & topRand ) {
	// Get Clasified data
	static float bestRand = 0;
	enableErrorControl;

	unsigned int * dPreclasified;
	cudaMalloc( &dPreclasified, hNumEntries * sizeof(unsigned int) );
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] rand malloc - %s\n", cudaGetErrorString( cuErr ));
	}
	cudaMemcpy( dPreclasified, dataStore->classes, hNumEntries * sizeof(unsigned int), cudaMemcpyHostToDevice );
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] Rand memcpy - %s\n", cudaGetErrorString( cuErr ));
	}

	float * dRandIndex;
	cudaMalloc( &dRandIndex, populationSize * sizeof(float) );
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] Rand malloc2 - %s\n", cudaGetErrorString( cuErr ));
	}
	
	// run comparision
	kernelCalculateRand<<< 1, populationSize >>>( dPreclasified, dRandIndex );
	cutilDeviceSynchronize();
	cuErr = cudaGetLastError();
	if ( cuErr != cudaSuccess ) {
		printf( "[E][cuda] After kernelRand - %s\n", cudaGetErrorString( cuErr ));
	}

	// display results
	float * hRandIndex = (float*)malloc( populationSize * sizeof(float) );
	cudaMemcpy( hRandIndex, dRandIndex, populationSize * sizeof(float), cudaMemcpyDeviceToHost );

	for ( int i = 1; i < populationSize; i++ ) {
		if ( topRand.min == 0 || hRandIndex[ i] < topRand.min ) {
			topRand.min = hRandIndex[ i];
		}
		if ( topRand.max == 0 || hRandIndex[ i] > topRand.max ) {
			topRand.max = hRandIndex[ i];
		}
		topRand.sum += hRandIndex[ i];
	}

	return errOk;
}
//====================================================================

void CleanResults( preciseResult & results ) {
	results.max = 0.0;
	results.mean = 0.0;
	results.min = 0.0;
	results.sum = 0.0;
}
//====================================================================

void CleanAlgResults( algResults & results ) {
	CleanResults( results.bdi );
	CleanResults( results.clusters );
	CleanResults( results.di );
	CleanResults( results.rand );
	CleanResults( results.time );
}
//====================================================================