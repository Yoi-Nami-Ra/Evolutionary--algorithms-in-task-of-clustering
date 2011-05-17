


//====================================================================
//== Includes
#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
//#include <math.h>
#include "globals.cuh"
#include "errors.cuh"
#include "clustering.cuh"
#include "distanceCalculator.cuh"
#include "dataLoader.cuh"
#include <cutil_inline.h>
#include <cutil_math.h>
#include <cuda.h>
#include <time.h>

//====================================================================
//== Types

//====================================================================
//== Constants and Globals
static float * dDistances;
static float * dNeighbours;
static unit * hPopulationPool;
static unit * dPopulationPool;
static unsigned int populationSize;

const char threadsPerBlock = 50;

//====================================================================
//== Declarations
texture<float, cudaTextureType1D, cudaReadModeElementType> texRefDistances;
texture<float, cudaTextureType1D, cudaReadModeElementType> texRefNeighbour;

// host
ErrorCode runAlgorithms( unsigned int steps );
// device
__global__ void randomPopulation( unsigned int popSize, unsigned char blockSize, unit * population, unsigned int numEntries );
__global__ void kernelMembershipAndDensity( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution );
__device__ float distance( unsigned int a, unsigned int b );
__device__ uint distanceIdx(uint x, uint y);
__global__ void kernelConnectivity( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution );
__global__ void kernelDisconnectivity( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution );
__global__ void kernelCorectness( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution );
__global__ void kernelSorting( float * fitnesResults, bool * dominanceMatrix, 
	unsigned int blocksPerSolution, unsigned int populationSize );
__global__ void kernelDominanceCount( bool * dominanceMatrix, unsigned int * dominanceCounts, unsigned int popSize );

//====================================================================
//== Functions

ErrorCode generateRandomPopulation( unsigned int popSize ) {

	populationSize = popSize;
	hPopulationPool = (unit*)malloc( popSize * sizeof(unit) );
	srand( time( 0 ));

	for ( int k = 0; k < popSize; k++ ) {
		// attributes
		hPopulationPool[k].attr.clusterMaxSize = rand() % MAX_CLUSTER_SIZE + 1;
		hPopulationPool[k].attr.numNeighbours = rand() % MAX_NEIGHBORS + 1;
		unsigned int clustersSum = MEDOID_VECTOR_SIZE;
		unsigned int proposal;
		bool proposalOk = false;
		for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
			do {
				proposalOk = true;
				proposal = rand() % numEntries();
				for ( int j = 0; j < i; j++ ) {
					if ( proposal == hPopulationPool[k].medoids[ j] ) {
						proposalOk = false;
						break;
					}
				}
			} while ( !proposalOk );
			hPopulationPool[k].medoids[ i] = proposal;

			if ( clustersSum > 0 ) {
				proposal = rand() % hPopulationPool[k].attr.clusterMaxSize + 1;
				if ( clustersSum < proposal ) {
					proposal += clustersSum - proposal;
				}
				clustersSum -= proposal;
				hPopulationPool[k].clusters[ i] = proposal;
			} else {
				hPopulationPool[k].clusters[ i] = 0;
			}
		} // for each medoid in vector
	} // for each member of population

	cudaMalloc( &dPopulationPool, popSize * sizeof(unit) );
	cudaMemcpy( dPopulationPool, hPopulationPool, popSize * sizeof(unit), cudaMemcpyHostToDevice );

	return errOk;
}
//====================================================================

ErrorCode runClustering( unsigned int popSize, unsigned int steps ) {
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
	
	unsigned int distancesSize = numEntries() * ( numEntries() -1 ) /2 * sizeof(float);
	cudaMalloc( &dDistances, distancesSize );
	cudaMemcpy( dDistances, getDistances(), distancesSize, cudaMemcpyHostToDevice );
	//   bind distances to texture
	cudaBindTexture( &offset, &texRefDistances, dDistances, &channelDesc, distancesSize );

	//   Allocate memory for distances	
	unsigned int neighbourSize = numEntries() * MAX_NEIGHBORS * sizeof(unsigned int);
	cudaMalloc( &dNeighbours, neighbourSize );
	cudaMemcpy( dNeighbours, getNeighbours(), neighbourSize, cudaMemcpyHostToDevice );
	//   bind distances to texture
	cudaBindTexture( &offset, &texRefNeighbour, dNeighbours, &channelDesc, neighbourSize );

	// generate startup population
	generateRandomPopulation( popSize );

	runAlgorithms( steps );

	return errOk;
}
//====================================================================

__global__ void randomPopulation( unsigned int popSize, unsigned char blockSize, unit * population, unsigned int numEntries ) {
	unsigned int index = blockIdx.x * blockSize + threadIdx.x;

	
}
//====================================================================

ErrorCode runAlgorithms( unsigned int steps ) {
	// Set ups
	float * hFitnesResults = (float*)malloc( populationSize * 4 * sizeof(float) );
	float * dFitnesResults = 0;
	char * dMembership = 0;
	unsigned int blocks = numEntries() / threadsPerBlock;
	while ( blocks * threadsPerBlock < numEntries()) {
		blocks++;
	}

	// populationSize x numObjectives x blocksPerSolution
	cudaMalloc( &dFitnesResults, populationSize * 4 * blocks * sizeof(float) );
	cudaMalloc( &dMembership, populationSize * numEntries() * sizeof(char) );

	// population dominations
	bool * dDominanceMatrix = 0;
	unsigned int * dDominanceCounts = 0;
	cudaMalloc( &dDominanceMatrix, populationSize * populationSize * sizeof(bool) );
	cudaMalloc( &dDominanceCounts, populationSize * sizeof(unsigned int) );

	dim3 dimGrid( blocks, populationSize );
	dim3 dimBlock( threadsPerBlock );

	dim3 dimGrid2( populationSize );
	dim3 dimBlock2( MEDOID_VECTOR_SIZE );

	unsigned int solutionsLeft = 0;
	bool * solutionsSelected = (bool*)malloc( populationSize * sizeof(bool) );
	unsigned int * solutionFronts = (unsigned int*)malloc( populationSize * ( populationSize + 1 ) * sizeof(unsigned int));
	unsigned int currFront;
	unsigned int currFrontSize;

	bool * hDominanceMatrix = (bool*)malloc( populationSize * populationSize * sizeof(bool));
	unsigned int * hDominanceCounts = (unsigned int*)malloc( populationSize * sizeof(unsigned int));

	for (int i = 0; i < steps; i++ ) {
		// membership and density phase
		kernelMembershipAndDensity<<<dimGrid, dimBlock>>>( dFitnesResults, dMembership, threadsPerBlock, dPopulationPool, numEntries(), blocks );
		cutilDeviceSynchronize();

		// connectivity phase
		kernelConnectivity<<<dimGrid, dimBlock>>>( dFitnesResults, dMembership, threadsPerBlock, dPopulationPool, numEntries(), blocks );
		cutilDeviceSynchronize();

		// sum up results for density and Connectivity
		// TODO:

		// disconnectivity phase
		kernelDisconnectivity<<<dimGrid2, dimBlock2>>>( dFitnesResults, dMembership, MEDOID_VECTOR_SIZE, dPopulationPool, numEntries(), populationSize );
		cutilDeviceSynchronize();

		// correctness phase
		kernelCorectness<<<dimGrid2, dimBlock2>>>( dFitnesResults, dMembership, MEDOID_VECTOR_SIZE, dPopulationPool, numEntries(), populationSize );
		cutilDeviceSynchronize();

		// sorting
		kernelSorting<<<populationSize, populationSize>>>( dFitnesResults, dDominanceMatrix, blocks, populationSize );
		cutilDeviceSynchronize();

		kernelDominanceCount<<<1, populationSize>>>( dDominanceMatrix, dDominanceCounts, populationSize );
		cutilDeviceSynchronize();

		cudaMemcpy( hDominanceMatrix, dDominanceMatrix, populationSize * populationSize * sizeof(bool), cudaMemcpyDeviceToHost );
		cudaMemcpy( hDominanceCounts, dDominanceCounts, populationSize * sizeof( unsigned int ), cudaMemcpyDeviceToHost );

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
			for ( j = 0; j < populationSize; j++ ) {
				if ( !solutionsSelected[ j] && hDominanceCounts[ j] == 0 ) {
					solutionFronts[ currFront * populationSize + (++currFrontSize)] = j;
					solutionsSelected[ j] = true;
					solutionsLeft--;
				}
			}
			solutionFronts[ currFront * populationSize] = currFrontSize;
			solutionsLeft -= currFrontSize;
			
			if ( solutionsLeft > 0 ) {
				// for each solution dominated by solution from this front - reduce domination count
				for ( j = 0; j < currFrontSize; j++ ) {
					for ( int k = 0; k < populationSize; k++ ) {
						if ( hDominanceMatrix[ solutionFronts[ currFront * populationSize + j + 1] * populationSize + k] ) {
							hDominanceCounts[ k] -= 1;
						}
					}
				}
			}

			// now for next front
			currFront++;
		}

		// selection
		solutionLeft = populationSize / 2; // select half size of population
		for ( j = 0; j < populationSize; j++ ) {
			solutionsSelected[ j] = false;
		}

		currFront = 0;
		while ( solutionLeft > 0 ) {
			// if we need more than the current front can offer
			if ( solutionLeft > solutionFronts[ currFront * populationSize] {
				for ( j = 0; j < solutionLeft > solutionFronts[ currFront * populationSize]; j++ ) {
					ssolutionSelected[ solutionFronts[ currFront * populationSize + j * 1]] = true;
				}
			} else {
				// this fron has more than we need
			}
		}
	}

	return errOk;
}
//====================================================================

__global__ void kernelMembershipAndDensity( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution ) {
	unsigned int solution = blockIdx.y;
	unsigned int rekord = blockIdx.x * threadsPerBlock + threadIdx.x;

	__shared__ unit thisSolution;
	__shared__ float density[ 256]; // shared table to hold density results for futher calculation

	if ( threadIdx.x == 0 ) {
		// first thread should load solution for this block to shared memory
		thisSolution = populationPool[ solution];
	}

	// Sync up threads
	__syncthreads();

	float prevDistance = distance( rekord, thisSolution.medoids[ 0] );
	float currDistance;
	unsigned int res = 0;

	for ( int i = 1; i < MEDOID_VECTOR_SIZE; i++ ) {
		currDistance = distance ( rekord, thisSolution.medoids[ 0] );
		if ( currDistance > prevDistance ) {
			prevDistance = currDistance;
			res = i;
		}
	}

	density[ threadIdx.x] = prevDistance;

	membership[ solution * numEntries + rekord] = res;

	// Sync up threads
	__syncthreads();	
	if ( threadIdx.x == 0 ) {
		// sum solutions from all threads in this block
		currDistance = 0;
		for ( int i = 0; i < threadsPerBlock; i++ ) {
			currDistance += density[ i];
		}

		// sum all solutions for this block
		fitnesResults[ solution * 4 * blocksPerSolution + 0 * blocksPerSolution + blockIdx.x] = currDistance;
	}

	// Sync once more
	__syncthreads();

	// now sum all results for this solution (from each block)
	if ( blockIdx.x == 0 && threadIdx.x == 0 ) {
		currDistance = 0;
		for ( int i = 0; i < blocksPerSolution; i++ ) {
			currDistance += fitnesResults[ solution * 4 * blocksPerSolution + 0 * blocksPerSolution + i];
		}
		fitnesResults[ solution * 4 * blocksPerSolution + 0 * blocksPerSolution + 0] = currDistance;
	}
}
//====================================================================

__device__ float distance( unsigned int a, unsigned int b ) {
	return tex1Dfetch( texRefDistances, distanceIdx( a, b ));
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
	return tex1Dfetch( texRefNeighbour, record * MAX_NEIGHBORS + num );
}
//====================================================================

__global__ void kernelConnectivity( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution ) {
	
	unsigned int solution = blockIdx.y;
	unsigned int record = blockIdx.x * threadsPerBlock + threadIdx.x;
	unit thisSolution = populationPool[ solution];
	unsigned int memberOf = membership[ solution * numEntries + record];
	unsigned int numOfNeighbours = thisSolution.attr.numNeighbours;

	__shared__ float blockResults[256];
	float result = 0;

	for ( int i = 0; i < numOfNeighbours; i++ ) {
		if ( memberOf == membership[ solution * numEntries + i] ) {
			result += 1.0f / (float)numOfNeighbours;
		}
	}

	blockResults[ threadIdx.x] = result;

	// Sync and sum solutions for this block
	__syncthreads();
	if ( threadIdx.x == 0 ) {
		// sum solutions from all threads in this block
		result = 0;
		for ( int i = 0; i < threadsPerBlock; i++ ) {
			result += blockResults[ i];
		}

		// sum all solutions for this block
		fitnesResults[ solution * 4 * blocksPerSolution + 1 * blocksPerSolution + blockIdx.x] = result;
	}

	// Sync once more
	__syncthreads();

	// now sum all results for this solution (from each block)
	if ( blockIdx.x == 0 && threadIdx.x == 0 ) {
		result = 0;
		for ( int i = 0; i < blocksPerSolution; i++ ) {
			result += fitnesResults[ solution * 4 * blocksPerSolution + 1 * blocksPerSolution + i];
		}
		fitnesResults[ solution * 4 * blocksPerSolution + 1 * blocksPerSolution + 0] = result;
	}

}
//====================================================================

__global__ void kernelDisconnectivity( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution ) {

	__shared__ unsigned int clusters[ MEDOID_VECTOR_SIZE];
	__shared__ unsigned int comparisions[ MEDOID_VECTOR_SIZE];
	__shared__ unsigned int counts[ MEDOID_VECTOR_SIZE];


	// label medoids to their clusters
	if ( threadIdx.x == 0 ) {
		for ( int i=0, j=0; i<MEDOID_VECTOR_SIZE && j<MEDOID_VECTOR_SIZE; i++ ) {
			for (int k=0; k<populationPool[ blockIdx.y].clusters[i]; k++,j++) {
				clusters[ j] = i;
			}
		}
	}

	__syncthreads();

	float currDistance  = 0;
	comparisions[ threadIdx.x] = 0;
	counts[ threadIdx.x] = 0;
	// For each medoid in the vector
	for ( unsigned int i=0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if ( i == threadIdx.x || clusters[ threadIdx.x] == clusters[ i] ) {
			// if medoid the same or same cluster - skip
			continue;
		}
		comparisions[ threadIdx.x]++;
		currDistance = distance( populationPool[ blockIdx.y].medoids[ threadIdx.x], populationPool[ blockIdx.y].medoids[ i] );
		// now find MND for this medoid
		for ( unsigned int j = 0; j < MEDOID_VECTOR_SIZE; j++ ) {
			if ( j == i || j == threadIdx.x ) {
				continue;
			}
			// counts if someone else is closer
			if ( distance( populationPool[ blockIdx.y].medoids[ threadIdx.x], populationPool[ blockIdx.y].medoids[ j] ) < currDistance ) {
				counts[ threadIdx.x]++;
			}
			if ( distance( i, j ) < currDistance ) {
				counts[ threadIdx.x]++;
			}
		}
	}

	__syncthreads();

	if ( threadIdx.x == 0 ) {
		float compars = 0;
		float count = 0;
		for ( int i = 0; i < MEDOID_VECTOR_SIZE; i++ ) {
			compars += comparisions[ i];
			count += counts[ i];
		}
		fitnesResults[ blockIdx.x * 4 * blocksPerSolution + 2 * blocksPerSolution + 0] = count/compars;
	}
}
//====================================================================

__global__ void kernelCorectness( float * fitnesResults, char * membership,
	unsigned int threadsPerBlock, unit * populationPool, unsigned int numEntries, unsigned int blocksPerSolution ) {

	__shared__ char checks[ MEDOID_VECTOR_SIZE];
	__shared__ unsigned int medoids[ MEDOID_VECTOR_SIZE];

	if ( threadIdx.x == 0 ) {
		memcpy( medoids, populationPool[ blockIdx.y].medoids, MEDOID_VECTOR_SIZE * sizeof(unsigned int) );
	}

	__syncthreads();

	unsigned int thisMedoid = medoids[ threadIdx.x];
	char count = 0;
	
	for ( int i=0; i < MEDOID_VECTOR_SIZE; i++ ) {
		if (i == threadIdx.x) {
			continue;
		}

		if ( thisMedoid == medoids[ i] ) {
			count++;
		}
	}

	checks[ threadIdx.x] = count;

	__syncthreads();

	if ( threadIdx.x == 0 )  {
		for ( int i = 0; i < threadsPerBlock; i++ ) {
			count += checks[ i];
		}
		fitnesResults[ blockIdx.x * 4 * blocksPerSolution + 3 * blocksPerSolution + 0] = count;
	}
}
//====================================================================

// <<< populationSize, populationSize >>>
__global__ void kernelSorting( float * fitnesResults, bool * dominanceMatrix, 
	unsigned int blocksPerSolution, unsigned int populationSize ) {

	__shared__ bool dominating[ MAX_POPULATION_SIZE];
	__shared__ float thisSolutionFitnesResults[4];

	if ( threadIdx.x == 0 ) {
		thisSolutionFitnesResults[ 0] = fitnesResults[ blockIdx.x * 4 * blocksPerSolution + 0 * blocksPerSolution + 0];
		thisSolutionFitnesResults[ 1] = fitnesResults[ blockIdx.x * 4 * blocksPerSolution + 1 * blocksPerSolution + 0];
		thisSolutionFitnesResults[ 2] = fitnesResults[ blockIdx.x * 4 * blocksPerSolution + 2 * blocksPerSolution + 0];
		thisSolutionFitnesResults[ 3] = fitnesResults[ blockIdx.x * 4 * blocksPerSolution + 3 * blocksPerSolution + 0];
	}

	__syncthreads();

	dominating[ threadIdx.x] = true;
	for ( int i = 0; i < 4 ;i++ ) {
		if ( fitnesResults[ threadIdx.x * 4 * blocksPerSolution + i * blocksPerSolution + 0] >=
			thisSolutionFitnesResults[ 0] ) {
				dominating[ threadIdx.x] = false;
		}
	}

	__syncthreads();

	if ( threadIdx.x == 0 ) {
		memcpy( &dominanceMatrix[ blockIdx.x * populationSize], dominating, populationSize * sizeof(float) );
	}
}
//====================================================================

// <<< 1, ppulationSize >>>
__global__ void kernelDominanceCount( bool * dominanceMatrix, unsigned int * dominanceCounts, unsigned int popSize ) {

	unsigned int count = 0;

	for ( int i = 0; i < popSize; i++ ) {
		count += dominanceMatrix[ i * popSize + threadIdx.x];
	}

	dominanceCounts[ threadIdx.x] = count;
}
//====================================================================

// <<< 1, (numSolutions, kryterions) >>>
__global__ void kernelFrontDensity ( unsigned int * front, unsigned int frontSize, unsigned int blocksPerSolution
	unsigned int toBeSelected,  ) {

	__shared__ float frontDensities [ MAX_POPULATION_SIZE * 4];

	unsigned int lesser;
	bool lesserFound = false;
	float lesserResult;
	unsigned int bigger;
	bool biggerFound = false;
	float biggerResult;

	float thisResult = fitnesResults[ front[ threadIdx.x] * 4 * blocksPerSolution + threadIdx.y * blocksPerSolution + 0];
	float currResult;

	for ( int i = 0; i < frontSize; i++ ) {
		if ( threadIdx.x == i ) {
			// skip if same
			continue;
		}

		currResult = fitnesResults[ front[ i] * 4 * blocksPerSolution + threadIdx.y * blocksPerSolution + 0];
		// check if lesser
		if ( thisResult > currResult ) {
			if ( !lesserFound ) {
				lesser = i;
				lesserFound = true;
				lesserResult = currResult;
			} else {
				if ( lesserResult < currResult ) {
					lesser = i;
					lesserFound = true;
					lesserResult = currResult;
				}
			}
		}

		// check if bigger
		if ( thisResult < currResult ) {
			if ( !biggerFound ) {
				bigger = i;
				biggerFound = true;
				biggerResult = currResult;
			} else {
				if ( biggerResult > currResult ) {
					bigger = i;					
					biggerResult = currResult;
				}
			}
		}
	} // for each solution in this front

	// is this edge solution ?
	if ( !lesserFound || !biggerFound ) {
		frontDensities[ frontSize * threadIdx.x + threadIdx.y] = -1;  
	} else {
		frontDensities[ frontSize * threadIdx.x + threadIdx.y] = biggerResult - lesserResult;
	}

}