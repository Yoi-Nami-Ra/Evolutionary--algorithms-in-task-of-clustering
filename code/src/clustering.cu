


//==============================================
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

//==============================================
//== Types

//==============================================
//== Constants and Globals
static float * dDistances;
static float * dNeighbours;
static unit * hPopulationPool;
static unit * dPopulationPool;
static unsigned int populationSize;

const char threadsPerBlock = 50;

//==============================================
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

//==============================================
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
//==============================================

__global__ void randomPopulation( unsigned int popSize, unsigned char blockSize, unit * population, unsigned int numEntries ) {
	unsigned int index = blockIdx.x * blockSize + threadIdx.x;

	
}
//==============================================

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

	dim3 dimGrid( blocks, populationSize );
	dim3 dimBlock( threadsPerBlock );

	for (int i = 0; i < steps; i++ ) {
		// membership and density phase
		kernelMembershipAndDensity<<<dimGrid, dimBlock>>>( dFitnesResults, dMembership, threadsPerBlock, dPopulationPool, numEntries(), blocks );
		cutilDeviceSynchronize();

		// connectivity phase
		kernelConnectivity<<<dimGrid, dimBlock>>>( dFitnesResults, dMembership, threadsPerBlock, dPopulationPool, numEntries(), blocks );
		cutilDeviceSynchronize();

		// disconnectivity phase
		// correctness phase
	}

	return errOk;
}
//==============================================

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
//==============================================

__device__ float distance( unsigned int a, unsigned int b ) {
	return tex1Dfetch( texRefDistances, distanceIdx( a, b ));
}
//==============================================

__device__ uint distanceIdx(uint x, uint y) {
	if ( y > x ) {
		return y * (y - 1) / 2 + x;
	} else {		
		return x * (x - 1) / 2 + y;
	}
}
//==============================================

__device__ unsigned int neighbour( unsigned int record, unsigned int num ) {
	return tex1Dfetch( texRefNeighbour, record * MAX_NEIGHBORS + num );
}
//==============================================

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
//==============================================