


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

//==============================================
//== Declarations
texture<float, cudaTextureType1D, cudaReadModeElementType> texRefDistances;
texture<float, cudaTextureType1D, cudaReadModeElementType> texRefNeighbour;

static const char blockSize = 16;

__global__ void randomPopulation( unsigned int popSize, unsigned char blockSize, unit * population, unsigned int numEntries );

//==============================================
//== Functions

ErrorCode generateRandomPopulation( unsigned int popSize ) {

	hPopulationPool = (unit*)malloc( popSize * sizeof(unit) );
	srand( time( 0 ));

	for ( int k = 0; k < popSize; k++ ) {
		// attributes
		hPopulationPool[k];
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
				if (( clustersSum - proposal ) < 0 ) {
					proposal += clustersSum - proposal;
				}
				clustersSum -= proposal;
				hPopulationPool[k].clusters[ i] = proposal;
			} else {
				hPopulationPool[k].clusters[ i] = 0;
			}
		} // for each medoid in vector
	} // for each member of population

	return errOk;
}

ErrorCode runClustering( unsigned int popSize ) {
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
	cudaBindTexture( &offset, &texRefDistances, dNeighbours, &channelDesc, neighbourSize );

	// generate startup population
	generateRandomPopulation( popSize );

	return errOk;
}
//==============================================

__global__ void randomPopulation( unsigned int popSize, unsigned char blockSize, unit * population, unsigned int numEntries ) {
	unsigned int index = blockIdx.x * blockSize + threadIdx.x;

	
}
//==============================================