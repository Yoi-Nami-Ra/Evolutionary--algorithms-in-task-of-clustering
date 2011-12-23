/**
 17/11/2011
 Jaroslaw Wojtasik

 noCuda

 clustering.h
 **/

#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "errors.h"
#include "dataLoader.h"

//==============================================
//== Globals


//==============================================
//== Types

/**
 * Additional attributes to describe population member.
 */
typedef struct {
	unsigned int numNeighbours;
	unsigned int clusterMaxSize;
} PopMemberAttributes;

/**
 * Dexcription of single population member.
 */
typedef struct {
	PopMemberAttributes attr;
	unsigned int * medoids; ///< list of medoids describing clusters [ medoid size]
	unsigned int * clusters; ///< list of cluster sizes [medoid size]
	unsigned int * clusterMembership; ///< cluster this medoid belongs to
} PopMember;

typedef struct {
    float min;
    float max;
    float mean;
    float sum;
    unsigned int count;
} Results;

/**
 * To hold solution results for each poppulation member.
 */
typedef struct {
	float densities; ///< sum of all densities
	unsigned int * recordMembership;
	unsigned int numOfClusters; ///< How many clusters there are
	float * clusterDensities; ///< densities for specific medoids
	float connectivity;
	float disconnectivity;
	float errors; ///< how much errors has been found
	float resultBDI; ///< Here result of BD index to be placed
	float resultDI; ///< Here result of D index to be placed
	float resultRand; ///< Here result of Rand index to be placed
} Solution;

/**
 * Holds data required for algorithms to run.
 */
typedef struct {
	unsigned int popSize;
	unsigned int evoSteps;
	DataStore *dataStore;
	unsigned int maxClusterSize;
	unsigned int medoidsVectorSize;
    unsigned int maxNeighbours;
	unsigned int crosFactor;
	PopMember * population;
	Solution * solutions;
	char * dominanceMatrix; ///< Describes which solution dominates which
	unsigned int * dominanceCounts; ///< 
	unsigned int blocksPerEntries;
    Results resultBDI; ///< Here result of BD index to be placed
	Results resultDI; ///< Here result of D index to be placed
	Results resultRand; ///< Here result of Rand index to be placed
} EvolutionProps;

/**
 * Data required to calculate
 */
typedef struct {
	EvolutionProps * props;
	unsigned int * front;
	unsigned int frontSize;
	float * densities;
} FrontDensities;

/**
 * Describes single child.
 */
typedef struct {
	unsigned int parent1;
	unsigned int parent2;
	unsigned int child;
	unsigned int factor;
} BreedDescriptor;

/**
 * Breeding properties for crossing function.
 */
typedef struct {
	BreedDescriptor * table;
	EvolutionProps * props;
} BreedingTable;
//==============================================
//== Functions

ErrorCode ConfigureAlgorithms(EvolutionProps * props);

/**
 * Just some defaults.
 */
ErrorCode DefaultProps(EvolutionProps * props, DataStore * dataStore);

/**
 * Starts the main process of eunning evolutional algorithms.
 * @param props	- [in, out] structure to hold required data and to return results.
 * 
 * @return		- error code if any.
 */
ErrorCode RunClustering( EvolutionProps * props );

#endif // CLUSTERING_H