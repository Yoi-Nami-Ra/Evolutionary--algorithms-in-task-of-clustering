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


#define MAX_CLUSTER_SIZE 2
#define MEDOID_VECTOR_SIZE 6
#define MAX_NEIGHBOURS 10

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
	unsigned int medoids[ MEDOID_VECTOR_SIZE];
	unsigned int clusters[ MEDOID_VECTOR_SIZE];
	unsigned int clusterMembership[ MEDOID_VECTOR_SIZE]; ///< cluster this medoid belongs to
} PopMember;

/**
 * To hold solution results for each poppulation member.
 */
typedef struct {
	float densities; ///< sum of all densities
	unsigned int * recordMembership;
	unsigned int numOfClusters; ///< How many clusters there are
	float clusterDensities[ MEDOID_VECTOR_SIZE]; ///< densities for specific medoids
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
	//unsigned int maxClusterSize;
	//unsigned int medoidsVectorSize;
	unsigned int crosFactor;
	PopMember * population;
	Solution * solutions;
	char * dominanceMatrix; ///< Describes which solution dominates which
	unsigned int * dominanceCounts; ///< 
	unsigned int blocksPerEntries;
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
	char crossTemplate[MEDOID_VECTOR_SIZE];
} BreedingTable;
//==============================================
//== Functions

/**
 * Just some defaults.
 */
ErrorCode GenerateDefaultProps( EvolutionProps * props );

/**
 * Starts the main process of eunning evolutional algorithms.
 * @param props	- [in, out] structure to hold required data and to return results.
 * 
 * @return		- error code if any.
 */
ErrorCode RunClustering( EvolutionProps * props );

#endif // CLUSTERING_H