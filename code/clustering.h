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


#define MAX_CLUSTER_SIZE 10
#define MEDOID_VECTOR_SIZE 20
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
	char clusters[ MEDOID_VECTOR_SIZE];
} PopMember;

/**
 * To hold solution results for each poppulation member.
 */
typedef struct {
	float densities; ///< sum of all densities
	unsigned int clusterMembership[ MEDOID_VECTOR_SIZE]; ///< cluster this medoid belongs to
	unsigned int * recordMembership;
	float connectivity;
	float disconnectivity;
	unsigned int errors; ///< how much errors has been found
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
	unsigned int blocksPerEntries;
} EvolutionProps;
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