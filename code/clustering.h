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
 * Holds data required for algorithms to run.
 */
typedef struct {
	unsigned int popSize;
	unsigned int evoSteps;
	DataStore *dataStore;
	unsigned int maxClusterSize;
	unsigned int medoidsVectorSize;
	unsigned int crosFactor;
	PopMember * population;
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