#include "globals.cuh"
#include "errors.cuh"

// =======================
// Types

#define MAX_NEIGHBORS 10
#define MAX_CLUSTER_SIZE 10
#define MEDOID_VECTOR_SIZE 20

typedef struct unitAttributes {
	unsigned int numNeighbours;
	unsigned int clusterMaxSize;
} unitAttributes;

typedef struct unit {
	unitAttributes attr;
	unsigned int medoids[ MEDOID_VECTOR_SIZE];
	char clusters[ MEDOID_VECTOR_SIZE];
} unit;

// =======================
// Functions

/*
 * 
 */
generateRandomPopulation();
