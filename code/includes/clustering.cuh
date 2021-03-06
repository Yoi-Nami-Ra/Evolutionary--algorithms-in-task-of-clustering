#include "globals.cuh"
#include "errors.cuh"
#include "dataLoader.cuh"

// =======================
// Types

#define MAX_CLUSTER_SIZE 1
#define MEDOID_VECTOR_SIZE 48
#define kMaxNeighboursToUSe 29

#define CROS_FACTOR 3
#define OBJECTIVES 4

typedef struct unitAttributes {
	unsigned int numNeighbours;
	unsigned int clusterMaxSize;
} unitAttributes;

typedef struct unit {
	unitAttributes attr;
	unsigned int medoids[ MEDOID_VECTOR_SIZE];
	char clusters[ MEDOID_VECTOR_SIZE];
} unit;

typedef struct index2d {
	unsigned int x;
	unsigned int y;
} index2d;

typedef struct breedDescriptor {
	unsigned int parent1;
	unsigned int parent2;
	unsigned int child;
	unsigned int factor;
} breedDescriptor;

typedef struct {
	float min;
	float max;
	float mean;
	float sum;
} preciseResult;

typedef struct {
	preciseResult rand;
	preciseResult bdi;
	preciseResult di;
	preciseResult clusters;
	preciseResult time;
} algResults;

// =======================
// Functions

/*
 * Generates Random Starting population
 */
ErrorCode generateRandomPopulation( unsigned int popSize );

/*
 * start and run clustering algorithm
 */
ErrorCode runClustering( unsigned int popSize, unsigned int steps, DataStore * dataStore, algResults * results );

/*
 * Calculated Davies-Bouldin index for all solutions left in Cuda memory.
 */
ErrorCode calculateBDI( preciseResult & topBDI, preciseResult & clusters );

ErrorCode calculateDI( preciseResult & topDi );

ErrorCode calculateRand( DataStore * dataStore, preciseResult & topRand ) ;

void CleanResults( preciseResult & results );

void CleanAlgResults( algResults & results );