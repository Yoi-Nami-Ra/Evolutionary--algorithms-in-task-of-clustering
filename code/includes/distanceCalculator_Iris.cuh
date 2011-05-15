// Header file for dataLoader

#ifndef DISTANCECALCULATOR_IRIS_CUH
#define DISTANCECALCULATOR_IRIS_CUH

#include "globals.cuh"

// Data types

// Functions
/*
 * Run Calculations.
 */
ErrorCode startCalculatingDistances();

/*
 * Get Array with calculated distances.
 */
const float* getDistances();

/*
 * Frees up memory
 */
ErrorCode releaseDistances();

/*
 * Gives array of neighbours
 */
const unsigned int* getNeighbours();

/*
 * release memory held by neighbours
 */
ErrorCode releaseNeighbours();

#endif //DISTANCECALCULATOR_IRIS_CUH