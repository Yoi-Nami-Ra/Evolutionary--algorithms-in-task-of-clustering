// Header file for dataLoader

#ifndef DISTANCECALCULATOR_IRIS_CUH
#define DISTANCECALCULATOR_IRIS_CUH

#include "globals.cuh"

// Data types

// Functions
/*
 * Run Calculations.
 */
ErrorCode StartCalculatingDistances();

/*
 * Get Array with calculated distances.
 */
float* GetDistances();

/*
 * Frees up memory
 */
ErrorCode ReleaseDistances();

#endif //DISTANCECALCULATOR_IRIS_CUH