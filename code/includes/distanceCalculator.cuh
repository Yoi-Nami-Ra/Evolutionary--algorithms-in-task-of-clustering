// Header file for distanceCalculator

#ifndef DISTANCECALCULATOR_CUH
#define DISTANCECALCULATOR_CUH


#include "globals.cuh"
#include "errors.cuh"

#ifdef DATABASE_IRIS
#include "distanceCalculator_Iris.cuh"
#endif // DATABASE_IRIS


// Functions

/*
 * Calculates distances betwen records of current database source
 */
ErrorCode CalculateDistances();


#endif //DISTANCECALCULATOR_CUH