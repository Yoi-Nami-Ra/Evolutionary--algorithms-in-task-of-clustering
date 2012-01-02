// Header file for distanceCalculator

#ifndef DISTANCECALCULATOR_CUH
#define DISTANCECALCULATOR_CUH


#include "globals.cuh"
#include "errors.cuh"
#include "dataLoader.cuh"

//==============================================
//== Globals

#define kMaxNeighbours 30


//==============================================
//== Types

//==============================================
//== Functions

/**
 * Calculates distances from a selected Loader.
 * @param num		- [in] index of selected loader.
 * @param dataStore	- [out] pointer to where buffer of distances should be stored.
 *
 * @return Error code if any
 */
ErrorCode GetCalculatedDistances( unsigned int num, DataStore * dataStore );

/**
 * Calculates proper index in the distances vector.
 * @param a	- [in] column or row
 * @param b - [in] column or row
 *
 * @return calculated index
 */
unsigned int DistanceVIdx( unsigned int a, unsigned int b );


#endif //DISTANCECALCULATOR_CUH