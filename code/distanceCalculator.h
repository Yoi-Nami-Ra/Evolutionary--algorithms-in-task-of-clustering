/**
 8/11/2011
 Jaroslaw Wojtasik

 noCuda

 distanceCalculator.h
 **/

#ifndef DISATNCECALCULATOR_H
#define DISTANCECALCULATOR_H

#include "dataLoader.h"

//==============================================
//== Types

/**
 * This struct is to be used for shipping parameters
 * to distance calculating loop.
 */
typedef struct {
	DataStore * dataStore;
	unsigned int gridSize;
	unsigned int blockSize;
} DistancesParams;

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

#endif // DISTANCECALCULATOR_H