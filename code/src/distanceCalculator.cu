// Module responsible for loading data

//==============================================
//== Includes
#include "globals.cuh"
#include "errors.cuh"
#include "distanceCalculator.cuh"
#ifdef DATABASE_IRIS
#include "distanceCalculator_Iris.cuh"
#endif // DATABASE_IRIS

//==============================================
//== Types

//==============================================
//== Constants

//==============================================
//== Functions

ErrorCode calculateDistances() {
	return startCalculatingDistances();
}

