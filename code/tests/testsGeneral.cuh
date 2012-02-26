/**
 25/02/2012
 Jaroslaw Wojtasik

 Cuda_test

 testGeneral.cuh
 **/

#ifndef DISTANCECALCULATOR_TEST_CUH
#define DISTANCECALCULATOR_TEST_CUH

#include "errors.cuh"

//==============================================
//== Globals

#define prepareTests uint testCounterP = 0; uint testCounter = 0; bool lastRes = true;

#define makeTest(x, y) do {															\
		testCounter++;																\
		if ( !lastRes ) break;														\
		bool lastRes = y;															\
		logMessage( " test: %u %s %s", testCounter, x, lastRes?"PASSED":"FAILED" );	\
		if ( lastRes ) {															\
			testCounterP++;															\
		} else {																	\
			goto exit;																\
		}																			\
	} while(0);

#define endTests exit:																\
	do {																			\
		logMessage( " passed %u out of %d tests", testCounterP, testCounter );		\
		logMessage( " %s in general", lastRes?"PASSED":"FAILED" );					\
	} while(0);

//==============================================
//== Types

//==============================================
//== Functions


#endif //DISTANCECALCULATOR_CUH