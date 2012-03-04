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

#define prepareTests(X) uint testCounterP = 0;			\
	uint testCounter = 0;								\
	bool lastRes = true;								\
	char * testName = X;								\
	logMessage( " === starting Tests: %s", testName );
	

#define makeTest(X, Y) do {																\
		testCounter++;																	\
		if ( !lastRes ) break;															\
		lastRes =  Y;																	\
		logMessage( " test: %u %s %s", testCounter, X, lastRes?"PASSED":"FAILED" );		\
		if ( lastRes ) {																\
			testCounterP++;																\
		}																				\
	} while(0);

#define endTests do {																					\
		logMessage( " passed %u out of %d tests", testCounterP, testCounter );				\
		logMessage( " === Tests: %s %s in general",testName, lastRes?"PASSED":"FAILED" );	\
	} while(0);

//==============================================
//== Types

//==============================================
//== Functions


#endif //DISTANCECALCULATOR_CUH