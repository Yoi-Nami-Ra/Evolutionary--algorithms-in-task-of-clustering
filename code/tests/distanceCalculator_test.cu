/**
 22/02/2012
 Jaroslaw Wojtasik

 Cuda_test

 distanceCalculator_test.cu
 **/

//============================================================================
//== Includes
#include <cutil_inline.h>
#include <cutil_math.h>
#include <cuda.h>
#include "testsGeneral.cuh"
#include "dataLoader.cuh"
#include "distanceCalculator.cuh"

//============================================================================
//== Globals
//============================================================================
//== Declarations

/*
 * Tests raw data binding to texture
 */
bool dTest01();

bool dTest02();

bool dTest03();

bool dTest04();

//============================================================================
//== Functions

void runDistancesTests() {
	prepareTests("Distances");

	makeTest( "Texture binding", dTest01() );

	makeTest( "Distances calculation", dTest02() );

	makeTest( "Distances binding", dTest03() );

	makeTest( "Neighbours calculation", dTest04() );

	endTests;
}
//----------------------------------------------------------------------------

bool dTest01() {
	return testRawTextures();
}
//----------------------------------------------------------------------------

bool dTest02() {
	return testDistanceCalculation();
}
//----------------------------------------------------------------------------

bool dTest03() {
	return testNeighbourCalculation();
}
//----------------------------------------------------------------------------

bool dTest04() {
	return testDistanceBinding();
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
