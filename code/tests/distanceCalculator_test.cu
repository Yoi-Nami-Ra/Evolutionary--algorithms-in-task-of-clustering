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
bool test01 ();

//============================================================================
//== Functions

void runDistancesTests() {
	prepareTests("Distances");

	makeTest( "Texture binding", test01() );

	endTests;
}
//----------------------------------------------------------------------------

bool test01 () {

	return testRawTexturesKernel();
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
