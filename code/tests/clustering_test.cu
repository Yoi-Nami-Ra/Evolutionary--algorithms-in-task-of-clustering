/**
 12/05/2012
 Jaroslaw Wojtasik

 Cuda_test

 clustering_test.cu
 **/

//============================================================================
//== Includes
#include <cutil_inline.h>
#include <cutil_math.h>
#include <cuda.h>
#include "testsGeneral.cuh"
#include "clustering.cuh"

//============================================================================
//== Globals
//============================================================================
//== Declarations

bool cTest01();
bool cTest02();

//============================================================================
//== Functions

void runClusteringTests() {
	prepareTests("Clustering");

	makeTest( "Test enviroment preparation", cTest02() );

	makeTest( "Membership and density", cTest01() );

	endTests;
}
//----------------------------------------------------------------------------

bool cTest01() {
	return testMembershipAndDensity();
}
//----------------------------------------------------------------------------

bool cTest02() {
	return testEnviroment();
}
//----------------------------------------------------------------------------


//============================================================================

//----------------------------------------------------------------------------
