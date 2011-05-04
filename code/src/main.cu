// Includes
#include <stdio.h>
#include <cutil_inline.h>
#include <shrQATest.h>

#include "dataLoader.cuh"

// Host code
int main(int argc, char** argv)
{
    shrQAStart(argc, argv);

	ErrorCode err;
	
	err = LoadData();

	return 0;
}