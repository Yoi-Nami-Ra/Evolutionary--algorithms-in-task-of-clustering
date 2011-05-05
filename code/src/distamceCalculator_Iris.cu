


// Module responsible for loading data

//==============================================
//== Includes
#include "globals.cuh"
#include "errors.cuh"
#include "distanceCalculator_Iris.cuh"
#include "dataLoader.cuh"

//==============================================
//== Types

//==============================================
//== Constants
float* sDistancesVector;

//==============================================
//== Functions
ErrorCode StartCalculatingDistances() {

	ErrorCode err = LoadData();
	if (err != errOk) {
		return err;
	}

	data = GetCurrDataStore();
	if (data == 0) {
		return 
	}

	return errOk;
}
//==============================================