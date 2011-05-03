

// Includes

#include "errors.cuh"

// Constants

static ErrorCode gLastErrorCode; //< Not a constant but global static

// Functions

ErrorCode SetError( ErrorCode error ) {
	return gLastErrorCode = error;
}
//==================================================================

ErrorCode GetLastErrorCode() {
	return gLastErrorCode;
}
//==================================================================
//== Log Functions

void logError(char* msg) {
	// TODO:
}
//==================================================================

void LogWarning(char* msg) {
	// TODO:
}
//==================================================================

void logDebug(char* msg) {
	//TODO:
}
//==================================================================