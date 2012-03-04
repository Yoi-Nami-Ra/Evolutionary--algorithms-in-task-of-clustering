/**
 03/10/2011
 Jaroslaw Wojtasik

 noCuda

 errors.c
 **/

#include "errors.cuh"
#include <stdio.h>
#include <string.h>

//==============================================
//== Globals

static ErrorCode gCurrError;
static char * gLogFile = 0;

//==============================================
//== Functions

/*
 * Allows to set a name for file where logs will be written.
 */
void SetLogFile( char * filename ) {
	if ( gLogFile ) free( gLogFile );

	gLogFile = strdup( filename );
}

const char * GetLogFile() {
	return gLogFile;
}

/*
 * Sets Error code to be retrived as the last one.
 */
ErrorCode SetLastErrorCode( ErrorCode error ) {
	gCurrError = error;
	return gCurrError;
}

/*
 * Returns the last set error code.
 */
ErrorCode GetLastErrorCode( void ) {
	return gCurrError;
}

char * ErrorDesc( ErrorCode errorType ) {
	char * res = NULL;
	switch ( errorType ) {
		case errOk: {
			res = "No Error";
		} break;
		case errGeneral: {
			res = "General Error";
		} break;
		case errFileNotFound: {
			res = "File Not Found";
		} break;
		case errFileCorupted: {
			res = "File Corupted";
		} break;
		case errFileWrite: {
			res = "File Write Error"; 
		} break;
		case errFileRead: {
			res = "File Read Error";
		} break;
		case errNoMemory: {
			res = "Out Of Memory";
		} break;
		case errDataNotReady: {
			res = "Data Not Ready";
		} break;
		case errNoData: {
			res = "No Data";
		} break;
		default: {
			res = "Unknown Error";
		}
	} // switch

	return res;
}
