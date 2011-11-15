/**
 03/10/2011
 Jaroslaw Wojtasik

 noCuda

 errors.h
 **/

#ifndef ERRORS_H
#define ERRORS_H

#include <stdio.h>

//==============================================
//== Types

typedef enum {
	errOk = 0,			// Everything under control
	errGeneral,			// When we can't say what's wrong
	errFileNotFound,	// Requestet file not found
	errFileCorupted,	// Formating of requested file is not as it was expected
	errFileWrite,		// We failed to write into file
	errFileRead,		// We failed to read from file
	errNoMemory,		// There is no enough memory
	errDataNotReady,	// LoadData didn't fail but still we receive nil
	errBadParam,		// Parameters we received are wrong
	errNoData,			// In place where data should be procesed we got empty
	errOutOfRange,		// We're asking for element pointing outside the range
	errWrongParameter,	// The parameter we got is wrong.
	errFailProcessData,	// Failed to process data.
} ErrorCode;

//==============================================
//== Functions

/*
 * Sets Error code to be retrived as the last one.
 */
ErrorCode SetLastErrorCode( ErrorCode error );

/*
 * Returns the last set error code.
 */
ErrorCode GetLastErrorCode();

/*
 * Returns short and simple error description
 */
char * ErrorDesc( ErrorCode errorType );

//---------------------------------------------------------
//== Log Functions

#define timeStamp

#define logMessage(A, ...) printf(timeStamp " " A "\n", __VA_ARGS__)

/*
 * Log given message as an error.
 */
#define logError(A, B, ...) logMessage("[E %s]\n%s:%d - %s\n" B, ErrorDesc(A), __FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)

/*
 * Log given message as a warning.
 */
#define LogWarning(A, ...) logMessage("[W] " A, __VA_ARGS__)

/*
 * Log given message as a debug entry
 */
#ifdef _DEBUG
#define logDebug(A, ...) logMessage("[D] " A, __VA_ARGS__)
#else
 #define logDebug(A, ...)
#endif //_DEBUG

#define checkAlloc(A) if(A==NULL) {\
	SetLastErrorCode( errNoMemory );\
	logError( errNoMemory, "A" );

#define reportError(A, B, ...) SetLastErrorCode(A); logError(A, B, __VA_ARGS__)

#define checkError

#endif //ERRORS_H