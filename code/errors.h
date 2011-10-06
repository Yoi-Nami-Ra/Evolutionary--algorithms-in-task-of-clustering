/**
 03/10/2011
 Jaroslaw Wojtasik

 noCuda

 errors.h
 **/

#ifndef ERRORS_H
#define ERRORS_H

//==============================================
//== Types

typedef enum ErrorCode {
	errOk = 0,			// Everything under control
	errGeneral,			// When we can't say what's wrong
	errFileNotFound,	// Requestet file not found
	errFileCorupted,	// Formating of requested file is not as it was expected
	errFileWrite,		// We failed to write into file
	errFileRead,		// We failed to read from file
	errNoMemory,		// There is no enough memory
	errDataNotReady,	// LoadData didn't fail but still we receive nil
	errNoData,			// In place where data should be procesed we got empty
} ErrorCode;

//==============================================
//== Functions

/*
 * Sets Error code to be retrived as the last one.
 */
ErrorCode SetLastError( ErrorCode error );

/*
 * Returns the last set error code.
 */
ErrorCode GetLastErrorCode();

//== Log Functions

#define timeStamp

#define logMessage(A, ...) printf(timeStamp " " A "\n", __VA_ARGS__)

/*
 * Log given message as an error.
 */
#define logError(A, B, ...) logMessage("[E %X] " B, A, __VA_ARGS__)

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


#endif //ERRORS_H