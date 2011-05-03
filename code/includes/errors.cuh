// Error codes
#ifndef ERRORS_CUH
#define ERRORS_CUH

// Types



enum ErrorCode {
	errOk = 0,			// Everything under control
	errGeneral,			// When we can't say what's wrong
	errFileNotFound,	// Requestet file not found
	errFileCorupted,	// Formating of requested file is not as it was expected
	errFileWrite,		// We failed to write into file
	errFileRead,		// We failed to read from file
	errNoMemory,		// There is no enough memory
};

//==============================================
//== Functions

/*
 * Sets Error code to be retrived as the last one.
 */
ErrorCode SetError( ErrorCode error );

/*
 * Returns the last set error code.
 */
ErrorCode GetLastErrorCode();

//== Log Functions

/*
 * Log given message as an error.
 */
void logError(char* msg);

/*
 * Log given message as a warning.
 */
void LogWarning(char* msg);

/*
 * Log given message as a debug entry
 */
void logDebug(char* msg);

#endif // ERRORS_CUH