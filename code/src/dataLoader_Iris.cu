// Module of data loader specialized in loading Iris database

#include <stdlib.h>
#include "errors.cuh"
#include "dataLoader_Iris.cuh"

// Constants

static const char* kdataFilePath = "./data/iris.data";
static const char* kConvFilePath = "./data/iris_converted.data";
static dataStore* gCurrDataStore = 0;

// Functions
// Local Declarations

/*
 * Load data in already converted state.
 */
ErrorCode LoadConverted();

/*
 * Convert the data.
 */
ErrorCode ConvertData();
// Local Definitions

bool CheckAlreadyConverted() {
	FILE* file = fopen( kConvFilePath, "r" );
	if ( file != 0 ) {
		fclose( file );
		return true;
	}
	return false;
}
//====================================================================

ErrorCode startLoadingData() {
	 ErrorCode err;
	 if ( CheckAlreadyConverted() ) {
		 err = LoadConverted();
	 } else {
		 err = ConvertData();
	 }

	 return err;
}
//====================================================================

 dataStore* GetCurrDataStore() {
	 return gCurrDataStore;
 }
 //====================================================================
 
ErrorCode LoadConverted() {	
	if ( gCurrDataStore == 0 ) {
		// prepare data store
		gCurrDataStore = (dataStore*)malloc( sizeof( dataStore ));
		if ( gCurrDataStore == 0 ) {
			logError("No memory for data store descriptor");
			return SetError( errNoMemory );
		}

		gCurrDataStore->info.numEntries = 150; // there are 150 entries in Iris database
		gCurrDataStore->dataVector = (dataEntry*)malloc( sizeof(dataEntry) * gCurrDataStore->info.numEntries );

		if ( gCurrDataStore->dataVector == 0 ) {
			logError("No memory for data vector");
			free( gCurrDataStore );
			return SetError( errNoMemory );
		}
	}

	FILE* dataFile = fopen( kConvFilePath, "r" );

	if ( dataFile == 0 ) {
		free( gCurrDataStore );
		free( gCurrDataStore->dataVector );
		logError("File with data Not found");
		return SetError( errFileNotFound );
	}

	size_t res;
	// numEntries
	res = fread( &gCurrDataStore->info.numEntries, sizeof(unsigned int), 1, dataFile );
	if ( res <= 0 ) {
		free( gCurrDataStore );
		free( gCurrDataStore->dataVector );
		logError("Failed to read converted data");
		fclose( dataFile );
		return SetError( errFileCorupted );	
	}
	// vector
	res = fread( gCurrDataStore->dataVector, sizeof(dataEntry), gCurrDataStore->info.numEntries, dataFile );
	if ( res <= 0 ) {
		free( gCurrDataStore );
		free( gCurrDataStore->dataVector );
		logError("Failed to read converted data");
		fclose( dataFile );
		return SetError( errFileCorupted );	
	}

	fclose( dataFile );


	return errOk;
}
//====================================================================

ErrorCode ConvertData() {
	if ( gCurrDataStore == 0 ) {
		// prepare data store
		gCurrDataStore = (dataStore*)malloc( sizeof( dataStore ));
		if ( gCurrDataStore == 0 ) {
			logError("No memory for data store descriptor");
			return SetError( errNoMemory );
		}

		gCurrDataStore->info.numEntries = 150; // there are 150 entries in Iris database
		gCurrDataStore->dataVector = (dataEntry*)malloc( sizeof(dataEntry) * gCurrDataStore->info.numEntries );

		if ( gCurrDataStore->dataVector == 0 ) {
			free( gCurrDataStore );
			logError("No memory for data vector");
			return SetError( errNoMemory );
		}
	}

	FILE* dataFile = fopen( kdataFilePath, "r" );
	if ( dataFile == 0 ) {
		free( gCurrDataStore->dataVector );
		free( gCurrDataStore );
		logError("File with data Not found");
		return SetError( errFileNotFound );
	}

	// load data into memory
	int index = 0;	
	char str[20];
	while( !feof( dataFile ) && ( index < gCurrDataStore->info.numEntries )) {
		fscanf( dataFile, "%3f,%3f,%3f,%3f,%s",
			&gCurrDataStore->dataVector[index].a,
			&gCurrDataStore->dataVector[index].b,
			&gCurrDataStore->dataVector[index].c,
			&gCurrDataStore->dataVector[index].d, str );
		index++;
	}
	fclose( dataFile );

	// save data into binary file
	dataFile = fopen( kConvFilePath, "w" );

	if ( dataFile == 0 ) {		
		LogWarning( "Failed opening converted data file for writting" );
		//return SetError( errFileWrite );		
		return errOk;
	}

	size_t res;
	// numEntries
	res = fwrite( &gCurrDataStore->info.numEntries, sizeof(unsigned int), 1, dataFile );
	if ( res <= 0 ) {
		free( gCurrDataStore->dataVector );
		free( gCurrDataStore );
		logError( "Failed to write down converted data" );
		fclose( dataFile );
		//return SetError( errFileWrite );	
		return errOk;
	}
	// vector
	res = fwrite( gCurrDataStore->dataVector, sizeof(dataEntry), gCurrDataStore->info.numEntries, dataFile );
	if ( res <= 0 ) {
		free( gCurrDataStore->dataVector );
		free( gCurrDataStore );
		logError("Failed to write down converted data");
		fclose( dataFile );
		//return SetError( errFileWrite );	
		return errOk;
	}

	fclose( dataFile );
	return errOk;
}
//====================================================================

ErrorCode releaseDataStore() {
	if ( gCurrDataStore ) {
		free( gCurrDataStore->dataVector );
		free( gCurrDataStore );
	}
	return errOk;
}
//====================================================================