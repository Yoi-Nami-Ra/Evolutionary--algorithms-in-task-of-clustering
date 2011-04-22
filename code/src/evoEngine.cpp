

#include "evoEngine.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>


struct Globals {
	DataPointDescriptor dataDescryptor;
	EvoError lastError;
	DataPoint * root;
	long int dataCount;
	DataPoint ** index;
	int populationSize;
	Chromosome ** currPopulation;

} globals;


//------------------------------------------------------------------------------
bool InitEvoEngine( DataPointDescriptor * descriptor ) {
	globals.dataDescryptor.dimensions = descriptor->dimensions;
	globals.dataDescryptor.fields = new AttrDataType[ descriptor->dimensions ];
	
	for ( int i; i < descriptor->dimensions; i++ ) {
		globals.dataDescryptor.fields[i] = descriptor->fields[i];
	}
	
	globals.root = NULL;
	globals.dataCount = 0;
}

//------------------------------------------------------------------------------
long long int ShortFromString( char * str ) {
	// TODO: crc
}

//------------------------------------------------------------------------------
bool AddDataObject ( DataPoint * data ) {	
	data->next = globals.root;
	data->next->prev = data;
	globals.root = data;
	data->prev = NULL;
	
	globals.dataCount++;
}

//------------------------------------------------------------------------------
int Distance( Datapoint * data1, DataPoint * data2 ) {
	int distance = 0;
	
	for (int i = 0; i < globals.dataDescryptor.dimensions; i++) {
		if ( globals.dataDescryptor.fields[i] == eShortText ) {
			if ( !strcmp( data1->fields[i].shortText, data2->fields[i].shortText ) ) {
				distance++;
			}
		} else {
			if ( data1->fields[i].integer == data2->fields[i].integer ) {
				distance++
			}
		}
	}
	return globals.dataDescryptor.dimensions - distance;
}

//------------------------------------------------------------------------------
bool PrepareData() {
	globals.index = new DataPoint*[ globals.dataCount ];

	DataPoint * curr = globals.root;
	int i = 0;
	while ( curr && i < globals.dataCount ) {
		index[ i++ ] = curr;
		curr = curr->next;
	}

	FillNearests();
}

//------------------------------------------------------------------------------
bool FillNearests() {
	for ( int i = 0; i < globals.dataCount; i++ ) {
		for ( int j = 0; j < globals.dataCount; j++ ) {
			if ( i != j ) {
				int distance = Distance( globals.index[i], globals.index[j] );
				for ( int k = 0; k <kNearestMax; k++ ) {
					if ( globals.index[ i ]->nearest[ k ].data == NULL ) {
						// If there is empty place, just put what we have
						globals.index[ i ]->nearest[ k ].data = globals.index[ j ];
						globals.index[ i ]->nearest[ k ].distance = distance;
						break;
					}

					if ( globals.index[ i ]->nearest[ k ].data == globals.index[ j ] ) {
						// so we won't put things twice in the same list
						break;
					}

					if ( distance < globals.index[ i ]->nearest[ k ].distance ) {
						for ( int l = kNearestMax - 1; l > k; l-- ) {
							globals.index[ i ].nearest[ l ] = globals.index[ i ].nearest[ l - 1 ];
						}
						globals.index[ i ]->nearest[ k ].data = globals.index[ j ];
						globals.index[ i ]->nearest[ k ].distance = distance;
					}
				} // for k
			}
		} // for j
	} // for i
}

//------------------------------------------------------------------------------
bool StartupGeneration() {	
	srand ( time(NULL) );
	int maxStep = globals.dataCount / kMaxClusters;
	int index = 0;
	Chromosome curr = NULL;
	
	globals.currPopulation = new Chromosome*[globals.populationSize];
	for ( int i = 0; i < globals.populationSize; i++ ) {
		globals.currPopulation[i] = curr = new Chromosome;
		index = 0;
		for ( int j = 0; j < kMaxClusters; j++ ) {
			curr.klustersweights[ j ] = (float)( rand() % 100 ) / 100.0;
			index += rand() % maxStep;

			if ( index > globals.dataCount ) {
				index -= globals.dataCount;
			}

			for ( int k = 0; k < globals.dataDescryptor.dimensions; k++ ) {
				curr.clusters[ j * globals.dataDescryptor.dimensions + k ].integer = globals.index[ index ]->fields[ k ].integer;
			}
		} // j
		curr.nearestAttr = rand() % kNearestMax;
	} // i

	// take twice as much, cross them, and there you have
}

//------------------------------------------------------------------------------
bool RateGeneration() {

}

//------------------------------------------------------------------------------
bool RateChromosme() {
}

//------------------------------------------------------------------------------
bool RateCluster() {
	
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
bool RunEvo( int population, int steps ) {

	globals.populationSize = population;
	// TODO: run the algo
	StartupGeneration();

	for ( int i = 0; i < steps; i++ ) {
		RateGeneration();

		Eliminate();

		NewGeneration();
	}
}

//------------------------------------------------------------------------------
bool CleanUp();