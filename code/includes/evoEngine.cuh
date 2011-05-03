// Header file for Evolutionary Algorith Engine

#define kFieldSize 8 // in bytes
#define kNearestMax 20 // maksimum number of nearest data poinst to be taken in consideration
#define kMaxClusters 10 // maksimum nuber of clusters




//||w|w|w|w|w||nnn|nnn|nnn|nnn|nnn||h||
// Chromoseme
struct Chromosome {
	float klustersweights[ KMaxClusters ];
	FieldContent * clusters; // KMaxClusters * dimensions
	char nearestAttr; // the amount of neares datapoints to take in account
	int ratting; // ratting score for this chromosome
}

// Possible errors
enum EvoError {
	errOk = 0, // no errors
};

// Types of data, recognized by engine
enum AttrDataType {
	eInteger,		// Decimal
	eReal,			// Real
	eIdSortable,	// IDs but still can calculate distance
	eIdComparable,	// IDs of limited values but only equality can be checked - Decimal
	eHash,			// Hash generated from bigger value
};

// Describes structure of a single data point in data storage to be clustered
struct DataPointDescriptor {
	int dimensions;
	AttrDataType * fields;
};

union FieldContent {
	double real;
	long long int integer;
	char * shortText;
	long long int longText;
};

struct NearestPair {
	DataPoint * data;
	int distance;
}

struct DataPoint {
	FieldContent * fields;
	DataPoint * next;
	DataPoint * prev;
	NearestPair nearest[ kNearestMax ];
};


//=======================
// Functions

EvoError GetLastError();

bool InitEvoEngine( DataPointDescriptor * descriptor );

long long int ShortFromString( char * str );

bool AddDataObject ( DataPoint * data );

// int Distance( Datapoint * data1, DataPoint * data2 );

bool PrepareData();

// bool FillNearests();

bool RunEvo( int population, int steps );

bool CleanUp();

