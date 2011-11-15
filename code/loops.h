/**
 11/11/2011
 Jaroslaw Wojtasik

 noCuda

 loops.h
 **/

#include "errors.h"

//==============================================
//== Types

typedef struct {
	unsigned int x;
	unsigned int y;
	unsigned int z;
} Point;

/**
 * Used to hold context for each iteration of the loop
 */
typedef struct {
	Point blockIdx; ///< Index of current block
	Point threadIdx; ///< Index of current thread
	void * params; ///< Here custom params to be stored
	ErrorCode err; ///< Error if any for single run of a loop;
} LoopContext;

typedef void (*KernelFunc)( LoopContext );

/**
 * USed to define loops size and dimension.
 */
typedef struct {
	Point blockSize; ///< threads dimensions
	Point gridSize; ///< blocks dimensions
	KernelFunc kernel;
	void * params;
} LoopDefinition;

//==============================================
//== Functions

/**
 * Runs a single loop using given definition;
 * @param loop	- [in] definition of the loop
 *
 * @return error code if any
 */
ErrorCode RunLoop( LoopDefinition loop );