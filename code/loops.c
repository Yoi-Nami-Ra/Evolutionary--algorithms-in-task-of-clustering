/**
 11/11/2011
 Jaroslaw Wojtasik

 noCuda

 loops.c
 **/

#include "loops.h"

//==============================================
//== Defines

#define LOOP(A,B,C) for ( context.A.C = 0; context.A.C < loop.B.C; context.A.C++ )
#define LOOP_BLOCK(A) LOOP(blockIdx, gridSize, A)
#define LOOP_THREAD(A) LOOP(threadIdx, blockSize, A)

//==============================================
//== Functions

ErrorCode RunLoop( LoopDefinition loop ) {
	LoopContext context;
	context.params = loop.params;
	context.blockIdx.x = context.blockIdx.y = context.blockIdx.z = 0;
	context.threadIdx.x = context.threadIdx.y = context.threadIdx.z = 0;
	context.err = errOk;
	
	LOOP_BLOCK(z) { 
		LOOP_BLOCK(y) { 
			LOOP_BLOCK(x) {
		LOOP_THREAD(z) { 
			LOOP_THREAD(y) { 
				LOOP_THREAD(x) {
			context.err = errOk; // reset error
			loop.kernel( context );
			if ( context.err != errOk ) {
				logError( context.err, "block(%u,%u,%u) thread(%u,%u,%u) returned with error.",
					context.blockIdx.x, context.blockIdx.y, context.blockIdx.z,
					context.threadIdx.x, context.threadIdx.y, context.threadIdx.z );
			}
		} } }
	} } }

	return context.err;
}
//---------------------------------------------