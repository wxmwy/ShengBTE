#include <stdio.h>
//#define DEBUG
#ifdef DEBUG
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#else
#define HANDLE_ERROR( err ) (err)
#endif
static void HandleError( cudaError_t err, const char *file, int line )
{

    if (err != cudaSuccess)
    {
    	fprintf(stderr, "ERROR: %s in %s at line %d (error-code %d)\n",
		cudaGetErrorString( err ), file, line, err );
	fflush(stdout);
	exit(-1);
    }
}
