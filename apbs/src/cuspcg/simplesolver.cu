#include <iostream>
#include "cuspcg_cuda.h"
#include <cusp/coo_matrix.h>
#include <cusp/print.h>
#include <cusp/krylov/gs.h>
#include <cusp/multiply.h>


template< typename V >
void 
buildCOOMatrix( CUSP_matrixtype &A, int nx, // A is a nx*nx array of Vs
	V *vA )
{
	int entry = 0;
	CUSP_hostmatrixtype hmtx; // host matrix
	for(int i = 0; i < nx; i++ )
	{
		for(int j = 0; j < nx; j++)
		{
			int idx = i * nx + j;
			if( vA[idx] != 0.0 )
			{
				hmtx.row_indices[entry] = i;
				hmtx.column_indices[entry] = j;
				hmtx.values[entry] = vA[idx];
				entry++;
			}
		}
	}
	A = hmtx;
}

template< typename V >
void 
buildCOOMatrix( CUSP_matrixtype &A, int numNonzeroes, // A is a nx*nx array of Vs
	V *vals, int *rows, int *cols )
{
	CUSP_hostmatrixtype hmtx; // host matrix
	for(int i = 0; i < numNonzeroes; i++ )
	{
		hmtx.row_indices[i] = rows[i];
		hmtx.column_indices[i] = cols[i];
		hmtx.values[i] = vals[i];
	}
	A = hmtx;
}

void* 
newCuspCOOMatrix(int rank, int numNonzeroes)
{
	CUSP_matrixtype *ret = new CUSP_matrixtype(rank,rank,numNonzeroes);

	return (void*)ret;
}

void 
deleteCuspCOOMatrix( void *m )
{
	CUSP_matrixtype *hm = static_cast<CUSP_matrixtype*>(m);
	delete hm;
}

void 
deleteCuspCOOHostMatrix( void *m )
{
	CUSP_hostmatrixtype *hm = static_cast<CUSP_hostmatrixtype*>(m);
	delete hm;
}

void* newCuspCOOHostMatrix(int rank, int numNonzeroes)
{
	CUSP_hostmatrixtype *ret = new CUSP_hostmatrixtype(rank,rank,numNonzeroes);

	return (void*)ret;
}

double 
getCuspCOOMatrixVal( void* matrixptr, int idx)
{
	CUSP_matrixtype* matrix = (CUSP_matrixtype*)matrixptr;
	int row = matrix->row_indices[idx] , col = matrix->column_indices[idx];
	double val = matrix->values[idx];
	std::cout << "(" << row << ", " << col << ") = " << val << std::endl;
	return val;
}

#if 0
void
setCuspCOOMatrixVal( void* matrixptr, int idx, int row, int col, double val )
{
	CUSP_matrixtype* matrix = (CUSP_matrixtype*)matrixptr;

	std::cout << "Setting matrix[" << row << "][" << col << "] = " << val << std::endl;
	
	matrix->row_indices[idx] = row;
	matrix->column_indices[idx] = col;
	matrix->values[idx] = (ValueType)val;
//printf("Set (%d,%d) <- %.3f.  idx is %d\n", row,col,val,idx);
}
#endif


void
setCuspCOOHostMatrixVal( void* matrixptr, int idx, int row, int col, double val )
{
	CUSP_hostmatrixtype* matrix = (CUSP_hostmatrixtype*)matrixptr;

	
	matrix->row_indices[idx] = row;
	matrix->column_indices[idx] = col;
	matrix->values[idx] = (ValueType)val;
//printf("Set (%d,%d) <- %.3f.  idx is %d\n", row,col,val,idx);
}


void*
newCUSParray(int numElems)
{
	return new CUSP_arraytype(numElems);
}

void
deleteCUSParray( void* a )
{
	CUSP_arraytype* ca = static_cast<CUSP_arraytype*>(a);
	delete ca;
}

void*
newCUSPhostarray(int numElems)
{
	return new CUSP_hostarraytype(numElems);
}

void deleteCUSPhostarray( void* a )
{
	CUSP_hostarraytype* ca = static_cast<CUSP_hostarraytype*>(a);
	delete ca;
}


void 
copyHostToDeviceMatrix( void* hmtx, void* dmtx )
{
	CUSP_hostmatrixtype *H = (CUSP_hostmatrixtype*)hmtx;
	CUSP_matrixtype *D = (CUSP_matrixtype*)dmtx;
	*D = *H;
}

void
setCUSParrayVal( void* arrayptr, int idx, double val )
{
	CUSP_arraytype *array = ((CUSP_arraytype*)arrayptr);
	array->operator[](idx) = val;
	
}

void
setCUSPhostarrayVal( void* arrayptr, int idx, double val )
{
	CUSP_hostarraytype *array = ((CUSP_hostarraytype*)arrayptr);
	array->operator[](idx) = val;
	
}

void 
copyArrayh2d( void* host, void* dev )
{
	CUSP_arraytype  *Dev = (CUSP_arraytype*) dev ; 
	CUSP_hostarraytype *Host = (CUSP_hostarraytype*) host;
	
	*Dev = *Host;
}

double
getCUSParrayVal( void* arrayptr, int idx)
{
	//CUSP_arraytype array = *((CUSP_arraytype*)arrayptr);
	//return array[idx];
	CUSP_arraytype *array = ((CUSP_arraytype*)arrayptr);
	return array->operator[](idx);
}

template< typename V >
void
buildCUSPArray( CUSP_arraytype &b, int nx, V *vB )
{
	for(int i = 0; i < nx; i++)
	{
		b[i] = vB[i];
	}
}

void
unpackCUSPArray( CUSP_arraytype *x, int nx, double *vx )
{
	CUSP_hostarraytype hx = *x;
	for(int i = 0; i < nx; i++)
	{
//		vx[i] = x->operator[](i);
		vx[i] = hx[i];
	}
}

/**
* Call the CUSP GS solver on a sparse matrix in COO form.  
* 
* @param numvals  The number of elements in the parallel arrays vals, cols,and rows.
* Also equals the number of nonzero elements in the sparse matrix.
*
* @param vals     Vals, cols, and rows are a set of parallel arrays.  They must hold the
* same number of elements, otherwise behavior is undefined.   Vals holds the nonzero
* elements of the sparse matrix.
*
* @param cols     Vals, cols, and rows are a set of parallel arrays.  They must hold the
* same number of elements, otherwise behavior is undefined.  Cols holds the column number
* of the corresponding element
*
* @param rows     Vals, cols, and rows are a set of parallel arrays.  They must hold the
* same number of elements, otherwise behavior is undefined.   Rows holds the row number of the
* corresponding element.
*/
void CUSP_GS( int* arows,
	      int* acols,
	      double* avals,
	      int numNonzeroes,
	      int xrank,
	      double *_x,
	      double *_b )
{
    // allocate device memory for CSR format
    int* device_I; cudaMalloc(&device_I, numNonzeroes * sizeof(int));
    int* device_J; cudaMalloc(&device_J, numNonzeroes * sizeof(int));
    double* device_V; cudaMalloc(&device_V, numNonzeroes * sizeof(double));
    // allocate device memory for x and y arrays
    double * device_x; cudaMalloc(&device_x, xrank * sizeof(double));
    double * device_b; cudaMalloc(&device_b, xrank * sizeof(double));

  // copy raw data from host to device
  cudaMemcpy(device_I, arows, numNonzeroes * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(device_J, acols, numNonzeroes * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(device_V, avals, numNonzeroes * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_x, _x, xrank * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_b, _b, xrank * sizeof(double), cudaMemcpyHostToDevice);

  
  // *NOTE* raw pointers must be wrapped with thrust::device_ptr!
  thrust::device_ptr<int> wrapped_device_I(device_I);
  thrust::device_ptr<int> wrapped_device_J(device_J);
  thrust::device_ptr<double> wrapped_device_V(device_V);
  thrust::device_ptr<double> wrapped_device_x(device_x);
  thrust::device_ptr<double> wrapped_device_b(device_b);
 
   // use array1d_view to wrap the individual arrays
  typedef typename cusp::array1d_view< thrust::device_ptr<int> > DeviceIndexArrayView;
  typedef typename cusp::array1d_view< thrust::device_ptr<double> > DeviceValueArrayView;

  DeviceIndexArrayView row_indices (wrapped_device_I, wrapped_device_I + numNonzeroes);
  DeviceIndexArrayView column_indices(wrapped_device_J, wrapped_device_J + numNonzeroes);
  DeviceValueArrayView values (wrapped_device_V, wrapped_device_V + numNonzeroes);
  DeviceValueArrayView x (wrapped_device_x, wrapped_device_x + xrank);
  DeviceValueArrayView b (wrapped_device_b, wrapped_device_b + xrank);


  // combine the three array1d_views into a coo_matrix_view
  typedef cusp::coo_matrix_view<DeviceIndexArrayView,
                                DeviceIndexArrayView,
                                DeviceValueArrayView> DeviceView;

  // construct a coo_matrix_view from the array1d_views
  DeviceView A(xrank, xrank, numNonzeroes, row_indices, column_indices, values);
    
  // set stopping criteria: iteration_limit = 100, relative_tolerance = 1e-5
  cusp::verbose_monitor<double> monitor(b, 100, 1e-3);


  // solve the linear system A * x = b with the Gauss-Seidel method
  cusp::krylov::gs(A, wrapped_device_x, wrapped_device_b, xrank, monitor);

  // copy the solution back to the host
  cudaMemcpy(_x, device_x, xrank * sizeof(double), cudaMemcpyDeviceToHost);

  // free device arrays
  cudaFree(device_I);
  cudaFree(device_J);
  cudaFree(device_V);
  cudaFree(device_x);
  cudaFree(device_b);

}




/**
* Call the CUSP CG solver on a sparse matrix in COO form.  
* 
* @param numvals  The number of elements in the parallel arrays vals, cols,and rows.
* Also equals the number of nonzero elements in the sparse matrix.
*
* @param vals     Vals, cols, and rows are a set of parallel arrays.  They must hold the
* same number of elements, otherwise behavior is undefined.   Vals holds the nonzero
* elements of the sparse matrix.
*
* @param cols     Vals, cols, and rows are a set of parallel arrays.  They must hold the
* same number of elements, otherwise behavior is undefined.  Cols holds the column number
* of the corresponding element
*
* @param rows     Vals, cols, and rows are a set of parallel arrays.  They must hold the
* same number of elements, otherwise behavior is undefined.   Rows holds the row number of the
* corresponding element.
*/
void CUSP_CG( void* p_cuspA,
	      void* p_cuspx,
	      void* p_cuspb,
	      int nx,
	      double *vals_x )
{

/*	CUSP_matrixtype cuspA( nx, nx, numvals_A );
	CUSP_arraytype cuspx( nx, 0 );
	CUSP_arraytype cuspb( nx, 1 );
	buildCOOMatrix( cuspA, numvals_A, vals_A, rows_A, cols_A );
	buildCUSPArray( cuspx, nx, vals_x );
	buildCUSPArray( cuspb, nx, vals_b );
*/
	CUSP_matrixtype *cuspA = (CUSP_matrixtype*)p_cuspA;
	CUSP_arraytype *cuspx = (CUSP_arraytype*)p_cuspx;
	CUSP_arraytype *cuspb = (CUSP_arraytype*)p_cuspb;

	_CUSP_CG( *cuspA, *cuspx, *cuspb );
	unpackCUSPArray( cuspx, nx, vals_x );
}
	

#if 0
void CUSP_CG( float* A, 
		float *x,
		int nx,  //  A is a nx * nx array of ints
		float* b	//  b is a 'nx' element array of ValueTypes
				// We solve for Ax=b
		)
{
	int numNonzeroes = 0; // count # of nonzeroes in A
	int nxnx = nx * nx;
	for(int i = 0; i < nxnx; i++ )
	{
		if( A[i] != 0.0 ) numNonzeroes++;
	}
	CUSP_matrixtype cuspA( nx, nx, numNonzeroes );
	CUSP_arraytype cuspx( nx, 0 );
	CUSP_arraytype cuspb( nx, 1 );
	buildCOOMatrix( cuspA, nx, A );
	buildCUSPArray( cuspx, nx, x );
	buildCUSPArray( cuspb, nx, b );

//	cusp::multiply( cuspA, cuspx, cuspb );
	

	_CUSP_CG( cuspA, cuspx, cuspb );
	unpackCUSPArray( cuspx, nx, x );

}
#endif

void 
setCuspCOOHostMatrixValues7x( void* vmtx, // CUSP_hostmatrixtype cast to void*
			      int* rows,
			      int* cols,
			      double* vals, // 
			      int numelems ) // number of elements in vals, rows, and cols (these are parallel arrays)
{
	CUSP_hostmatrixtype* matrix = (CUSP_hostmatrixtype*)vmtx;
	for(int idx = 0; idx < numelems; idx++)
	{
		matrix->row_indices[idx] = rows[idx];
		matrix->column_indices[idx] = cols[idx];
		matrix->values[idx] = vals[idx];
	}
}
			

void _CUSP_CG( CUSP_matrixtype &A, CUSP_arraytype &x, CUSP_arraytype &b )
{
    // set stopping criteria:
    //  iteration_limit    = 100
    //  relative_tolerance = 1e-3
//    cusp::verbose_monitor<ValueType> monitor(b, 1000, 1e-3);
    cusp::default_monitor<ValueType> monitor(b, 1000, 1e-4);
    cusp::identity_operator<ValueType, MemorySpace> M(A.num_rows, A.num_rows);
#if DEBUG
	std::cout << "Trying to solve: " << std::endl;
	cusp::print( A );
std::cout << "this is x: " << std::endl;
cusp::print(x);
#endif

    // solve the linear system A * x = b with the Conjugate Gradient method
    cusp::krylov::cg(A, x, b, monitor, M);
//cusp::multiply(A,x,b);
//    cusp::krylov::cg(A, x, b);
	
}


int solvertest(void)
{
    // create an empty sparse matrix structure (HYB format)
    //cusp::hyb_matrix<int, ValueType, MemorySpace> A;
   CUSP_matrixtype A;

    // create a 2d Poisson problem on a 10x10 mesh
    cusp::gallery::poisson5pt(A, 10, 10);

    // allocate storage for solution (x) and right hand side (b)
    cusp::array1d<ValueType, MemorySpace> x(A.num_rows, 0);
    cusp::array1d<ValueType, MemorySpace> b(A.num_rows, 1);

/*
    // set stopping criteria:
    //  iteration_limit    = 100
    //  relative_tolerance = 1e-3
    cusp::verbose_monitor<ValueType> monitor(b, 100, 1e-3);

    // set preconditioner (identity)
    cusp::identity_operator<ValueType, MemorySpace> M(A.num_rows, A.num_rows);

    // solve the linear system A * x = b with the Conjugate Gradient method
    cusp::krylov::cg(A, x, b, monitor, M);

*/
    
    _CUSP_CG( A,x,b );
    return 0;
}

