#ifndef __CUSP_CG_H
#define __CUSP_CG_H

#ifdef __CUDACC__
#include <cusp/hyb_matrix.h>
#include <cusp/gallery/poisson.h>
#include <cusp/krylov/cg.h>
#endif

// where to perform the computation
 typedef cusp::device_memory MemorySpace;
// which floating point type to use
 typedef float ValueType;

 typedef cusp::coo_matrix<int, ValueType, MemorySpace> CUSP_matrixtype;
 typedef cusp::coo_matrix<int, ValueType, cusp::host_memory> CUSP_hostmatrixtype;
//typedef cusp::hyb_matrix<int, ValueType, MemorySpace> CUSP_matrixtype;
typedef cusp::array1d<ValueType, MemorySpace> CUSP_arraytype;
typedef cusp::array1d<ValueType, cusp::host_memory> CUSP_hostarraytype;

void _CUSP_CG( CUSP_matrixtype &A, CUSP_arraytype &x, CUSP_arraytype &b );

#include "cuspcg.h"


#define eps  0.00001
#define ABS(a)  ( ((a) > 0.0) ? (a) : -(a) )
#define NEQ(a,b)  (  (abs(b-a) < eps) ? true : false )
#endif
