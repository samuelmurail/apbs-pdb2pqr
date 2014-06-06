#pragma once 

#ifdef __cplusplus
extern "C"{
#endif
        void CUSP_CG( void* p_cuspA,
              void* p_cuspx,
              void* p_cuspb,
              int nx,
              double *vals_x );
        void setCUSParrayVal( void* array, int idx, double val );
        void setCUSPhostarrayVal( void* array, int idx, double val );
	void copyArrayh2d( void* host, void* dev );
        double getCUSParrayVal( void* array, int idx);
        void* newCUSParray(int size);
        void deleteCUSParray( void* a );
        void* newCUSPhostarray(int size);
        void deleteCUSPhostarray( void* a );
        void* newCuspCOOMatrix(int rank, int numNonzeroes);
        void deleteCuspCOOMatrix( void *m );
        void* newCuspCOOHostMatrix(int rank, int numNonzeroes);
        void deleteCuspCOOHostMatrix( void *m );
        //void setCuspCOOMatrixVal( void* matrixptr, int idx, int row, int col, double val);
	void setCuspCOOHostMatrixVal( void* matrixptr, int idx, int row, int col, double val );
	void setCuspCOOHostMatrixValues7x( void* vmtx, int* rows, int* cols, double* vals, int numelems ) ;
        void CUSP_GS( int* arows, int* acols, double* avals, int numNonzeroes, int xrank, double *_x, double *_b );
	void copyHostToDeviceMatrix( void* hmtx, void* dmtx );


#ifdef __cplusplus
}
#endif
