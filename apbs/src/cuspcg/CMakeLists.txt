cmake_minimum_required (VERSION 2.8.8)


SET ( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} )
SET ( CUSP_HOME "/usr/local/cusp/cusplibrary" )

find_package(CUDA)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-I${CUSP_HOME};-O0 -g;-gencode arch=compute_20,code=sm_20)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-I${CUSP_HOME};-O3;-Xcompiler -fopenmp;-gencode arch=compute_20,code=sm_20)
cuda_compile(SS_O simplesolver.cu)

#SET ( CMAKE_EXE_LINKER_FLAGS  "-L/usr/local/cuda/lib64 -lcudart") 
#add_executable(ss ssdriver.C  ${SS_O} )
#add_library( cuspcg ssdriver.C ${SS_O} )
add_library( cuspcg cuspcgtest.C ${SS_O})
