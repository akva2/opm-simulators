# defines that must be present in config.h for our headers
set (opm-simulators_CONFIG_VAR
  HAVE_MPI
  COMPILE_GPU_BRIDGE
  HAVE_AVX2_EXTENSION
  HAVE_CUDA
  HAVE_OPENCL
  HAVE_OPENCL_HPP
  HAVE_AMGCL
  HAVE_AMGX
  HAVE_VEXCL
  HAVE_ROCALUTION
  HAVE_ROCSPARSE
  HAVE_SUITESPARSE_UMFPACK
  HAVE_DAMARIS
  HAVE_HDF5
  HAVE_HYPRE
  HAVE_DUNE_ISTL
  HAVE_DUNE_COMMON
  USE_HIP
  FLOW_INSTANTIATE_FLOAT
)

find_package(Boost COMPONENTS date_time REQUIRED)
find_package(dune-common REQUIRED)
find_package(dune-istl REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(SuiteSparse COMPONENTS UMFPACK REQUIRED)
find_package(opm-grid REQUIRED)
find_package(fmt)
find_package(HDF5)
find_package(MPI)
find_package(SuperLU)

if(TARGET opmsimulators)
else()
  if(USE_GPU_BRIDGE)
    find_package(rocalution)
    find_package(rocblas)
    find_package(rocsparse)
  endif()
  if(USE_DAMARIS_LIB)
    find_package(Damaris 1.9)
  endif()
endif()
