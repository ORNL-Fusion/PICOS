#-------------------------------------------------------------------------------
#  Picos Library
#-------------------------------------------------------------------------------
add_library (picos)
target_compile_features (picos

                         PUBLIC

                         cxx_std_23
)
target_compile_definitions (picos

                            PUBLIC

                            $<$<BOOL:${USE_HDF5_DOUBLE}>:HDF5_DOUBLE>
                            $<$<BOOL:${USE_HDF5_DOUBLE}>:HDF5_FLOAT>
)
target_sources (picos

                PRIVATE

                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fieldSolve.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/outputHDF5.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rfOperator.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/collisionOperator.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpi_main.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PIC.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/units.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initDistribution.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/particleBC.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/types.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parallel_random.cpp>
)
target_precompile_headers (picos

                           PRIVATE

                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initDistribution.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/particleBC.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/types.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fieldSolve.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/outputHDF5.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rfOperator.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/collisionOperator.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initialize.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpi_main.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PIC.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/units.h>>
                           $<$<BOOL:${USE_PCH}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parallel_random.hpp>>
)
target_link_libraries (picos

                       PUBLIC

                       sanitizer
                       armadillo
                       $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
                       MPI::MPI_CXX
                       Bessel::Bessel
                       hdf5::hdf5
                       hdf5::hdf5_cpp
)
set_target_properties (picos PROPERTIES
                       MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}"
                       MPIEXEC_NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
                       MPIEXEC_MAX_NUMPROCS "${MPIEXEC_MAX_NUMPROCS}"
)

#-------------------------------------------------------------------------------
#  Picos Executable
#-------------------------------------------------------------------------------
add_executable (xpicos)
target_sources (xpicos

                PRIVATE

                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/main.cpp>
)
target_link_libraries (xpicos

                       PRIVATE

                       picos
)
