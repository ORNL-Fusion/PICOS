add_executable (xpicos)

target_compile_features (xpicos

                         PRIVATE

                         cxx_std_23
)
target_compile_definitions (xpicos

                            PRIVATE

                            $<$<BOOL:${USE_HDF5_DOUBLE}>:HDF5_DOUBLE>
                            $<$<BOOL:${USE_HDF5_DOUBLE}>:HDF5_FLOAT>
)
target_sources (xpicos

                PRIVATE

                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fieldSolve.cpp>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/main.cpp>
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
)

target_precompile_headers (xpicos

                           PRIVATE

                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initDistribution.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/particleBC.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/types.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fieldSolve.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/outputHDF5.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rfOperator.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/collisionOperator.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/initialize.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpi_main.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PIC.h>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/units.h>
)

target_link_libraries (xpicos

                       PRIVATE

                       sanitizer
                       armadillo
                       $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
                       MPI::MPI_CXX
                       Bessel::Bessel
                       hdf5::hdf5
                       hdf5::hdf5_cpp
)
