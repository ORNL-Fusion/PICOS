#-------------------------------------------------------------------------------
#  Picos Executable
#-------------------------------------------------------------------------------
add_executable (xunit_tests)
target_sources (xunit_tests

                PRIVATE

                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/unit_tests.cpp>
)
target_link_libraries (xunit_tests

                       PRIVATE

                       picos
)

set (PICOS_TEST_MPI_RANKS 2)
add_test (NAME picos_unit_tests
          COMMAND $<TARGET_PROPERTY:picos,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:picos,MPIEXEC_NUMPROC_FLAG> ${PICOS_TEST_MPI_RANKS} xunit_tests
)
set_tests_properties (picos_unit_tests
                      PROPERTIES
                      PROCESSORS ${PICOS_TEST_MPI_RANKS}
                      TIMEOUT 60
                      ENVIRONMENT "OMP_NUM_THREADS=2"
)
if (${USE_PCH})
    target_precompile_headers (xunit_tests REUSE_FROM picos)
endif ()
