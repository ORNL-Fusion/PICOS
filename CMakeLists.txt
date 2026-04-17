cmake_minimum_required (VERSION 3.14)

project (picos CXX)

#-------------------------------------------------------------------------------
#  Build Options
#-------------------------------------------------------------------------------
option (USE_PCH "Enable the use of precompiled headers" ON)
option (USE_HDF5_DOUBLE "Save hdf5 files as double. When off uses float" OFF)

#-------------------------------------------------------------------------------
#  Set the cmake module path.
#-------------------------------------------------------------------------------
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#-------------------------------------------------------------------------------
#  Setup build types
#-------------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()
set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
              Debug
              Release
              MinSizeRel
              RelWithDebInfo
)

#-------------------------------------------------------------------------------
#  Sanitizer options
#-------------------------------------------------------------------------------
add_library (sanitizer INTERFACE)
target_compile_options (sanitizer
                        INTERFACE
                        $<$<CONFIG:Sanitized>:-g>
)

macro (register_sanitizer_option name)
    string (TOUPPER ${name} upper_name)

    option (SANITIZE_${upper_name} "Enable the ${name} sanitizer" OFF)

    target_compile_options (sanitizer
                            INTERFACE
                            $<$<BOOL:${SANITIZE_${upper_name}}>:-fsanitize=${name}>
    )
    target_link_options (sanitizer
                         INTERFACE
                         $<$<BOOL:${SANITIZE_${upper_name}}>:-fsanitize=${name}>
    )
endmacro ()

register_sanitizer_option (address)
register_sanitizer_option (leak)
register_sanitizer_option (memory)
register_sanitizer_option (thread)
register_sanitizer_option (undefined)
register_sanitizer_option (float-divide-by-zero)

#-------------------------------------------------------------------------------
#  Setup access method.
#-------------------------------------------------------------------------------
option (USE_SSH "Use ssh to access git repos" OFF)
if (${USE_SSH})
    set (URL_PROTO "git@")
    set (URL_SEP ":")
else ()
    set (URL_PROTO "https://")
    set (URL_SEP "/")
endif ()

#-------------------------------------------------------------------------------
#  Define a macro to register new projects.
#-------------------------------------------------------------------------------
include (FetchContent)
find_package (Git)

function (register_project reg_name dir url default_tag)
    set (BUILD_TAG_${dir} ${default_tag} CACHE STRING "Name of the tag to checkout.")

    FetchContent_Declare (
        ${reg_name}
        GIT_REPOSITORY ${url}
        GIT_TAG origin/${BUILD_TAG_${dir}}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${dir}
    )

    FetchContent_MakeAvailable (${reg_name})

    if (${GIT_FOUND})
#  By default cmake clones projects in a headless state. After the repo is
#  cloned checkout the requested tag so repo is in a working state.
        execute_process (
            COMMAND ${GIT_EXECUTABLE} checkout ${BUILD_TAG_${dir}}
            WORKING_DIRECTORY ${${reg_name}_SOURCE_DIR}
        )

#  Add a taraget to pull the latest version before building. Note dependency is
#  registered in the sub project CMakeList.txt. Not sure how this should handle
#  multiple targets in a project yet. Name must match the subproject pull_
#  dependency.
        add_custom_target (
            pull_${reg_name}
            ALL
            COMMAND ${GIT_EXECUTABLE} pull
            WORKING_DIRECTORY ${${reg_name}_SOURCE_DIR}
        )
    endif ()
endfunction ()

register_project (armadillo
                  ARMADILLO
                  ${URL_PROTO}gitlab.com${URL_SEP}conradsnicta/armadillo-code.git
                  15.2.x
)
add_dependencies (armadillo pull_armadillo)

#-------------------------------------------------------------------------------
#  External dependencies.
#-------------------------------------------------------------------------------
find_package (OpenMP)
find_package (MPI REQUIRED)
find_package (Bessel)
find_package (HDF5 MODULE REQUIRED COMPONENTS C CXX)

#-------------------------------------------------------------------------------
#  Setup testing
#-------------------------------------------------------------------------------
enable_testing ()

#-------------------------------------------------------------------------------
#  Setup targets
#-------------------------------------------------------------------------------
add_subdirectory (picosFILES)
