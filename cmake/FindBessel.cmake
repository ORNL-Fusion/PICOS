#[==[
Provides the following variables:

  * `Bessel::Bessel`: A target to use with `target_link_libraries`.
#]==]

include (CheckCXXSourceCompiles)

set (_BESSEL_OLD_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
set (CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++17")
check_cxx_source_compiles ("
#include <cmath>

int main() {
    volatile double x = std::cyl_bessel_j(1,1.0);
    return x == 0.0;
}
" HAS_BESSEL)
set (CMAKE_REQUIRED_FLAGS "${_BESSEL_OLD_CMAKE_REQUIRED_FLAGS}")

add_library (Bessel INTERFACE)

if (HAS_BESSEL)
    target_compile_definitions (Bessel

                                INTERFACE

                                HAS_STD_BESSEL
                                CYL_BESSEL_J=std::cyl_bessel_j
    )
else ()
    register_project (boost_math
                      BOOST_MATH
                      ${URL_PROTO}github.com${URL_SEP}boostorg/math.git
                      develop
    )
    add_dependencies (boost_math pull_boost_math)

    target_link_libraries (Bessel INTERFACE Boost::math)
    target_compile_definitions (Bessel

                                INTERFACE

                                CYL_BESSEL_J=boost::math::cyl_bessel_j
    )
endif ()
add_library (Bessel::Bessel ALIAS Bessel)
