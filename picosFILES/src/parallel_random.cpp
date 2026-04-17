#ifdef _OPENMP
#include <omp.h>
#endif

#include "parallel_random.hpp"

//------------------------------------------------------------------------------
///  @brief A random number instance class.
///
///  When not using openmp default to one thread.
///
///  @returns The number of threads used.
//------------------------------------------------------------------------------
size_t picos::random::threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

//------------------------------------------------------------------------------
///  @brief Get the current thread number.
///
///  When not using openmp defaults to zero.
///
///  @note This must be called from inside a parallel section.
///
///  @returns The thread number.
//------------------------------------------------------------------------------
size_t picos::random::thread() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

//------------------------------------------------------------------------------
///  @brief Compute a seed offset for a thread number and comm rank.
///
///  @returns rank*num_thread + thread_number
//------------------------------------------------------------------------------
size_t picos::random::offset(const size_t t) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank*threads() + t;
}

