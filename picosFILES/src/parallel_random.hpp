#ifndef PARALLEL_RANDOM_HPP
#define PARALLEL_RANDOM_HPP

#include <mpi.h>

#include <cmath>
#include <concepts>
#include <random>
#include <vector>

/// Name space for PICOS.
namespace picos {
/// Name space for random numbers.
    namespace random {
/// Distribution concept currently only suports uniform real distributions.
        template<typename T>
        concept distribution = std::same_as<T, std::uniform_real_distribution<float>> ||
                               std::same_as<T, std::uniform_real_distribution<double>> ||
                               std::same_as<T, std::uniform_int_distribution<short>>;

///  Base type concept for random distributions.
        template<typename T>
        concept base = std::floating_point<T> ||
                       std::same_as<T, short>;

//------------------------------------------------------------------------------
///  @brief Get the number of threads used.
///
///  When not using openmp default to one thread.
///
///  @returns The number of threads used.
//------------------------------------------------------------------------------
        size_t threads();

//------------------------------------------------------------------------------
///  @brief Get the current thread number.
///
///  When not using openmp defaults to zero.
///
///  @note This must be called from inside a parallel section.
///
///  @returns The thread number.
//------------------------------------------------------------------------------
        size_t thread();

//------------------------------------------------------------------------------
///  @brief Compute a seed offset for a thread number and comm rank.
///
///  @returns rank*num_thread + thread_number
//------------------------------------------------------------------------------
        size_t offset(const size_t t);

//------------------------------------------------------------------------------
///  @brief A random number instance class.
///
///  @tparam T    Base type.
///  @tparam D    Distribution type.
///  @tparam LOW  Lower bound of the distribution.
///  @tparam HIGH Upper bound of the distribution.
//------------------------------------------------------------------------------
        template<base T, distribution D, T LOW, T HIGH>
        class instance {
        private:
///  Random distribution function.
            D dist;
///  Random engine.
            std::mt19937_64 engine;

        public:
//------------------------------------------------------------------------------
///  @brief An instance constructor.
///
///  When not using openmp default to one thread.
///
///  @param[in] seed   Seed value for the random engine.
///  @param[in] offser Seed offset to the seed.
//------------------------------------------------------------------------------
            instance(const std::uint_fast64_t seed,
                     const std::uint_fast64_t offset) :
            engine(seed + offset), dist(LOW, HIGH) {}

//------------------------------------------------------------------------------
///  @brief Call operator for the distribution.
///
///  When not using openmp defaults to zero.
///
///  @returns A random number from the distribution.
//------------------------------------------------------------------------------
            T operator()() {
                return dist(engine);
            }
        };

//------------------------------------------------------------------------------
///  @brief Factory function to construct the parallel random instances.
///
///  @param[in] seed The inital seed.
///  @returns A vector of parallel instances.
//------------------------------------------------------------------------------
        template<base T, distribution D, T LOW, T HIGH>
        std::vector<instance<T, D, LOW, HIGH>>
        instances(std::uint_fast64_t seed) {
            MPI_Bcast(&seed, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
            std::vector<instance<T, D, LOW, HIGH>> rands;
            for (size_t i = 0, ie = threads(); i < ie; i++) {
                rands.emplace_back(seed, offset(i));
            }
            return rands;
        }

//##############################################################################
//  Unit Tests
//##############################################################################
//------------------------------------------------------------------------------
///  @brief Compute the auto correlation for a specific offset.
///
///  @tparam T Base type.
///
///  @param[in] sequence Random sequence.
///  @param[in] offset   Offset of the correlation.
///  @returns The autocorrelation for a given offset.
//------------------------------------------------------------------------------
        template<std::floating_point T>
        T autocorrelation(const std::vector<T> &sequence,
                          const size_t offset) {
            T result = 0.0;
            for (size_t i = 0, ie = sequence.size() - offset; i < ie; i++) {
                result += sequence[i]*sequence[offset + i];
            }
            return result/static_cast<T> (sequence.size() - offset);
        }

//------------------------------------------------------------------------------
///  @brief MPI type.
///
///  @tparam T Base type.
///
///  @returns The MPI type number.
//------------------------------------------------------------------------------
        template<std::floating_point T>
        MPI_Datatype to_MPI_type() {
            if constexpr (std::is_same<T, float> ()) {
                return MPI_FLOAT;
            } else {
                return MPI_DOUBLE;
            }
        }

//------------------------------------------------------------------------------
///  @brief Run unit tests for random numbers.
///
///  @tparam T Base type.
//------------------------------------------------------------------------------
        template<std::floating_point T>
        void test() {
            constexpr T one = static_cast<T> (1.0);
            constexpr T none = -one;
            typedef instance<T, uniform_real_distribution<T>, none, one> uniform_random;
            std::vector<uniform_random> randoms = instances<T, uniform_real_distribution<T>, none, one> (0);

            const size_t batch_size = 10000;

//  Define a buffer so each rank and thread can operator on batch_size elements.
            int num_ranks;
            MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
            std::vector<T> result(num_ranks*randoms.size()*batch_size, 0.0);

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#pragma omp parallel default(shared)
            for (size_t i = rank*randoms.size()*batch_size, ie = i + batch_size; i < ie; i++) {
                result[i] = randoms[thread()]();
            }

            MPI_Allreduce(MPI_IN_PLACE, result.data(), result.size(),
                          to_MPI_type<T> (), MPI_SUM, MPI_COMM_WORLD);

            const T base = autocorrelation(result, 0)*static_cast<T>(0.05);

            // Representative short, intra-stream, and cross-stream lags are enough to
            // catch correlated generators. Testing thousands of full-sequence lags made
            // this unit test quadratic in the global sample count. With MPI ranks each
            // also inheriting every OpenMP thread, that could consume all CPUs for hours.
            if (rank == 0) {
                const std::vector<size_t> offsets {
                    1, 2, 3, 5, 8, 13, 21, 55,
                    batch_size/2, batch_size - 1, batch_size, batch_size + 1
                };
                for (const size_t i : offsets) {
                    const T test = autocorrelation(result, i);
                    if (std::abs(test) > base) {
                        std::cerr << "Auto correlation failure at lag " << i << ". "
                                  << base << " " << test << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

#endif
