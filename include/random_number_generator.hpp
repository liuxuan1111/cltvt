#pragma once
#include <preliminaries.hpp>
#include <random>

namespace cltvt
{
    class StandardNormalGenerator
    {
    public:
        StandardNormalGenerator(const size_t seed = DEFAULT_RNG_SEED);

        void seed(const size_t seed = DEFAULT_RNG_SEED)
        {
            m_rng.seed(seed);
            m_seed = seed;
        }

        void reset() { seed(m_seed); }
        void populate_standard_normals(std::vector<double>& rn_out, const size_t size);

    private:
        size_t m_seed;
        std::mt19937_64 m_rng;
        std::normal_distribution<double> m_dist;
    };
}