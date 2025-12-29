#include <random_number_generator.hpp>

namespace cltvt
{
    StandardNormalGenerator::StandardNormalGenerator(const size_t seed) : m_seed(seed)
    {
        m_rng = std::mt19937_64(m_seed);
        m_dist = std::normal_distribution<double>(0.0, 1.0);
    }

    void StandardNormalGenerator::populate_standard_normals(std::vector<double>& rn_out, const size_t size)
    {
        rn_out.resize(0);
        rn_out.reserve(size);
        for (size_t i = 0; i < size; ++i)
            rn_out.push_back(m_dist(m_rng));
    }
}