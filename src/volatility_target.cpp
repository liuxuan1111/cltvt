#include <volatility_target.hpp>
#include <random_number_generator.hpp>
#include <cmath>

namespace cltvt
{
    VolatilityTarget::VolatilityTarget(
        const BlackScholesPtr& sde,
        const double lamb,
        const size_t num_time_steps,
        const double target_volatility,
        const double tenor,
        const double init_var,
        const double init_level
    )
        :
        m_sde(sde),
        m_lamb(lamb),
        m_num_time_steps(num_time_steps),
        m_target_vol(target_volatility),
        m_tenor(tenor),
        m_init_var(init_var),
        m_init_level(init_level),
        m_dt(tenor / num_time_steps)
    {
        ASSERT(m_lamb > 0.0 && m_lamb < 1.0, "0.0 < lamb < 1.0 must be true (lamb=" + std::to_string(m_lamb) + ")");
        ASSERT(m_target_vol > 0.0, "target_volatility must be positive");
        ASSERT(m_tenor > 0.0, "tenor must be positive");
        ASSERT(m_num_time_steps > 1, "num_time_steps > 1 must be true");
        ASSERT(m_init_var > 1e-12, "init_var must be positive");
        ASSERT(m_init_level > 1e-12, "init_level must be positive");
    }

    double VolatilityTarget::lambda() const
    {
        return m_lamb;
    }

    double VolatilityTarget::target_volatility() const
    {
        return m_target_vol;
    }

    double VolatilityTarget::tenor() const
    {
        return m_tenor;
    }

    double VolatilityTarget::init_var() const
    {
        return m_init_var;
    }

    double VolatilityTarget::init_level() const
    {
        return m_init_level;
    }

    size_t VolatilityTarget::num_time_steps() const
    {
        return m_num_time_steps;
    }

    double VolatilityTarget::rebalance_time_step() const
    {
        return m_dt;
    }

    double VolatilityTarget::compute_vt_level(const std::vector<double>& stock_path) const
    {
        ASSERT(stock_path.size() == m_num_time_steps + 1, "stock_path size should be num_time_step + 1");
        double level = m_init_level;
        double var = m_init_var;
        for (size_t i = 1; i < stock_path.size(); ++i)
        {
            const double ret = stock_path[i] / stock_path[i - 1] - 1.0;
            const double w = m_target_vol / std::sqrt(var);
            level *= 1.0 + (1.0 - w) * m_sde->discount_rate() * m_dt + w * ret;
            var = m_lamb * var + (1.0 - m_lamb) * ret * ret / m_dt;
        }
        return level;
    }

    void VolatilityTarget::simulate_vt_levels(std::vector<double>& vt_levels, const size_t num_samples, const size_t seed) const
    {
        vt_levels.resize(0);
        vt_levels.reserve(num_samples);
        StandardNormalGenerator rng(seed);
        const std::vector<double> dtimes(m_num_time_steps, m_dt);
        std::vector<double> random_normals;
        std::vector<double> stock_path;
        for (size_t i = 0; i < num_samples; ++i)
        {
            rng.populate_standard_normals(random_normals, m_num_time_steps);
            m_sde->populate_path(stock_path, dtimes, random_normals);
            const double level = compute_vt_level(stock_path);
            vt_levels.push_back(level);
        }
    }
}