#pragma once
#include <preliminaries.hpp>
#include <black_scholes.hpp>

namespace cltvt
{
    class VolatilityTarget
    {
    public:
        VolatilityTarget(
            const BlackScholesPtr& sde,
            const double lamb,
            const size_t num_time_steps,
            const double target_volatility,
            const double tenor,
            const double init_var,
            const double init_level
        );

        double lambda() const;

        double target_volatility() const;

        double tenor() const;

        double init_var() const;

        double init_level() const;

        size_t num_time_steps() const;

        double rebalance_time_step() const;

        double compute_vt_level(const std::vector<double>& stock_path) const;

        void simulate_vt_levels(std::vector<double>& vt_levels, const size_t num_samples, const size_t seed = DEFAULT_RNG_SEED) const;

    private:
        BlackScholesPtr m_sde;
        double m_lamb;
        double m_target_vol;
        double m_tenor;
        double m_init_var;
        double m_init_level;
        size_t m_num_time_steps;
        double m_dt;
    };
}