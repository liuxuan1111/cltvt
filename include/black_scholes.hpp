#pragma once
#include <preliminaries.hpp>
#include <vector>
#include <memory>

namespace cltvt
{
    class BlackScholes;
    typedef std::shared_ptr<BlackScholes> BlackScholesPtr;

    class BlackScholes
    {
    public:
        BlackScholes(
            const double discount_rate,
            const double repo_rate,
            const double volatility,
            const double init_level = 1.0
        );

        static BlackScholesPtr create(
            const double discount_rate,
            const double repo_rate,
            const double volatility,
            const double init_level = 1.0
        );

        double discount_rate() const;

        double repo_rate() const;

        double volatility() const;

        double init_level() const;

        double get_call_price(const double strike, const double tenor) const;

        double get_put_price(const double strike, const double tenor) const;

        double get_vega(const double strike, const double tenor) const;

        double get_call_rho(const double strike, const double tenor) const;

        double get_put_rho(const double strike, const double tenor) const;

        void populate_path(
            std::vector<double>& stock_path, 
            const std::vector<double>& dtimes,
            const std::vector<double>& random_normals
        ) const;

        void simulate_stock_levels(
            std::vector<double>& stock_levels,
            const std::vector<double>& dtimes,
            const size_t num_samples,
            const size_t seed = DEFAULT_RNG_SEED
        ) const;

    private:
        double m_discount_rate;
        double m_repo_rate;
        double m_volatility;
        double m_init_level;
    };
}