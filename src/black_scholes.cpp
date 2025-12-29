#include <black_scholes.hpp>
#include <random_number_generator.hpp>
#include <special_functions.hpp>
#include <iterator>

namespace cltvt
{
    BlackScholes::BlackScholes(
        const double discount_rate,
        const double repo_rate,
        const double volatility,
        const double init_level
    )
        :
        m_discount_rate(discount_rate),
        m_repo_rate(repo_rate),
        m_volatility(volatility),
        m_init_level(init_level)
    {
        ASSERT(m_volatility > 1e-12, "volatility must be positive");
        ASSERT(m_init_level > 1e-12, "init_level must be positive");
    }

    BlackScholesPtr BlackScholes::create(
        const double discount_rate,
        const double repo_rate,
        const double volatility,
        const double init_level
    )
    {
        return BlackScholesPtr(new BlackScholes(discount_rate, repo_rate, volatility, init_level));
    }

    double BlackScholes::discount_rate() const
    {
        return m_discount_rate;
    }

    double BlackScholes::repo_rate() const
    {
        return m_repo_rate;
    }

    double BlackScholes::volatility() const
    {
        return m_volatility;
    }

    double BlackScholes::init_level() const
    {
        return m_init_level;
    }

    double BlackScholes::get_call_price(const double strike, const double tenor) const
    {
        const double forward = m_init_level * std::exp((m_discount_rate - m_repo_rate) * tenor);
        const double discount_factor = std::exp(-m_discount_rate * tenor);
        const double total_vol = m_volatility * std::sqrt(tenor);
        const double d1 = std::log(forward / strike) / total_vol + 0.5 * total_vol;
        const double d2 = d1 - total_vol;
        return discount_factor * (forward * normal_cdf(d1) - strike * normal_cdf(d2));
    }

    double BlackScholes::get_put_price(const double strike, const double tenor) const
    {
        const double forward = m_init_level * std::exp((m_discount_rate - m_repo_rate) * tenor);
        const double discount_factor = std::exp(-m_discount_rate * tenor);
        const double total_vol = m_volatility * std::sqrt(tenor);
        const double d1 = std::log(strike / forward) / total_vol + 0.5 * total_vol;
        const double d2 = d1 - total_vol;
        return discount_factor * (strike * normal_cdf(-d2) - forward * normal_cdf(-d1));
    }

    double BlackScholes::get_vega(const double strike, const double tenor) const
    {
        const double vol_bump = 0.001;
        const BlackScholes bs_up(m_discount_rate, m_repo_rate, m_volatility + vol_bump, m_init_level);
        const double price = get_call_price(strike, tenor);
        const double price_up = bs_up.get_call_price(strike, tenor);
        return (price_up - price) / vol_bump;
    }

    double BlackScholes::get_call_rho(const double strike, const double tenor) const
    {
        const double repo_bump = 0.01 * m_repo_rate;
        const BlackScholes bs_bumped(m_discount_rate, m_repo_rate + repo_bump, m_volatility, m_init_level);
        const double price = get_call_price(strike, tenor);
        const double price_bumped = bs_bumped.get_call_price(strike, tenor);
        return (price - price_bumped) / repo_bump;
    }

    double BlackScholes::get_put_rho(const double strike, const double tenor) const
    {
        const double repo_bump = 0.01 * m_repo_rate;
        const BlackScholes bs_bumped(m_discount_rate, m_repo_rate + repo_bump, m_volatility, m_init_level);
        const double price = get_put_price(strike, tenor);
        const double price_bumped = bs_bumped.get_call_price(strike, tenor);
        return (price - price_bumped) / repo_bump;
    }

    void BlackScholes::populate_path(
        std::vector<double>& stock_path, 
        const std::vector<double>& dtimes,
        const std::vector<double>& random_normals
    ) const
    {
        ASSERT(dtimes.size() == random_normals.size(), "dtimes and  random_normals must have same size");
        stock_path.resize(0);
        stock_path.reserve(random_normals.size() + 1);
        std::vector<double>::const_iterator it_dt = dtimes.cbegin();
        std::vector<double>::const_iterator it_rn = random_normals.cbegin();
        const double rho = m_discount_rate - m_repo_rate;
        double lev = m_init_level;
        stock_path.push_back(lev);
        for (; it_rn != random_normals.cend(); ++it_dt, ++it_rn)
        {
            const double dt = *it_dt;
            const double z = *it_rn;
            lev *= std::exp((rho - 0.5 * m_volatility * m_volatility) * dt + m_volatility * std::sqrt(dt) * z);
            stock_path.push_back(lev);
        }
    }

    void BlackScholes::simulate_stock_levels(
        std::vector<double>& stock_levels,
        const std::vector<double>& dtimes,
        const size_t num_samples,
        const size_t seed
    ) const
    {
        stock_levels.resize(0);
        stock_levels.reserve(num_samples);
        StandardNormalGenerator rng(seed);
        std::vector<double> random_normals;
        std::vector<double> stock_path;
        for (size_t i = 0; i < num_samples; ++i)
        {
            rng.populate_standard_normals(random_normals, dtimes.size());
            populate_path(stock_path, dtimes, random_normals);
            stock_levels.push_back(stock_path.back());
        }
    }
}