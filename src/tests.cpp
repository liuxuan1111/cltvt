#include <tests.hpp>
#include <volatility_target.hpp>
#include <special_functions.hpp>
#include <integration.hpp>
#include <algorithm>
#include <fstream>

namespace cltvt
{
    double sample_mean(const std::vector<double>& vec)
    {
        if (vec.empty())
            return 0.0;

        double s = 0.0;
        for (double x : vec)
            s += x;
        return s / vec.size();
    }

    double sample_std(const std::vector<double>& vec)
    {
        if (vec.empty())
            return std::nan("");

        const double m = sample_mean(vec);
        double s = 0.0;
        for (double x : vec)
            s += x * x;
        return std::sqrt(s / vec.size() - m * m);
    }

    double multiplier_U(const double lambda)
    {
        auto f = [lambda](const double t) {
            return 1.0 / std::sqrt(q_pochhammer(-t * t, lambda));
        };

        return std::sqrt(2.0 / PI / (1.0 - lambda)) * integrate(f, 0.0, 20.0);
    }

    double multiplier_V(const double lambda)
    {
        auto f = [lambda](const double t) {
            return 1.0 / std::sqrt(q_pochhammer(-t, lambda));
        };

        return 0.5 / (1.0 - lambda) * integrate(f, 0.0, 20.0);
    }

    #define BEGIN_TEST(test_name) std::cout << "Running " + std::string(test_name) + "..." << std::endl;
    #define END_TEST(test_name) std::cout << "Test results saved to tests/" + std::string(test_name) + ".csv\n" << std::endl;

    void test_multiplier_U_bounds()
    {
        BEGIN_TEST("test_multiplier_U_bounds");

        std::vector<double> lamb_vec;
        double lamb = 0.7;
        while (lamb < 1.0)
        {
            lamb_vec.push_back(lamb);
            lamb += 0.02;
        }

        std::vector<double> vals;
        std::vector<double> upper_bounds;
        std::vector<double> lower_bounds;
        for (const double lamb : lamb_vec)
        {
            const double val = multiplier_U(lamb);
            const double one_by_lamb = 1.0 / lamb;
            const double upper_bound = std::sqrt(std::pow(one_by_lamb, 1.25) * std::log(one_by_lamb) / (one_by_lamb - 1.0))
                / (1.0 - std::exp(-2.0 * PI * PI / std::log(one_by_lamb)));
            const double lower_bound = std::sqrt(std::pow(one_by_lamb, 1.2) * std::log(one_by_lamb) / (one_by_lamb - 1.0));
            vals.push_back(val);
            upper_bounds.push_back(upper_bound);
            lower_bounds.push_back(lower_bound);
            std::cout << "lambda=" << lamb << ", U=" << val << ", upper_bound=" << upper_bound << ", lower_bound=" << lower_bound << std::endl;
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_multiplier_U_bounds.csv");
        outfile << "lambda,U,upper_bound,lower_bound\n";
        for (size_t i = 0; i < lamb_vec.size(); ++i)
        {
            outfile << lamb_vec[i] << "," << vals[i] << "," << upper_bounds[i] << "," << lower_bounds[i] << "\n";
        }
        outfile.close();

        END_TEST("test_multiplier_U_bounds");
    }

    void test_multiplier_V_bounds()
    {
        BEGIN_TEST("test_multiplier_V_bounds");

        std::vector<double> lamb_vec;
        double lamb = 0.7;
        while (lamb < 1.0)
        {
            lamb_vec.push_back(lamb);
            lamb += 0.02;
        }

        std::vector<double> vals;
        std::vector<double> upper_bounds;
        std::vector<double> lower_bounds;
        for (const double lamb : lamb_vec)
        {
            const double val = multiplier_V(lamb);
            const double one_by_lamb = 1.0 / lamb;
            const double upper_bound = std::pow(one_by_lamb, 1.5) * std::log(one_by_lamb) / (one_by_lamb - 1.0);
            const double lower_bound = std::pow(one_by_lamb, 1.45) * std::log(one_by_lamb) / (one_by_lamb - 1.0);
            vals.push_back(val);
            upper_bounds.push_back(upper_bound);
            lower_bounds.push_back(lower_bound);
            std::cout << "lambda=" << lamb << ", V=" << val << ", upper_bound=" << upper_bound << ", lower_bound=" << lower_bound << std::endl;
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_multiplier_V_bounds.csv");
        outfile << "lambda,V,upper_bound,lower_bound\n";
        for (size_t i = 0; i < lamb_vec.size(); ++i)
        {
            outfile << lamb_vec[i] << "," << vals[i] << "," << upper_bounds[i] << "," << lower_bounds[i] << "\n";
        }
        outfile.close();

        END_TEST("test_multiplier_V_bounds");
    }

    void test_vt_volatility(const size_t num_samples)
    {
        BEGIN_TEST("test_vt_volatility");

        const double discount_rate = 0.05;
        const double rho = 0.03;
        const double volatility = 0.5;
        const double target_volatility = 0.2;
        const double tenor = 1.0;
        const double init_var = 0.02;
        const double init_stock_level = 1.0;
        const double init_vt_level = 1.0;
        const double repo_rate = discount_rate - rho;

        const std::vector<size_t> num_time_steps { 1000, 2000, 5000, 10000, 50000 };
        std::vector<double> lamb_vec;
        double lamb = 0.7;
        while (lamb < 1.0)
        {
            lamb_vec.push_back(lamb);
            lamb += 0.05;
        }
        if (lamb_vec.back() < 0.97)
            lamb_vec.push_back(0.97);

        const BlackScholesPtr sde = BlackScholes::create(discount_rate, repo_rate, volatility, init_stock_level);
        std::vector<std::vector<double>> vt_vol_array(0);
        std::vector<std::vector<double>> limit_vol_array(0);
        vt_vol_array.reserve(num_time_steps.size());
        limit_vol_array.reserve(num_time_steps.size());
        std::vector<double> log_vt_levels;
        for (const size_t num_steps : num_time_steps)
        {
            std::vector<double> vols;
            std::vector<double> limit_vols;
            vols.reserve(lamb_vec.size());
            limit_vols.reserve(lamb_vec.size());
            for (const double lamb : lamb_vec)
            {
                VolatilityTarget vt(sde, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
                vt.simulate_vt_levels(log_vt_levels, num_samples);
                for (std::vector<double>::iterator itv = log_vt_levels.begin(); itv != log_vt_levels.end(); ++itv)
                    *itv = std::log(*itv / vt.init_level());
                const double vol = sample_std(log_vt_levels) / std::sqrt(tenor);
                const double limit_vol = target_volatility * std::sqrt(multiplier_V(lamb));
                vols.push_back(vol);
                limit_vols.push_back(limit_vol);
                std::cout << "N=" << num_steps << ", lamb=" << lamb << ", vt_vol=" << vol << ", limit_vol=" << limit_vol << std::endl;
            }
            vt_vol_array.push_back(vols);
            limit_vol_array.push_back(limit_vols);
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_vt_volatility.csv");
        outfile << "N,lambda,vt_vol,limit_vol\n";
        for (size_t i_num_step = 0; i_num_step < num_time_steps.size(); ++i_num_step)
        {
            for (size_t i_lamb = 0; i_lamb < lamb_vec.size(); ++i_lamb)
            {
                outfile << num_time_steps[i_num_step] << "," << lamb_vec[i_lamb] << "," 
                    << vt_vol_array[i_num_step][i_lamb] << "," << limit_vol_array[i_num_step][i_lamb] << "\n";
            }
        }
        outfile.close();

        END_TEST("test_vt_volatility");
    }

    void test_vt_volatility_simultaneous_limit(const size_t num_samples)
    {
        BEGIN_TEST("test_vt_volatility_simultaneous_limit");

        const double discount_rate = 0.05;
        const double rho = 0.03;
        const double volatility = 0.5;
        const double target_volatility = 0.2;
        const double tenor = 1.0;
        const double init_var = volatility * volatility;
        const double init_stock_level = 1.0;
        const double init_vt_level = 1.0;
        const double repo_rate = discount_rate - rho;

        const std::vector<size_t> num_time_steps { 1000, 2000, 5000, 10000, 50000 };
        std::vector<double> lamb_vec;
        double lamb = 0.7;
        while (lamb < 1.0)
        {
            lamb_vec.push_back(lamb);
            lamb += 0.05;
        }
        if (lamb_vec.back() < 0.97)
            lamb_vec.push_back(0.97);
        if (lamb_vec.back() < 0.99)
            lamb_vec.push_back(0.99);

        const BlackScholesPtr sde = BlackScholes::create(discount_rate, repo_rate, volatility, init_stock_level);
        std::vector<std::vector<double>> vt_vol_array(0);
        vt_vol_array.reserve(num_time_steps.size());
        std::vector<double> log_vt_levels;
        for (const size_t num_steps : num_time_steps)
        {
            std::vector<double> vols;
            std::vector<double> lim_vols;
            vols.reserve(lamb_vec.size());
            for (const double lamb : lamb_vec)
            {
                VolatilityTarget vt(sde, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
                vt.simulate_vt_levels(log_vt_levels, num_samples);
                for (std::vector<double>::iterator itv = log_vt_levels.begin(); itv != log_vt_levels.end(); ++itv)
                    *itv = std::log(*itv / vt.init_level());
                const double vol = sample_std(log_vt_levels) / std::sqrt(tenor);
                vols.push_back(vol);
                std::cout << "N=" << num_steps << ", lamb=" << lamb << ", vt_vol=" << vol << ", target_vol=" << target_volatility << std::endl;
            }
            vt_vol_array.push_back(vols);
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_vt_volatility_simultaneous_limit.csv");
        outfile << "N,lambda,vt_vol,target_vol\n";
        for (size_t i_num_step = 0; i_num_step < num_time_steps.size(); ++i_num_step)
        {
            for (size_t i_lamb = 0; i_lamb < lamb_vec.size(); ++i_lamb)
            {
                outfile << num_time_steps[i_num_step] << "," << lamb_vec[i_lamb] << ","
                    << vt_vol_array[i_num_step][i_lamb] << "," << target_volatility << "\n";
            }
        }
        outfile.close();

        END_TEST("test_vt_volatility_simultaneous_limit");
    }

    void test_vt_volatility_limit_along_path_1(const size_t num_samples)
    {
        BEGIN_TEST("test_vt_volatility_limit_along_path_1");

        const double discount_rate = 0.05;
        const double rho = 0.03;
        const double volatility = 0.5;
        const double target_volatility = 0.2;
        const double tenor = 1.0;
        const double init_var = 0.02;
        const double init_stock_level = 1.0;
        const double init_vt_level = 1.0;
        const double repo_rate = discount_rate - rho;
        const double kappa = 1.0;

        const std::vector<size_t> num_time_steps { 1000, 2000, 5000, 10000, 50000 };
        std::vector<double> lamb_vec;
        for (const size_t num_steps : num_time_steps)
            lamb_vec.push_back(1.0 - 1.0 / (num_steps * num_steps));

        const BlackScholesPtr sde = BlackScholes::create(discount_rate, repo_rate, volatility, init_stock_level);
        std::vector<double> vt_vols;
        std::vector<double> log_vt_levels;
        for (size_t i = 0; i < num_time_steps.size(); ++i)
        {
            const size_t num_steps = num_time_steps[i];
            const double lamb = lamb_vec[i];
            VolatilityTarget vt(sde, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
            vt.simulate_vt_levels(log_vt_levels, num_samples);
            for (std::vector<double>::iterator itv = log_vt_levels.begin(); itv != log_vt_levels.end(); ++itv)
                *itv = std::log(*itv / vt.init_level());
            const double vol = sample_std(log_vt_levels) / std::sqrt(tenor);
            vt_vols.push_back(vol);
            std::cout << "N=" << num_steps << ", lamb=" << lamb << ", v0=" << init_var << ", stock_vol=" << volatility 
                << ", target_vol=" << target_volatility << ", vt_vol=" << vol << std::endl;
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_vt_volatility_limit_along_path_1.csv");
        outfile << "N,lambda,v0,stock_vol,target_vol,vt_vol\n";
        for (size_t i = 0; i < num_time_steps.size(); ++i)
            outfile << num_time_steps[i] << "," << lamb_vec[i] << "," << init_var << "," << volatility 
                << "," << target_volatility << "," << vt_vols[i] << "\n";
        outfile.close();

        END_TEST("test_vt_volatility_limit_along_path_1");
    }

    void test_vt_volatility_limit_along_path_2(const size_t num_samples)
    {
        BEGIN_TEST("test_vt_volatility_limit_along_path_2");

        const double discount_rate = 0.05;
        const double rho = 0.03;
        const double volatility = 0.5;
        const double target_volatility = 0.2;
        const double tenor = 1.0;
        const double init_var = 0.02;
        const double init_stock_level = 1.0;
        const double init_vt_level = 1.0;
        const double repo_rate = discount_rate - rho;
        const double kappa = 1.0;

        const std::vector<size_t> num_time_steps { 1000, 2000, 5000, 10000, 50000 };
        std::vector<double> lamb_vec;
        for (const size_t num_steps : num_time_steps)
            lamb_vec.push_back(1.0 - std::log(num_steps) / std::sqrt(num_steps));

        const BlackScholesPtr sde = BlackScholes::create(discount_rate, repo_rate, volatility, init_stock_level);
        std::vector<double> vt_vols;
        std::vector<double> log_vt_levels;
        for (size_t i = 0; i < num_time_steps.size(); ++i)
        {
            const size_t num_steps = num_time_steps[i];
            const double lamb = lamb_vec[i];
            VolatilityTarget vt(sde, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
            vt.simulate_vt_levels(log_vt_levels, num_samples);
            for (std::vector<double>::iterator itv = log_vt_levels.begin(); itv != log_vt_levels.end(); ++itv)
                *itv = std::log(*itv / vt.init_level());
            const double vol = sample_std(log_vt_levels) / std::sqrt(tenor);
            vt_vols.push_back(vol);
            std::cout << "N=" << num_steps << ", lamb=" << lamb << ", v0=" << init_var << ", stock_vol=" << volatility
                << ", target_vol=" << target_volatility << ", vt_vol=" << vol << std::endl;
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_vt_volatility_limit_along_path_2.csv");
        outfile << "N,lambda,v0,stock_vol,target_vol,vt_vol\n";
        for (size_t i = 0; i < num_time_steps.size(); ++i)
            outfile << num_time_steps[i] << "," << lamb_vec[i] << "," << init_var << "," << volatility
                << "," << target_volatility << "," << vt_vols[i] << "\n";
        outfile.close();

        END_TEST("test_vt_volatility_limit_along_path_2");
    }

    void test_vt_pricing(const size_t num_samples)
    {
        BEGIN_TEST("test_vt_pricing");

        const double discount_rate = 0.05;
        const double rho = 0.03;
        const double volatility = 0.5;
        const double target_volatility = 0.2;
        const double tenor = 1.0;
        const double init_var = 0.02;
        const double init_stock_level = 1.0;
        const double init_vt_level = 1.0;
        const double repo_rate = discount_rate - rho;

        const std::vector<size_t> num_time_steps { 1000, 2000, 5000, 10000, 50000 };
        std::vector<double> lamb_vec;
        double lamb = 0.7;
        while (lamb < 1.0)
        {
            lamb_vec.push_back(lamb);
            lamb += 0.05;
        }
        if (lamb_vec.back() < 0.97)
            lamb_vec.push_back(0.97);

        const BlackScholesPtr sde = BlackScholes::create(discount_rate, repo_rate, volatility, init_stock_level);
        std::vector<std::vector<double>> mc_vt_price_array(0);
        std::vector<std::vector<double>> bs_limit_price_array(0);
        mc_vt_price_array.reserve(num_time_steps.size());
        bs_limit_price_array.reserve(num_time_steps.size());
        std::vector<double> vt_levels;
        std::vector<double> stock_levels;
        for (const size_t num_steps : num_time_steps)
        {
            std::vector<double> mc_vt_prices;
            std::vector<double> bs_limit_prices;
            mc_vt_prices.reserve(lamb_vec.size());
            bs_limit_prices.reserve(lamb_vec.size());
            for (const double lamb : lamb_vec)
            {
                VolatilityTarget vt(sde, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
                vt.simulate_vt_levels(vt_levels, num_samples);
                double mc_vt_price = 0.0;
                for (const double level : vt_levels)
                {
                    mc_vt_price += std::max(level - vt.init_level(), 0.0);
                }
                mc_vt_price *= std::exp(-discount_rate * tenor);
                mc_vt_price /= vt_levels.size();

                const double limit_vol = target_volatility * std::sqrt(multiplier_V(lamb));
                const double limit_repo = multiplier_U(lamb) * target_volatility / volatility * repo_rate;
                const BlackScholesPtr limit_bs = BlackScholes::create(discount_rate, limit_repo, limit_vol, vt.init_level());
                const double bs_limit_price = limit_bs->get_call_price(vt.init_level(), tenor);

                mc_vt_prices.push_back(mc_vt_price);
                bs_limit_prices.push_back(bs_limit_price);
                std::cout << "N=" << num_steps << ", lamb=" << lamb << ", mc_vt_price=" 
                    << mc_vt_price << ", bs_limit_price=" << bs_limit_price << std::endl;
            }
            mc_vt_price_array.push_back(mc_vt_prices);
            bs_limit_price_array.push_back(bs_limit_prices);
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_vt_pricing.csv");
        outfile << "N,lambda,mc_vt_price,bs_limit_price\n";
        for (size_t i_num_step = 0; i_num_step < num_time_steps.size(); ++i_num_step)
        {
            for (size_t i_lamb = 0; i_lamb < lamb_vec.size(); ++i_lamb)
            {
                outfile << num_time_steps[i_num_step] << "," << lamb_vec[i_lamb] << "," << mc_vt_price_array[i_num_step][i_lamb] 
                    << "," << bs_limit_price_array[i_num_step][i_lamb] << "\n";
            }
        }
        outfile.close();

        END_TEST("test_vt_pricing");
    }

    void test_vt_vega(const size_t num_samples)
    {
        BEGIN_TEST("test_vt_vega");

        const double discount_rate = 0.05;
        const double rho = 0.03;
        const double volatility = 0.5;
        const double target_volatility = 0.2;
        const double tenor = 1.0;
        const double init_var = 0.02;
        const double init_stock_level = 1.0;
        const double init_vt_level = 1.0;
        const double repo_rate = discount_rate - rho;

        const std::vector<size_t> num_time_steps { 1000, 2000, 5000, 10000, 50000 };
        std::vector<double> lamb_vec;
        double lamb = 0.7;
        while (lamb < 1.0)
        {
            lamb_vec.push_back(lamb);
            lamb += 0.05;
        }
        if (lamb_vec.back() < 0.97)
            lamb_vec.push_back(0.97);

        const double vol_bump = 0.001;
        const BlackScholesPtr sde = BlackScholes::create(discount_rate, repo_rate, volatility, init_stock_level);
        const BlackScholesPtr sde_bumped = BlackScholes::create(discount_rate, repo_rate, volatility + vol_bump, init_stock_level);
        std::vector<std::vector<double>> mc_vt_vega_array(0);
        std::vector<std::vector<double>> bs_limit_vega_array(0);
        mc_vt_vega_array.reserve(num_time_steps.size());
        bs_limit_vega_array.reserve(num_time_steps.size());
        std::vector<double> vt_levels;
        std::vector<double> stock_levels;
        for (const size_t num_steps : num_time_steps)
        {
            std::vector<double> mc_vt_vegas;
            std::vector<double> bs_limit_vegas;
            mc_vt_vegas.reserve(lamb_vec.size());
            bs_limit_vegas.reserve(lamb_vec.size());
            for (const double lamb : lamb_vec)
            {
                VolatilityTarget vt(sde, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
                vt.simulate_vt_levels(vt_levels, num_samples);
                double mc_vt_price = 0.0;
                for (const double level : vt_levels)
                {
                    mc_vt_price += std::max(level - vt.init_level(), 0.0);
                }
                mc_vt_price *= std::exp(-discount_rate * tenor);
                mc_vt_price /= vt_levels.size();

                VolatilityTarget vt_bumped(sde_bumped, lamb, num_steps, target_volatility, tenor, init_var, init_vt_level);
                vt_bumped.simulate_vt_levels(vt_levels, num_samples);
                double mc_vt_price_bumped = 0.0;
                for (const double level : vt_levels)
                {
                    mc_vt_price_bumped += std::max(level - vt.init_level(), 0.0);
                }
                mc_vt_price_bumped *= std::exp(-discount_rate * tenor);
                mc_vt_price_bumped /= vt_levels.size();

                const double mc_vt_vega = (mc_vt_price_bumped - mc_vt_price) / vol_bump;

                const double limit_vol = target_volatility * std::sqrt(multiplier_V(lamb));
                const double limit_repo = multiplier_U(lamb) * target_volatility / volatility * repo_rate;
                const BlackScholesPtr limit_bs = BlackScholes::create(discount_rate, limit_repo, limit_vol, vt.init_level());
                const double bs_limit_rho = limit_bs->get_call_rho(vt.init_level(), tenor);
                const double bs_limit_vega = limit_repo / volatility * bs_limit_rho;

                mc_vt_vegas.push_back(mc_vt_vega);
                bs_limit_vegas.push_back(bs_limit_vega);
                std::cout << "N=" << num_steps << ", lamb=" << lamb << ", mc_vt_vega="
                    << mc_vt_vega << ", bs_limit_vega=" << bs_limit_vega << std::endl;
            }
            mc_vt_vega_array.push_back(mc_vt_vegas);
            bs_limit_vega_array.push_back(bs_limit_vegas);
        }

        std::ofstream outfile;
        outfile.open(root_dir() + "/tests/test_vt_vega.csv");
        outfile << "N,lambda,mc_vt_vega,bs_limit_vega\n";
        for (size_t i_num_step = 0; i_num_step < num_time_steps.size(); ++i_num_step)
        {
            for (size_t i_lamb = 0; i_lamb < lamb_vec.size(); ++i_lamb)
            {
                outfile << num_time_steps[i_num_step] << "," << lamb_vec[i_lamb] << "," << mc_vt_vega_array[i_num_step][i_lamb]
                    << "," << bs_limit_vega_array[i_num_step][i_lamb] << "\n";
            }
        }
        outfile.close();

        END_TEST("test_vt_vega");
    }

}
