#pragma once

namespace cltvt
{
    void test_multiplier_U_bounds();

    void test_multiplier_V_bounds();

    void test_vt_volatility(const size_t num_samples = 100000);

    void test_vt_volatility_simultaneous_limit(const size_t num_samples = 100000);

    void test_vt_volatility_limit_along_path_1(const size_t num_samples = 100000);

    void test_vt_volatility_limit_along_path_2(const size_t num_samples = 100000);

    void test_vt_pricing(const size_t num_samples = 100000);

    void test_vt_vega(const size_t num_samples = 100000);

}
