#include <tests.hpp>
#include <special_functions.hpp>
#include <iostream>

int main()
{
    using namespace cltvt;

    test_multiplier_U_bounds();

    test_multiplier_V_bounds();

    test_vt_volatility();

    test_vt_volatility_simultaneous_limit();

    test_vt_volatility_limit_along_path_1();

    test_vt_volatility_limit_along_path_2();

    test_vt_pricing();

    test_vt_vega();

    return 0;
}