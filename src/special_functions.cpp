#include <preliminaries.hpp>
#include <special_functions.hpp>
#include <cmath>

namespace cltvt
{
    double normal_cdf(const double x)
    {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2)));
    }

    double q_pochhammer(const double a, const double q, const int n)
    {
        ASSERT(std::abs(q) < 1, "abs(q) < 1 must be true");

        if (std::abs(a) < 1e-12)
            return 1.0;

        int n_to_use = n;
        if (n < 0)
        {
            const double eps = 1e-8;
            const double abs_q = std::abs(q);
            const double nd = std::log(0.5 * eps * (1.0 - abs_q) / std::abs(a)) / std::log(abs_q);
            n_to_use = (int)std::ceil(nd);
        }

        double s = 0.0;
        double qk = 1.0;
        for (int k = 0; k < n_to_use; ++k)
        {
            s += std::log(1.0 - a * qk);
            qk *= q;
        }
        return std::exp(s);
    }
}