#include <integration.hpp>

namespace cltvt
{
    double integrate(const Function& f, const double a, const double b, const size_t N)
    {
        ASSERT(a < b, "a < b must be true");
        ASSERT(N > 0, "N > 0 must be true");
        
        const double dx = (b - a) / N;
        double x0 = a;
        double x1 = x0 + dx;
        double s = 0.0;
        while (x1 < b)
        {
            const double x = 0.5 * (x0 + x1);
            s += f(x);
            x0 += dx;
            x1 += dx;
        }
        s *= dx;
        return s;
    }
}