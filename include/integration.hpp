#pragma once
#include <preliminaries.hpp>

namespace cltvt
{
    double integrate(const Function& f, const double a, const double b, const size_t N = 5000);
}