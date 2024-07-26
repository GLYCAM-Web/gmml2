#ifndef INCLUDES_CODEUTILS_METROPOLISCRITERION_HPP
#define INCLUDES_CODEUTILS_METROPOLISCRITERION_HPP

#include <cmath>

namespace monte_carlo
{
    inline bool accept_via_metropolis_criterion(double acceptance, double value)
    {
        return (value < 0) || (std::exp(-value) > acceptance);
    }
} // namespace monte_carlo
#endif
