#include "includes/CodeUtils/random.hpp"

#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <utility>
#include <numeric>
#include <random>

size_t codeUtils::weightedRandomIndex(pcg32& rng, const std::vector<double>& weights)
{
    double sum    = vectorSum(0.0, weights);
    double target = uniformRandomDoubleWithinRange(rng, 0.0, sum);
    double accum  = 0.0;
    for (size_t n = 0; n < weights.size(); n++)
    {
        accum += weights[n];
        if (target < accum)
        {
            return n;
        }
    }
    return weights.size() - 1;
}

std::vector<size_t> codeUtils::weightedRandomOrder(pcg32& rng, std::vector<double> weights)
{
    std::vector<size_t> result;
    result.reserve(weights.size());
    std::vector<size_t> indices = indexVector(weights);

    while (!indices.empty())
    {
        double sum       = std::accumulate(weights.begin(), weights.end(), 0.0);
        double threshold = uniformRandomDoubleWithinRange(rng, 0.0, sum);
        for (size_t n = 0; n < indices.size(); n++)
        {
            threshold -= weights[n];
            if (n == indices.size() - 1 || threshold < 0.0)
            {
                result.push_back(indices[n]);
                indices.erase(indices.begin() + n);
                weights.erase(weights.begin() + n);
                break;
            }
        }
    }

    return result;
}

double codeUtils::normalDistributionRandomDoubleWithCutoff(pcg32& rng, double lowerCutoff, double upperCutoff)
{
    double num = normalDistributionRandomDouble(rng);
    // number is within cutoff
    if ((num >= lowerCutoff) && (num <= upperCutoff))
    {
        return num;
    }
    // number outside of requested range, fall back to uniform distribution
    else
    {
        return uniformRandomDoubleWithinRange(rng, lowerCutoff, upperCutoff);
    }
}
