#include "includes/CodeUtils/random.hpp"

#include "includes/CodeUtils/containers.hpp"

#include <vector>
#include <utility>
#include <numeric>
#include <random>

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
