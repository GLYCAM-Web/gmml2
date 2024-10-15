#ifndef INCLUDES_CODEUTILS_RANDOM_HPP
#define INCLUDES_CODEUTILS_RANDOM_HPP

#include "includes/External_Libraries/PCG/pcg_random.h"

#include <vector>
#include <algorithm>
#include <random>

namespace codeUtils
{
    inline uint64_t generateRandomSeed()
    {
        std::random_device rdev;
        return (uint64_t(rdev()) << 32) | rdev();
    }

    template<class T> int uniformRandomVectorIndex(pcg32& rng, const std::vector<T>& vec)
    {
        std::uniform_int_distribution<> distr(0, (vec.size() - 1));
        return distr(rng);
    }

    template<class T> T uniformRandomVectorEntry(pcg32& rng, const std::vector<T>& vec)
    {
        return vec.at(uniformRandomVectorIndex(rng, vec));
    }

    inline double uniformRandomDoubleWithinRange(pcg32& rng, double lower, double upper)
    {
        std::uniform_real_distribution<> distr(lower, upper);
        return distr(rng);
    }

    inline double normalDistributionRandomDouble(pcg32& rng)
    {
        std::normal_distribution<> distr(0.0, 1.0);
        return distr(rng);
    }

    template<class T> std::vector<T> shuffleVector(pcg32& rng, std::vector<T> vec)
    {
        std::shuffle(vec.begin(), vec.end(), rng);
        return vec;
    }

    template<class T> size_t randomIndex(pcg32& rng, std::vector<T>& vec)
    {
        if (vec.empty())
        {
            throw std::runtime_error("can't generate random index on empty vector");
        }
        std::uniform_int_distribution<> dist(0, vec.size() - 1);
        return dist(rng);
    }

    template<class T> std::vector<T> sample(pcg32& rng, size_t count, std::vector<T>& vec)
    {
        std::vector<T> result;
        result.reserve(count);
        std::sample(vec.begin(), vec.end(), std::back_inserter(result), count, rng);
        return result;
    }

    std::vector<size_t> weightedRandomOrder(pcg32& rng, std::vector<double> weights);
    double normalDistributionRandomDoubleWithCutoff(pcg32& rng, double lowerCutoff, double upperCutoff);
} // namespace codeUtils
#endif
