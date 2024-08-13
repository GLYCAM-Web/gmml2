#ifndef INCLUDES_CODEUTILS_METROPOLISCRITERION_HPP
#define INCLUDES_CODEUTILS_METROPOLISCRITERION_HPP

#include "includes/External_Libraries/PCG/pcg_random.h" // $GEMSHOME/gmml/includes is in the makefile
#include <cmath>
#include <stdio.h>
#include <random>

namespace monte_carlo
{
    inline double get_random_acceptance_probability()
    {
        // I did this when testing and calling multiple times, as time wasn't progressing between calls. So Fast!
        // srand ((time(NULL) + rand())); // initialise a random seed for rand(). Otherwise it's always the same.
        // std::cout << (rand() % 100) << "rand()\n";
        // double r = (rand() % 100); // get a number between 1 and 100
        // return (r / 100); // get it between 0 and 1

        // Seed with a real random value, if available
        pcg_extras::seed_seq_from<std::random_device> metropolis_seed_source;
        // Make a random number engine
        pcg32 rng_engine(metropolis_seed_source);
        std::uniform_real_distribution<> real_number_distribution(0, 1); // define the range
        return real_number_distribution(rng_engine);
    }

    inline bool accept_via_metropolis_criterion(double change_in_overlap)
    {
        if (change_in_overlap < 0)
        {
            // std::cout << "ACCEPTED as better: " << change_in_overlap << "\n";
            return true;
        }
        double r = get_random_acceptance_probability();
        double p = exp(-change_in_overlap);
        if (p > r)
        {
            // std::cout << "ACCEPTED: " << change_in_overlap << " p: " << p << " r: " << r << "\n";
            return true;
        }
        // std::cout << "REJECTED: " << change_in_overlap << " p: " << p << " r: " << r << "\n";
        return false;
    }
} // namespace monte_carlo
#endif
