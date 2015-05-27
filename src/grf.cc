#include "grf.hh"

// returns a function that returns gaussian distributed
// noise, mean = 0, sigma = 1.0.
std::function<double ()> gaussian_white_noise(
        unsigned long seed)
{
    std::shared_ptr<std::mt19937>
        random(new std::mt19937(seed));
    std::shared_ptr<std::normal_distribution<double>>
        normal(new std::normal_distribution<double>);

    return [random, normal] () -> double
        { return (*normal)(*random); };
}
