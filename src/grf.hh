#pragma once
#include "fft.hh"
#include "boxparam.hh"
#include "polygon.hh"
#include <functional>
#include <random>
#include "cgal_base.hh"

extern std::function<double ()> gaussian_white_noise(
        unsigned long seed);

template <typename Function>
vector_ptr<double> generate_potential(
        BoxParam const &box,
        Function power_spectrum)
{
    FFT3 fft(box);

    auto noise = gaussian_white_noise(23);
    for (std::complex<double> &z : fft.real_space)
        z = noise();

    fft.forward_transform();
    for (size_t i = 0; i < box.size; ++i)
        fft.fourier_space[i] *= sqrt(power_spectrum(box.k_space<Vector>(i)));

    fft.backward_transform();
    auto phi = make_vector_ptr<double>(box.size);
    for (size_t i = 0; i < box.size; ++i)
        (*phi)[i] = fft.real_space[i].real();

    return phi;
}
