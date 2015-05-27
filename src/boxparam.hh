#pragma once
#include <cmath>
#include <cstdlib>

// a struct to collects the paramaters of the box that most
// functions will need.
struct BoxParam {
    unsigned N;         // number of grid points in a row
    size_t   size;      // number of points in the entire box
    double   L;         // physical size of the box
    double   res;       // resolution of the box
    double   res_sqr;   // square of the resolution

    BoxParam(unsigned N_, double L_):
        N(N_),
        size(N*N*N),
        L(L_),
        res(L/N),
        res_sqr(res*res) {}

    // BoxParam::point takes an index into an array and
    // returns the grid point belonging to that index
    template <typename Point>
    Point point(size_t i) const
    {
        int x = i % N;
        int y = (i / N) % N;
        int z = i / (N*N);

        return Point(x * res, y * res, z * res);
    }

    size_t idx(int i[3]) const
    {
        return (i[0]%N) + (i[1]%N) * N + (i[2]%N) * N*N;
    }

    double wave_number(int i) const
    {
        return (int(i) > int(N)/2 ? int(i) - int(N) : int(i)) * (2*M_PI)/L;
    }

    template <typename Point>
    Point k_space(size_t i) const
    {
        int x = i % N;
        int y = (i / N) % N;
        int z = i / (N*N);

        return Point(
            wave_number(x),
            wave_number(y),
            wave_number(z));
    }
};
