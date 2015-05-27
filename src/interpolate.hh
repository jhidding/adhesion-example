#pragma once
#include "polygon.hh"
#include "boxparam.hh"

template <typename Vector>
std::function<double (Vector const &)> interpolate(
    BoxParam const &box,
    vector_ptr<double> f)
{
    return [box, f] (Vector const &x_) -> double
    {
        static int const block[8][3] = {
            {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
            {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}
        };

        Vector  x = x_ / box.res;
        int    origin[3];
        double A[2][3];

        for (unsigned k = 0; k < 3; ++k)
        {
            origin[k] = (int)(x[k]);
            A[1][k] = x[k] - origin[k];
            A[0][k] = 1 - A[1][k];
        }

        double v = 0.0;

        for (unsigned i = 0; i < 8; ++i)
        {
            double z = 1;
            int p[3];

            for (unsigned k = 0; k < 3; ++k)
            {
                z *= A[block[i][k]][k];
                p[k] = origin[k] + block[i][k];
            }

            v += (*f)[box.idx(p)] * z;
        }

        return v;
    };
}
