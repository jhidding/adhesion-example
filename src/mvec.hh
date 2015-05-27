#pragma once
#include <array>

namespace J
{
    template <typename T, unsigned R>
    class mVec: public std::array<T, R>
    {
    public:
        using std::array<T,R>::array;
    };

    class Point: public mVec<double, 3>
    {
    public:
        Point(double x, double y, double z):
            mVec<double, 3>({x, y, z}) {}

        double x() const { return (*this)[0]; }
        double y() const { return (*this)[0]; }
        double z() const { return (*this)[0]; }
    };
}
