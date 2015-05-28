#pragma once

#include "polygon.hh"

template <typename K>
class Plane
{
    using Point   = typename K::Point_3;
    using Vector  = typename K::Vector_3;

    Point  origin;
    Vector normal; // pointing to oriented positive side

public:
    Sphere(Point const &p, double r):
        origin(p), radius_squared(r*r) {}

    int oriented_side(Point const &p) const
    {
        double d = normal * (p - origin);

        if (d < 0)
            return -1;

        if (d > 0)
            return +1;

        return 0;
    }

    std::shared_ptr<Point> intersect(Point const &a, Point const &b) const
    {
        if (oriented_side(a) * oriented_side(b) >= 0)
            return std::shared_ptr<Point>();

        Vector v = b - a;
        double t = normal * (origin - a) / (normal * v);
        return std::make_shared<Point>(a + t*v);
    }
};
