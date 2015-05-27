#pragma once

#include "polygon.hh"
#include <algorithm>

template <typename Point, typename Surface>
PolygonPair<Point> split_polygon(
        Polygon<Point> const &P,
        Surface const &S,
        bool closed = true)
{
    vector_ptr<Point> pts;
    std::vector<unsigned> orig, r1, r2;

    std::tie(pts, orig) = P;

    auto is_below = [&pts, &S] (unsigned i) -> bool {
        return (S.oriented_side((*pts)[i]) == -1);
    };

    if (closed)
        orig.push_back(orig.front());

    auto i = orig.begin();
    auto j = i; ++j;
    bool below = is_below(*i);

    while (j != orig.end())
    {
        if (below != is_below(*j)) // surface crossed
        {
            auto q = S.intersect((*pts)[*i], (*pts)[*j]);

            if (q)
            {
                r1.push_back(pts->size());
                r2.push_back(pts->size());
                pts->push_back(*q);
                std::swap(r1, r2);
            }

            below = not below;
        }
        else
        {
            r1.push_back(*j);
            ++i; ++j;
        }
    }

    if (below)
        return PolygonPair<Point>(pts, r1, r2);
    else
        return PolygonPair<Point>(pts, r2, r1);
}
