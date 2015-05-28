#include "lloyd.hh"

using namespace Lloyd;

std::function<double ()> Lloyd::uniform_noise(
        unsigned long seed)
{
    auto random = std::make_shared<std::mt19937>(seed);
    auto dist   = std::make_shared<std::uniform_real_distribution<double>>();

    return [random, dist] () -> double
        { return (*dist)(*random); };
}

vector_ptr<Point> Lloyd::make_glass(size_t N, unsigned m)
{
    Iso_cuboid domain(0, 0, 0, 1.0, 1.0, 1.0);
    auto points = make_vector_ptr<Point>(N);
    auto noise  = uniform_noise(0);
    for (Point &p : *points)
        p = Point(noise(), noise(), noise());

    for (unsigned i = 0; i < m; ++i)
    {
        PDT T(points->begin(), points->end(), domain);
        auto centroids = make_vector_ptr<Point>();

        if (not T.is_triangulation_in_1_sheet())
            throw "need a one-sheet triangulation.";

        for (auto p  = T.vertices_begin();
                  p != T.vertices_end();
                  ++p)
            centroids->push_back(T.dual_centroid(p));

        points = centroids;
    }

    return points;
}
