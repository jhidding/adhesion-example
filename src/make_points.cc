#include "boxparam.hh"
#include "polygon.hh"
#include "interpolate.hh"
#include "lloyd.hh"

#include "cgal_base.hh"

#include <iostream>

vector_ptr<Weighted_point> make_points(
        BoxParam const &box,
        vector_ptr<double> phi,
        double D)
{
    BoxParam glass_tile(box.N/8, box.L);
    std::cerr << "\n\t Crafting glass ... ";
    auto glass = Lloyd::make_glass(512, 10);
    std::cerr << "[done]: " << glass->size() << " points in glass.\n";
    auto points = make_vector_ptr<Weighted_point>();
    auto f = interpolate<Vector>(box, phi);

    std::cerr << "\t Interpolating potential ... ";
    Point origin(0, 0, 0);
    for (size_t i = 0; i < glass_tile.size; ++i)
    {
        Vector tile_origin = glass_tile.point<Vector>(i);
        for (Point const &p : *glass)
        {
            Vector v = tile_origin + (p - origin) * glass_tile.res;
            points->push_back(Weighted_point(origin + v, 2 * D * f(v)));
        }
    }
    std::cerr << "[done] . ";

    return points;
}
