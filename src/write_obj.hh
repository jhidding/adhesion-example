#pragma once

#include "polygon.hh"

template <typename Point>
void write_obj(std::ostream &out, DecoratedMesh<Point,double> const &M)
{
    auto points = std::get<0>(M);
    auto polygons = std::get<1>(M);

    for (Point const &p : *points)
    {
        out << "v " << p << " 1.0\n";
    }
    out << "\n";

    std::vector<double> val;
    double min = 1e6, max = 0.0;
    for (auto const &p : *polygons)
    {
        double a = p.get_info();
        if (a < min) min = a;
        if (a > max) max = a;
    }

    unsigned i = 0;
    for (auto const &p : *polygons)
    {
        double a = p.get_info();
        out << "vt " << (a - min)/(max-min) << " 0\n";
        ++i;
    }
    out << "\n";

    i = 1;
    for (auto const &p : *polygons)
    {
        out << "f";
        for (unsigned j : p)
            out << " " << j+1 << "/" << i;
        out << "\n";
        ++i;
    }
}
