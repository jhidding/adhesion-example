#pragma once

#include "polygon.hh"

template <typename Point>
void write_ply(std::ostream &out, PolygonMesh<Point> const &M)
{
    vector_ptr<Point> points;
    vector_ptr<std::vector<unsigned>> polygons;
    std::tie(points, polygons) = M;

    out << "ply\n"
        << "format ascii 1.0\n"
        << "element vertex " << points->size() << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n"
        << "element face " << polygons->size() << "\n"
        << "property list uchar int vertex_index\n"
        << "end_header\n";

    for (Point const &p : *points)
    {
        out << p << "\n";
    }

    for (auto const &p : *polygons)
    {
        out << p.size();
        for (unsigned i : p)
            out << " " << i;
        out << "\n";
    }
}
