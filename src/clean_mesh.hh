#pragma once

#include "polygon.hh"
#include <map>

template <typename Point>
PolygonMesh<Point> clean_mesh(
        PolygonMesh<Point> const &M)
{
    vector_ptr<Point> source_pts;
    vector_ptr<std::vector<unsigned>> source_vertices;
    std::tie(source_pts, source_vertices) = M;

    std::map<unsigned,unsigned> vertex_map;
    auto clean_pts = make_vector_ptr<Point>();
    auto clean_vertices = make_vector_ptr<std::vector<unsigned>>();

    for (auto const &v : *source_vertices)
    {
        std::vector<unsigned> p;

        for (unsigned i : v)
        {
            if (vertex_map.count(i) == 0)
            {
                vertex_map[i] = clean_pts->size();
                clean_pts->push_back((*source_pts)[i]);
            }

            p.push_back(vertex_map[i]);
        }

        clean_vertices->push_back(p);
    }

    return PolygonMesh<Point>(clean_pts, clean_vertices);
}

template <typename Point, typename Info>
DecoratedMesh<Point, Info> clean_mesh(
        DecoratedMesh<Point, Info> const &M)
{
    using PolygonData = Decorated<std::vector<unsigned>, Info>;

    auto source_pts   = std::get<0>(M);
    auto source_verts = std::get<1>(M);

    std::map<unsigned,unsigned> vertex_map;
    auto clean_pts   = make_vector_ptr<Point>();
    auto clean_verts = make_vector_ptr<PolygonData>();

    for (auto const &v : *source_verts)
    {
        PolygonData p(v.get_info());

        for (unsigned i : v)
        {
            if (vertex_map.count(i) == 0)
            {
                vertex_map[i] = clean_pts->size();
                clean_pts->push_back((*source_pts)[i]);
            }

            p.push_back(vertex_map[i]);
        }

        clean_verts->push_back(p);
    }

    return DecoratedMesh<Point, Info>(clean_pts, clean_verts);
}
