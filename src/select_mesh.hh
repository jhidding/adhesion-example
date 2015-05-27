#pragma once

#include "split_polygon.hh"
#include "clean_mesh.hh"

template <typename Point, typename Surface>
PolygonMesh<Point> select_mesh(
        PolygonMesh<Point> const &M,
        Surface const &S)
{
    vector_ptr<Point> pts;
    vector_ptr<std::vector<unsigned>> source_vertices;
    std::tie(pts, source_vertices) = M;

    auto target_vertices = make_vector_ptr<std::vector<unsigned>>();

    for (std::vector<unsigned> const &v : *source_vertices)
    {
        auto split = split_polygon(Polygon<Point>(pts, v), S);
        auto below = std::get<1>(split);
        if (below.size() > 0)
            target_vertices->push_back(below);
    }

    return clean_mesh(PolygonMesh<Point>(pts, target_vertices));
}

template <typename Point, typename Surface, typename Info>
DecoratedMesh<Point,Info> select_mesh(
        DecoratedMesh<Point, Info> const &M,
        Surface const &S)
{
    using PolygonData = Decorated<std::vector<unsigned>, Info>;

    auto pts = std::get<0>(M);
    auto source = std::get<1>(M);
    auto target = make_vector_ptr<PolygonData>();

    for (PolygonData const &v : *source)
    {
        Info info = v.get_info();
        auto split = split_polygon(Polygon<Point>(pts, v), S);
        auto below = std::get<1>(split);
        if (below.size() > 0)
            target->push_back(PolygonData(below, info));
    }

    return clean_mesh(DecoratedMesh<Point,Info>(pts, target));
}
