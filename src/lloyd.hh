#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

#include "polygon.hh"

namespace Lloyd
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Periodic_3_triangulation_traits_3<K> GT;
    typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT> PDT;
    typedef PDT::Cell_handle    Cell_handle;
    typedef PDT::Vertex_handle  Vertex_handle;
    typedef PDT::Locate_type    Locate_type;
    typedef PDT::Point          Point;
    typedef PDT::Iso_cuboid     Iso_cuboid;

    extern std::function<double ()> uniform_noise(unsigned long seed);
    extern vector_ptr<Point> make_glass(size_t N, unsigned m);
}
