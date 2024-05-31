#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using RT = CGAL::Regular_triangulation_3<K>;

using Vector = K::Vector_3;
using Point = RT::Bare_point;
using Weighted_point = RT::Weighted_point;
using Segment = RT::Segment;
using Weight = double;

