#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Vector_3                                    Vector;
typedef Traits::Weighted_point                              Weighted_point;
typedef CGAL::Regular_triangulation_3<Traits>               Rt;
