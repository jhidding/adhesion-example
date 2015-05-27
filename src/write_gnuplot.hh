#pragma once

#include <iostream>

// writes nodes and filaments to a text file suitable for Gnuplot
template <typename Rt>
void write_gnuplot(
        Rt const &T,
        double l_th,
        std::ostream &out)
{
    // write nodes
    for (auto h = T.finite_cells_begin(); h != T.finite_cells_end(); ++h)
    {
        int cnt = 0;
        for (unsigned i = 1; i < 4; ++i)
            for (unsigned j = 0; j < i; ++j)
                if (T.segment(h, i, j).squared_length() > l_th)
                    ++cnt;

        if (cnt == 6) // all edges of the cell exceed threshold
        {
            // output the Eulerian location of the node and its mass.
            out << T.dual(h) << " "
                << T.tetrahedron(h).volume() << std::endl;
        }
    }

    out << "\n\n";
    // write filaments
    for (auto f = T.finite_facets_begin(); f != T.finite_facets_end(); ++f)
    {
        double A = T.triangle(*f).squared_area();

        auto c1 = f->first;  // cell-handle to which facet f belongs
        auto c2 = T.mirror_facet(*f).first;

        if (T.is_infinite(c1) or T.is_infinite(c2))
            continue;

        int i = f->second;  // index of vertex oposite facet
        int cnt = 0;
        for (unsigned j = 1; j < 4; ++j)
            if (j != i) for (unsigned k = 0; k < j; ++k)
                if (k != i && T.segment(c1, j, k).squared_length() > l_th)
                    ++cnt;

        if (cnt == 3)
        {
            out << T.dual(c1) << " "
                << T.dual(c2) - T.dual(c1) << " " << A << std::endl;
        }
    }
}
