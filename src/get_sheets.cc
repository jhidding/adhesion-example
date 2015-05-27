#include "polygon.hh"
#include "cgal_base.hh"

DecoratedMesh<Point, double> get_sheets(
        Rt const &T,
        double l_th)
{
    using PolygonData = Decorated<std::vector<unsigned>, double>;

    std::map<Rt::Cell_handle,unsigned> cell_index;
    auto dual_vertices = make_vector_ptr<Point>();
    auto polygons      = make_vector_ptr<PolygonData>();

    auto stash = [&T, &cell_index, dual_vertices] (
            Rt::Cell_handle const &h) -> unsigned
    {
        if (cell_index.count(h) == 0)
        {
            cell_index[h] = dual_vertices->size();
            dual_vertices->push_back(T.dual(h));
        }

        return cell_index[h];
    };

    for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e)
    {
        std::vector<unsigned> P;

        double l = T.segment(*e).squared_length();
        if (l < l_th) continue;

        auto first = T.incident_cells(*e), c = first;
        bool ok = true;
        do {
            if (T.is_infinite(++c))
            {
                ok = false;
                break;
            }

            P.push_back(stash(c));
        } while (c != first);

        if (ok) polygons->push_back(PolygonData(P, sqrt(l)));
    }

    return DecoratedMesh<Point, double>(dual_vertices, polygons);
}
