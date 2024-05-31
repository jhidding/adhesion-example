/* adhesion_example.cc

    adapted from CGAL example: Triangulation_3/regular_3.cpp

    to compile:
        g++ -std=c++11 -O3 -frounding-math -lm -lrt -lCGAL \
            -lgmp -lboost_thread -lmpfr -lfftw3            \
            adhesion_example.cc lloyd.cc -o ae
 */

#include "cgal_base.hh"

#include <algorithm>
#include <numeric>
#include <iterator>

#include <iostream>
#include <fstream>

#include "boxparam.hh"
#include "grf.hh"
#include "polygon.hh"
#include "write_obj.hh"
#include "write_ply.hh"
#include "sphere.hh"
#include "select_mesh.hh"

// takes the potential function and generates the weighted points
// belonging to this potential. Returns a pointer to a vector of
// Weighted_point.
extern vector_ptr<Weighted_point> make_points(
        BoxParam const &box,
        vector_ptr<double> phi,
        double D);

extern DecoratedMesh<Point, double> get_sheets(
        RT const &T,
        double l_th);

vector_ptr<double> from_file(BoxParam const &box, std::string const &fn)
{
    auto v = make_vector_ptr<double>(box.size);
    std::ifstream fi(fn);
    fi.read(reinterpret_cast<char *>(v->data()), box.size * sizeof(double));
    fi.close();
    return v;
}

template <typename Sel>
void write_selection(
        std::string const &fn,
        DecoratedMesh<Point, double> const &mesh,
        Sel const &selector)
{
    std::cerr << "\n\t" << fn << "\t\t";
    std::ofstream fo(fn);
    write_obj(fo, select_mesh(mesh, selector));
    fo.close();
}

vector_ptr<double> flat_potential(BoxParam const &box)
{
    return make_vector_ptr<double>(box.size, 1.0);
}

/*
void save_dual_example(RT const &T, Point const probe)
{
    using PolygonData = std::vector<unsigned>;

    std::map<RT::Cell_handle,unsigned> cell_index;
    std::map<RT::Vertex_handle,unsigned> vertex_index;
    auto vertices      = make_vector_ptr<Point>();
    auto polygons      = make_vector_ptr<PolygonData>();

    auto stash_dual = [&T, &cell_index, vertices] (
            RT::Cell_handle const &h) -> unsigned
    {
        if (cell_index.count(h) == 0)
        {
            cell_index[h] = vertices->size();
            vertices->emplace_back(T.dual(h));
        }

        return cell_index[h];
    };

    auto stash_delaunay = [&T, &vertex_index, vertices] (
            RT::Vertex_handle const &h) -> unsigned
    {
        if (vertex_index.count(h) == 0)
        {
            vertex_index[h] = vertices->size();
            vertices->emplace_back(T.point(h));
        }

        return vertex_index[h];
    };

    RT::Cell_handle c = T.locate(probe);
    auto v = c->vertex(0);
    std::vector<RT::Edge> edges;
    std::vector<RT::Facet> facets;
    T.incident_edges(c->vertex(0), std::back_inserter(edges));
    T.incident_facets(c->vertex(0), std::back_inserter(facets));

    for (auto &e : edges)
    {
        std::vector<unsigned> P;

        auto first = T.incident_cells(e), c = first;
        bool ok = true;
        do {
            if (T.is_infinite(++c))
            {
                ok = false;
                break;
            }

            P.emplace_back(stash_dual(c));
        } while (c != first);

        if (ok) polygons->emplace_back(P);
    }

    PolygonMesh<Point> M1(vertices, polygons);
    std::ofstream fo1("voronoi.ply");
    write_ply(fo1, M1);
    fo1.close();

    vertices->clear();
    polygons->clear();
    cell_index.clear();

    for (auto &f : facets)
    {
        std::vector<unsigned> P;

        for (unsigned i = 0; i < 4; ++i)
        {
            if (i == f.second) continue;

            P.push_back(stash_delaunay(f.first->vertex(i)));
        }

        polygons->push_back(P);
    }

    PolygonMesh<Point> M2(vertices, polygons);
    std::ofstream fo2("delaunay.ply");
    write_ply(fo2, M2);
    fo2.close();
} */

int main() {
    // parameters of the box
    BoxParam box(8, 10.0);

    // growing mode solution
    double   D    = 1.0;

    // threshold for structure detection
    double   l_th = 0.0; //5.0 * box.res_sqr;

    // the power spectrum of the potential
    auto power_spectrum = [] (Vector const &k) -> double
    {
        double k_sqr = k*k;
        double g = exp(-k_sqr * 0.001 * M_PI * M_PI);
        return 20000*g*(k_sqr == 0 ? 0 : pow(k_sqr, -3.0));
    };

    std::cerr << "Generating potential ... ";

    // Get a potential function, in this case white noise
    // To get prettier results, use a proper GRF
    //auto phi = generate_potential(box, power_spectrum);
    auto phi = flat_potential(box);
    //auto phi = read_potential(box);
    //auto phi = from_file(box, "i2456-128.pot.bin");
    std::cerr << "[done]\n";

    // create a grid with the box specifications and add
    // the weights given by the potential
    std::cerr << "Making points ... ";
        auto pts = make_points(box, phi, D);
    std::cerr << " [done]\n";

    // insert all points in a row.
    std::cerr << "Computing the triangulation ... ";
        RT T;
        T.insert(pts->begin(), pts->end());
    std::cerr << " [done]\n";

    std::cerr << "Writing files ... ";
        auto sheets_mesh = get_sheets(T, l_th);

    // save_dual_example(T, Point(5,5,5));

	Sphere<K> HOME(Point(5,5,5), 5.0);
	write_selection("example-mesh.obj", sheets_mesh, HOME);

    /*  Sphere<K> HOME(Point(90,90,90), 30.0);
        write_selection("home30.obj", sheets_mesh, HOME);

        Sphere<K> A262(Point(134, 69, 89), 30.0);
        write_selection("perseus30.obj", sheets_mesh, A262); */
    std::cerr << "[done]\n";

    std::cerr << "Bye!\n";
    return 0;
}
