// #include <pybind11/pybind11.h> // should be added once we start working on the python part
#include "clipper2/clipper.h"
#include "rapidjson/document.h"
#include "homog2d/homog2d.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <thread>
#include <future>

#include <chrono> // temporay to allow performance measurements

#include "const.h"
#include "structs.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace std;
using namespace rapidjson;
using namespace Clipper2Lib;

PathsD h2dpoints_to_pathsd(const vector<h2d::Point2d> &points)
{
    string s;
    for (h2d::Point2d p : points)
    {
        s += to_string(p.getX()) + " ";
        s += to_string(p.getY()) + " ";
    }
    PathsD p;
    p.push_back(MakePathD(s));
    return p;
}

PathsD soilpolygon_to_pathsd(SoilPolygon &soilpolygon)
{
    string s;
    PathsD p = h2dpoints_to_pathsd(soilpolygon.points);
    return p;
}

// double arr[] = {x1, top, x2, top, x2, z2, x3, z3, x1, z1};
//         string s_cut_slice;
//         for (double d : arr)
//         {
//             s_cut_slice += to_string(d) + " ";
//         }

//         PathsD cut_slice;
//         cut_slice.push_back(MakePathD(s_cut_slice));

/*
This function will calculate the Bishop safety factor for the given
model, centerpoint of the slope circle and z coordinate of the tangent line
*/
void sf_bishop(const int i, const BishopModel &model, double mx, double mz, double z_tangent, double *sf)
{
    // get the highest point on the surface
    double top = 0.0;
    for (auto pt : model.surface.getPts())
    {
        if (pt.getY() > top)
        {
            top = pt.getY();
        }
    }

    // create the circle
    double radius = mz - z_tangent;
    h2d::Circle circle(mx, mz, radius);

    // get the intersections between the circle and the polyline
    auto res = circle.intersects(model.surface);
    // only proceed if we have at least 2 intersections
    if (!res())
    {
        *sf = -9999;
        return;
    }
    if (res.size() < 2)
    {
        *sf = -9999;
        return;
    }

    auto pts = res.get();
    double entry_point = pts[0].getX();
    double exit_point = pts[1].getX();

    cout << "x coordinate entry point: " << entry_point << endl;
    cout << "x coordinate exit point : " << exit_point << endl;

    // # divide the space between the entry and exit point into evenly distributed points
    vector<double> slice_coordinates = {};
    double dl = exit_point - entry_point;
    for (int i = 0; i < NUM_S; ++i)
    {
        double x = entry_point + dl * double(i) / (NUM_S - 1);
        slice_coordinates.push_back(x);
    }

    // slices = []
    for (int i = 1; i < slice_coordinates.size(); i++)
    {
        double x1 = slice_coordinates[i - 1];                        // left side of the slice
        double x2 = slice_coordinates[i];                            // right side of the slice
        double x3 = (x1 + x2) / 2.0;                                 // mid point of the slice
        double z1 = mz - sqrt(pow(radius, 2.0) - pow(mx - x1, 2.0)); // z value at bottom left of the slice
        double z2 = mz - sqrt(pow(radius, 2.0) - pow(mx - x2, 2.0)); // z value at bottom right of the slice
        double z3 = mz - sqrt(pow(radius, 2.0) - pow(mx - x3, 2.0)); // z value at bottom mid of the slice

        // ready to create the slice polygon (clockwise)
        double arr[] = {x1, top, x2, top, x2, z2, x3, z3, x1, z1};
        string s_cut_slice;
        for (double d : arr)
        {
            s_cut_slice += to_string(d) + " ";
        }

        PathsD cut_slice;
        cut_slice.push_back(MakePathD(s_cut_slice));

        // now get all the intersections with the soilpolygons above the phreatic line
        for (SoilPolygon spg : model.soilpolygons_above_pl)
        {
            PathsD p = soilpolygon_to_pathsd(spg);
        }

        // cout << "x slice left            : " << x1 << endl;
        // cout << "x slice middle          : " << x3 << endl;
        // cout << "x slice right           : " << x2 << endl;
        // cout << "bottom slice left       : " << z1 << endl;
        // cout << "bottom slice middle     : " << z3 << endl;
        // cout << "bottom slice right      : " << z2 << endl;
    }

    int j = 1;

    //     surface_points = (
    //         self.calculation_model.soilprofile2.get_surface_within_limits(x1, x2) # als je de surface hebt moet dit makkelijk in cpp kunnen
    //     )

    //     polygon = surface_points + [(x2, z2), (x3, z3), (x1, z1)]

    //     (
    //         soilpolygons_above_phreatic_layer,
    //         soilpolygons_below_phreatic_layer,
    //     ) = self.calculation_model.soilprofile2.get_slice( # dit moet met clipper
    //         polygon, self.calculation_model.phreatic_line
    //     )

    //     base_alpha = math.atan2((mz - z3), (mx - x3)) - 0.5 * math.pi
    //     b = x2 - x1  # width of the slice
    //     base_L = b / math.cos(base_alpha)  # length at bottom of slice
    //     u = self.calculation_model.u_at(x3, z3) # dit kun je makkelijk berekenen met de phreatic line

    //     # calculate weight
    //     W = 0
    //     for spg in soilpolygons_above_phreatic_layer:
    //         soil = self.calculation_model.soilcollection.get(spg.soilcode)
    //         W += spg.area * soil.y_dry
    //     for spg in soilpolygons_below_phreatic_layer:
    //         soil = self.calculation_model.soilcollection.get(spg.soilcode)
    //         W += spg.area * soil.y_sat

    //     soilcode = self.calculation_model.soilprofile2.soilcode_at(x3, z3) # dit moet ook makkelijk kunnen
    //     if soilcode == "":
    //         raise ValueError(
    //             f"No soil found at the given coordinates.. fix me please :-/"
    //         )
    //     soil = self.calculation_model.soilcollection.get(soilcode)
    //     c = soil.cohesion
    //     phi = soil.friction_angle

    //     slices.append(
    //         (i, b, W, -1.0 * base_alpha, base_L, u, c, (phi / 180.0) * math.pi)
    //     )

    // M = np.array(slices)
    // # M[:,0] = index
    // # M[:,1] = b (width of slice)
    // # M[:,2] = W (weight of slice)
    // # M[:,3] = alpha
    // # M[:,4] = L
    // # M[:,5] = u
    // # M[:,6] = c
    // # M[:,7] = phi (in radians)
    // denom = np.sum(M[:, 2] * np.sin(M[:, 3]))
    // cl = np.sum(M[:, 6] * M[:, 4])

    // # assume first sf
    // sf = 1.0

    // iteration = 0
    // while 1:
    //     N1 = M[:, 6] * M[:, 4] * np.sin(M[:, 3])
    //     N2 = M[:, 5] * M[:, 4] * np.sin(M[:, 3] * np.tan(M[:, 7]))
    //     N3 = np.cos(M[:, 3]) + (np.sin(M[:, 3]) * np.tan(M[:, 7])) / sf
    //     N = (M[:, 2] - (N1 - N2) / sf) / N3

    //     fos = (cl + np.sum((N - M[:, 5] * M[:, 4]) * np.tan(M[:, 7]))) / denom

    //     if abs(sf - fos) < 0.01:
    //         break

    //     sf = (sf + fos) / 2.0

    //     iteration += 1
    //     if iteration > 20:
    //         return None
    //         # raise ValueError(f"Calculation not converging...")

    // return sf
}

/*
This function will parse a json string from the typical leveelogic calculation model
and create a BishopModel structure to pass to the actual calculation
*/
BishopModel parse_bishop_model(const string &json)
{
    // create a document using rapidjson
    Document document;
    document.Parse(json.c_str());
    assert(document.IsObject());

    // PARSE SOILCOLLECTION
    vector<Soil> soils = {};
    Value &v_soilcollection = document["soilcollection"]["soils"];
    assert(v_soilcollection.IsArray());
    for (SizeType i = 0; i < v_soilcollection.Size(); i++)
    {
        soils.push_back(Soil{
            v_soilcollection[i]["code"].GetString(),
            v_soilcollection[i]["color"].GetString(),
            v_soilcollection[i]["y_dry"].GetDouble(),
            v_soilcollection[i]["y_sat"].GetDouble(),
            v_soilcollection[i]["cohesion"].GetDouble(),
            v_soilcollection[i]["friction_angle"].GetDouble(),
        });
    }

    // PARSE SOILPROFILE 2
    // keep track of the limits
    double xmin = 1e9;
    double xmax = -1e9;
    double zmin = 1e9;
    double zmax = -1e9;
    vector<SoilPolygon> soilpolygons = {};
    Value &v_soilprofile2 = document["soilprofile2"]["soilpolygons"];
    assert(v_soilprofile2.IsArray());

    for (SizeType i = 0; i < v_soilprofile2.Size(); i++)
    {
        Value &v_points = v_soilprofile2[i]["points"];
        vector<h2d::Point2d> points;
        for (SizeType i = 0; i < v_points.Size(); i++)
        {
            double x = v_points[i][0].GetDouble();
            if (x < xmin)
                xmin = x;
            if (x > xmax)
                xmax = x;
            double z = v_points[i][1].GetDouble();
            if (z < zmin)
                zmin = z;
            if (z > zmax)
                zmax = z;
            points.push_back(h2d::Point2d{x, z});
        }
        string soilcode = v_soilprofile2[i]["soilcode"].GetString();

        soilpolygons.push_back(SoilPolygon{points, soilcode});
    }

    // PARSE PHREATIC LINE
    vector<h2d::Point2d> plpoints = {};
    Value &v_phreatic_line = document["phreatic_line"];
    assert(v_phreatic_line.IsArray());
    for (SizeType i = 0; i < v_phreatic_line.Size(); i++)
    {
        double x = v_phreatic_line[i][0].GetDouble();
        double z = v_phreatic_line[i][1].GetDouble();
        plpoints.push_back(h2d::Point2d{x, z});
    }
    h2d::OPolyline phreatic_line = h2d::OPolyline{plpoints};

    // PARSE SURFACE LINE
    vector<h2d::Point2d> spoints = {};
    Value &v_surface_line = document["surface_line"];
    assert(v_surface_line.IsArray());
    for (SizeType i = 0; i < v_surface_line.Size(); i++)
    {
        double x = v_surface_line[i][0].GetDouble();
        double z = v_surface_line[i][1].GetDouble();
        spoints.push_back(h2d::Point2d{x, z});
    }
    h2d::OPolyline surface = h2d::OPolyline{spoints};

    // PARSE BISHOP SEARCH GRID
    Value &v_bishop_search_grid = document["bishop_search_grid"];
    double left = v_bishop_search_grid["x_left"].GetDouble();
    double bottom = v_bishop_search_grid["z_bottom"].GetDouble();
    double width = v_bishop_search_grid["width"].GetDouble();
    double height = v_bishop_search_grid["height"].GetDouble();
    double tangents_top = v_bishop_search_grid["tangents_top"].GetDouble();
    double tangents_bottom = v_bishop_search_grid["tangents_bottom"].GetDouble();
    double minimum_slip_plane_length = v_bishop_search_grid["minimum_slip_plane_length"].GetDouble();

    BishopSearchGrid bishop_search_grid = BishopSearchGrid{
        left,
        bottom,
        width,
        height,
        tangents_top,
        tangents_bottom,
        minimum_slip_plane_length,
    };

    // Extra step, split soilpolygons in above and below phreatic line
    // create the placeholders for the result
    vector<SoilPolygon> soilpolygons_above_pl = {};
    vector<SoilPolygon> soilpolygons_below_pl = {};

    // split the soilpolygons in those above and below the phreatic line, saves calculation time
    if (phreatic_line.size() > 0) // do we have a phreatic line?
    {
        // polygon above the phreatic line
        vector<h2d::Point2d> points_above_pl = {};
        points_above_pl.push_back(h2d::Point2d(xmin, zmax + 1.0));
        points_above_pl.push_back(h2d::Point2d(xmax, zmax + 1.0));
        for (int i = phreatic_line.size() - 1; i >= 0; i--)
        {
            points_above_pl.push_back(phreatic_line.getPoint(i));
        }
        PathsD polygon_above_pl = h2dpoints_to_pathsd(points_above_pl);

        // cut out the soillayers above the phreatic line
        for (SoilPolygon spg : soilpolygons)
        {
            PathsD polygon = soilpolygon_to_pathsd(spg);
            PathsD solution = Intersect(polygon_above_pl, polygon, FillRule::NonZero);

            if (solution.size() > 0)
            {
                for (PathD path : solution)
                {
                    SoilPolygon new_soilpolygon = {};
                    new_soilpolygon.soilcode = spg.soilcode;
                    for (PointD p : path)
                    {
                        new_soilpolygon.points.push_back(p);
                    }
                    soilpolygons_above_pl.push_back(new_soilpolygon);
                }
            }
            int j = 1;
        }

        // polygon below the phreatic line
        vector<h2d::Point2d> points_below_pl = {};
        for (h2d::Point2d p : phreatic_line.getPts())
        {
            points_below_pl.push_back(p);
        }
        points_below_pl.push_back(h2d::Point2d(xmax, zmin - 1.0));
        points_below_pl.push_back(h2d::Point2d(xmin, zmin - 1.0));
        PathsD polygon_below_pl = h2dpoints_to_pathsd(points_below_pl);

        // cut out the soillayers below the phreatic line
        for (SoilPolygon spg : soilpolygons)
        {
            PathsD polygon = soilpolygon_to_pathsd(spg);
            PathsD solution = Intersect(polygon_below_pl, polygon, FillRule::NonZero);

            if (solution.size() > 0)
            {
                for (PathD path : solution)
                {
                    SoilPolygon new_soilpolygon = {};
                    new_soilpolygon.soilcode = spg.soilcode;
                    for (PointD p : path)
                    {
                        new_soilpolygon.points.push_back(p);
                    }
                    soilpolygons_below_pl.push_back(new_soilpolygon);
                }
            }
            int j = 1;
        }
    }
    else
    { // there is no phreatic line so everything is above the phreatic line
        for (SoilPolygon spg : soilpolygons)
        {
            soilpolygons_above_pl.push_back(spg);
        }
    }

    return BishopModel{
        soils,
        soilpolygons_above_pl,
        soilpolygons_below_pl,
        bishop_search_grid,
        phreatic_line,
        surface,
    };
}

vector<double> calculate_bishop() // will become calculate_bishop(const string &json)
{
    // for now we skip the given string and read a test file
    ifstream file("../test/model.json");
    stringstream buffer{};
    buffer << file.rdbuf();
    string json = buffer.str();

    // get the model from the string
    BishopModel model = parse_bishop_model(json);

    // DEBUG
    model.print();
    // END DEBUG

    double x = model.bishop_search_grid.left;
    double z = model.bishop_search_grid.bottom;
    double dx = model.bishop_search_grid.width / double(NUM_X - 1);
    double dz = model.bishop_search_grid.height / double(NUM_Z - 1);
    double dt = (model.bishop_search_grid.tangents_top - model.bishop_search_grid.tangents_bottom) / double(NUM_T - 1);

    // temporary code to measure performance
    auto start = std::chrono::high_resolution_clock::now();

    // iterate over the possible slope circle locations, multithreaded
    array<thread, NUM_T * NUM_X * NUM_Z> threads;
    array<double, NUM_T * NUM_X * NUM_Z> sfs;
    int i = 0;
    for (int nx = 0; nx < NUM_X; ++nx)
    {
        for (int nz = 0; nz < NUM_Z; ++nz)
        {
            for (int nt = 0; nt < NUM_T; ++nt)
            {
                double x = model.bishop_search_grid.left + nx * dx;
                double z = model.bishop_search_grid.bottom + nz * dz;
                double t = model.bishop_search_grid.tangents_bottom + nt * dt;

                threads[i] = thread(sf_bishop, i, model, x, z, t, &sfs[i]);
                threads[i].join(); // remove this to enable threading again

                ++i;

                nt = 100; // remove this to allow all situations to be calculated
                nz = 100;
                nx = 100;
            }
        }
    }

    // for (auto &t : threads) // uncomment this to add threading again
    // {
    //     t.join();
    // }

    // temporary code to measure performance
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Elapsed time " << elapsed.count() << " ms\n";

    // for (auto &sf : sfs)
    // {
    //     std::cout << "sf =" << sf << endl;
    // }

    return {};
}

int main()
{
    vector<double> sfs = calculate_bishop(); // will become calculate_bishop("jsonstring");
}

// THE NEXT CODE IS COMMENTED SINCE I USE THIS AS A CPP PROJECT FIRST
// ONCE THIS IS DONE IT WILL BE PART OF THE PYTHON LIBRARY CODE
// namespace py = pybind11;

// PYBIND11_MODULE(pstab, m)
// {
//     m.doc() = R"pbdoc(
//         Pybind11 example plugin
//         -----------------------

//         .. currentmodule:: pstab

//         .. autosummary::
//            :toctree: _generate

//            bishop
//     )pbdoc";

//     m.def("bishop", &bishop, R"pbdoc(
//         Calculate the Bishop safety factor from the given model

//         Some other explanation about the bishop function.
//     )pbdoc");

// #ifdef VERSION_INFO
//     m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
// #else
//     m.attr("__version__") = "dev";
// #endif
// }
