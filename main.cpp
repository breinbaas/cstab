#include "clipper2/clipper.h"
#include "rapidjson/document.h"
#include "homog2d/homog2d.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <pthread.h>
#include <future>
#include <cmath>

#include <chrono> // temporary to allow performance measurements

#include "const.h"
#include "structs.h"

using namespace std;
using namespace rapidjson;
using namespace Clipper2Lib;

/*
Function to create a clipper path from homog2d points
*/
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

/*
Function to convert a soilpolygon to a clipper path
*/
PathsD soilpolygon_to_pathsd(SoilPolygon &soilpolygon)
{
    string s;
    PathsD p = h2dpoints_to_pathsd(soilpolygon.points);
    return p;
}

/*
Function to calculate the water pressure at the given x, z location
*/
double u_at(double x, double z, const h2d::OPolyline &phreatic_line)
{
    for (int i = 1; i < phreatic_line.size(); i++)
    {
        if (phreatic_line.getPoint(i - 1).getX() <= x && x <= phreatic_line.getPoint(i).getX())
        {
            double x1 = phreatic_line.getPoint(i - 1).getX();
            double x2 = phreatic_line.getPoint(i).getX();
            double z1 = phreatic_line.getPoint(i - 1).getY();
            double z2 = phreatic_line.getPoint(i).getY();
            double zt = z1 + (x - x1) / (x2 - x1) * (z2 - z1);
            return (zt - z) * VOL_WEIGHT_WATER;
        }
    }
    return 0.0;
}

/*
Function to get the soil from a given soilcode
*/
Soil get_soil(const string &soilcode, const vector<Soil> &soils)
{
    for (Soil soil : soils)
    {
        if (soil.code.compare(soilcode) == 0)
        {
            return soil;
        }
    }
    throw "received unknown soilcode";
}

/*
Function to get a soilcode from a given x,z point
*/
string get_soil_at(const double x, const double z, const vector<SoilPolygon> &soilpoylgons)
{
    for (SoilPolygon spg : soilpoylgons)
    {

        h2d::Point2d pt = h2d::Point2d(x, z);
        if (pt.isInside(h2d::CPolyline(spg.points)))
        {
            return spg.soilcode;
        }
    }
    throw "point not in polygons or on the edge of a polygon";
}

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
        *sf = 9999;
        return;
    }
    if (res.size() < 2)
    {
        *sf = 9999;
        return;
    }

    auto pts = res.get();
    double entry_point = pts[0].getX();
    double exit_point = pts[1].getX();

    // cout << "x coordinate entry point: " << entry_point << endl;
    // cout << "x coordinate exit point : " << exit_point << endl;

    // # divide the space between the entry and exit point into evenly distributed points
    vector<double> slice_coordinates = {};
    double dl = exit_point - entry_point;
    for (int i = 0; i < NUM_S; ++i)
    {
        double x = entry_point + dl * double(i) / (NUM_S - 1);
        slice_coordinates.push_back(x);
    }

    // variables to store during calculation
    vector<int> M0 = {};    // index
    vector<double> M1 = {}; // b (width of slice)
    vector<double> M2 = {}; // W (weight of slice)
    vector<double> M3 = {}; // alpha
    vector<double> M4 = {}; // L
    vector<double> M5 = {}; // u
    vector<double> M6 = {}; // c
    vector<double> M7 = {}; // phi
    // FOR EACH SLICE
    for (int i = 1; i < slice_coordinates.size(); i++)
    {
        double x1 = slice_coordinates[i - 1];                        // left side of the slice
        double x2 = slice_coordinates[i];                            // right side of the slice
        double x3 = (x1 + x2) / 2.0;                                 // mid point of the slice
        double z1 = mz - sqrt(pow(radius, 2.0) - pow(mx - x1, 2.0)); // z value at bottom left of the slice
        double z2 = mz - sqrt(pow(radius, 2.0) - pow(mx - x2, 2.0)); // z value at bottom right of the slice
        double z3 = mz - sqrt(pow(radius, 2.0) - pow(mx - x3, 2.0)); // z value at bottom mid of the slice

        // cout << endl;
        // cout << "SLICE " << i << endl;
        // cout << "x slice left            : " << x1 << endl;
        // cout << "x slice middle          : " << x3 << endl;
        // cout << "x slice right           : " << x2 << endl;
        // cout << "bottom slice left       : " << z1 << endl;
        // cout << "bottom slice middle     : " << z3 << endl;
        // cout << "bottom slice right      : " << z2 << endl;

        // ready to create the slice polygon (clockwise)
        double arr[] = {x1, top, x2, top, x2, z2, x3, z3, x1, z1};
        vector<h2d::Point2d> slice_points = {h2d::Point2d(x1, top), h2d::Point2d(x2, top), h2d::Point2d(x2, z2), h2d::Point2d(x3, z3), h2d::Point2d(x1, z1)};
        PathsD cut_slice = h2dpoints_to_pathsd(slice_points);

        // now get all the intersections with the soilpolygons above the phreatic line
        // GET ALL SOILPOLYGONS THAT INTERSECT WITH THE SLICE
        vector<SoilPolygon> soilpolygons_above_pl = {};
        for (SoilPolygon spg : model.soilpolygons_above_pl)
        {
            PathsD polygon = soilpolygon_to_pathsd(spg);
            PathsD solution = Intersect(cut_slice, polygon, FillRule::NonZero);

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
        }

        vector<SoilPolygon> soilpolygons_below_pl = {};
        for (SoilPolygon spg : model.soilpolygons_below_pl)
        {
            PathsD polygon = soilpolygon_to_pathsd(spg);
            PathsD solution = Intersect(cut_slice, polygon, FillRule::NonZero);

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
        }

        // cout << "#soilpolygons above pl  : " << soilpolygons_above_pl.size() << endl;
        // cout << "#soilpolygons below pl  : " << soilpolygons_below_pl.size() << endl;

        double base_alpha = atan2((mz - z3), (mx - x3)) - 0.5 * PI;
        double b = x2 - x1;                           // width of the slice
        double base_L = b / cos(base_alpha);          // length at bottom of slice
        double u = u_at(x3, z3, model.phreatic_line); // waterpressure
        double W = 0;

        for (SoilPolygon spg : soilpolygons_above_pl)
        {
            h2d::CPolyline polygon = h2d::CPolyline(spg.points);
            try
            {
                Soil soil = get_soil(spg.soilcode, model.soils);
                W += polygon.area() * soil.y_dry;
            }
            catch (const char *e)
            {
                *sf = 9999;
                return;
            }
        }
        for (SoilPolygon spg : soilpolygons_below_pl)
        {
            h2d::CPolyline polygon = h2d::CPolyline(spg.points);
            try
            {
                Soil soil = get_soil(spg.soilcode, model.soils);
                W += polygon.area() * soil.y_sat;
            }
            catch (const char *e)
            {
                *sf = 9999;
                return;
            }
        }

        // cout << "base_alpha              : " << base_alpha << endl;
        // cout << "b                       : " << b << endl;
        // cout << "base_L                  : " << base_L << endl;
        // cout << "u                       : " << u << endl;
        // cout << "W                       : " << W << endl;

        // soil at the bottom of the slice
        string soilcode = get_soil_at(x3, z3, model.soilpolygons);
        Soil soil = get_soil(soilcode, model.soils);

        double c = soil.cohesion;
        double phi = soil.friction_angle;
        // cout << "soilcode                : " << soilcode << endl;
        // cout << "c                       : " << c << endl;
        // cout << "phi                     : " << phi << endl;

        M0.push_back(i);
        M1.push_back(b);
        M2.push_back(W);
        M3.push_back(-1.0 * base_alpha);
        M4.push_back(base_L);
        M5.push_back(u);
        M6.push_back(c);
        M7.push_back((phi / 180.0) * PI);
    }

    double denom = 0.0;
    for (int i = 0; i < M2.size(); i++)
    {
        denom += M2[i] * sin(M3[i]);
    }
    // cout << "denom                   : " << denom << endl;
    double cl = 0.0;
    for (int i = 0; i < M6.size(); i++)
    {
        cl += M6[i] * M4[i];
    }
    // cout << "cl                      : " << cl << endl;

    double isf = 1.0;
    int iteration = 0;

    while (1)
    {
        if (iteration >= MAX_I || isnan(isf))
        {
            isf = 9999.0;
            break;
        }

        vector<double> N = {};
        vector<double> N1 = {};
        vector<double> N2 = {};
        vector<double> N3 = {};

        // cout << "M" << endl;
        // for (int i = 0; i < M6.size(); i++)
        // {
        //     cout << M0[i] << " | " << M1[i] << " | " << M2[i] << " | " << M3[i] << " | " << M4[i] << " | " << M5[i] << " | " << M6[i] << " | " << M7[i] << endl;
        // }

        // cout << "N" << endl;
        for (int i = 0; i < M6.size(); i++)
        {
            // cout << M5[i] << " || " << M4[i] << " || " << M3[i] << " || " << M7[i] << endl;
            N1.push_back(M6[i] * M4[i] * sin(M3[i]));
            N2.push_back(M5[i] * M4[i] * sin(M3[i]) * tan(M7[i]));
            N3.push_back(cos(M3[i]) + (sin(M3[i]) * tan(M7[i])) / isf);
            N.push_back((M2[i] - (N1[i] - N2[i]) / isf) / N3[i]);

            // cout << i << " | " << N[i] << " | " << N1[i] << " | " << N2[i] << " | " << N3[i] << " | " << endl;
        }

        double d = 0.0;
        for (int i = 0; i < N.size(); i++)
        {
            d += (N[i] - M5[i] * M4[i]) * tan(M7[i]);
        }
        // cout << "d                       : " << d << endl;

        double fos = (cl + d) / denom;
        // cout << "iteration               : " << iteration << endl;
        // cout << "fos                     : " << fos << endl;
        if (abs(isf - fos) < 0.005)
            break;

        isf = (isf + fos) / 2.0;

        // cout << "isf                     : " << isf << endl;
        ++iteration;
    }

    // cout << "iteration               : " << iteration << endl;
    // cout << "sf                      : " << isf << endl;
    *sf = isf;
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

    if (!document.IsObject())
    {
        throw "Error reading document";
    }

    // PARSE SOILCOLLECTION
    vector<Soil> soils = {};
    Value &v_soilcollection = document["soilcollection"]["soils"];
    if (!v_soilcollection.IsArray())
    {
        throw "Error reading soilcollection";
    }

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
    if (!v_soilprofile2.IsArray())
    {
        throw "Error reading soilprofile2 soilpolygons";
    }
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
    if (!v_phreatic_line.IsArray())
    {
        throw "Error reading phreatic line";
    }
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
    if (!v_surface_line.IsArray())
    {
        throw "Error reading surface line";
    }
    for (SizeType i = 0; i < v_surface_line.Size(); i++)
    {
        double x = v_surface_line[i][0].GetDouble();
        double z = v_surface_line[i][1].GetDouble();
        spoints.push_back(h2d::Point2d{x, z});
    }
    h2d::OPolyline surface = h2d::OPolyline{spoints};

    // PARSE BISHOP SEARCH GRID
    Value &v_bishop_search_grid = document["bishop_search_grid"];
    double left = v_bishop_search_grid["left"].GetDouble();
    double bottom = v_bishop_search_grid["bottom"].GetDouble();
    double width = v_bishop_search_grid["width"].GetDouble();
    double height = v_bishop_search_grid["height"].GetDouble();
    int num_x = v_bishop_search_grid["num_x"].GetInt();
    int num_z = v_bishop_search_grid["num_z"].GetInt();
    int num_t = v_bishop_search_grid["num_tangent"].GetInt();
    double tangents_top = v_bishop_search_grid["tangents_top"].GetDouble();
    double tangents_bottom = v_bishop_search_grid["tangents_bottom"].GetDouble();
    double minimum_slip_plane_length = v_bishop_search_grid["minimum_slip_plane_length"].GetDouble();

    BishopSearchGrid bishop_search_grid = BishopSearchGrid{
        left,
        bottom,
        width,
        height,
        num_x,
        num_z,
        tangents_top,
        tangents_bottom,
        num_t,
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
        soilpolygons,
        soilpolygons_above_pl,
        soilpolygons_below_pl,
        bishop_search_grid,
        phreatic_line,
        surface,
    };
}

/*
Function that will calculate the Bishop slope stability safety factor based on
the given json string. The input json should comply with the leveelogic
CalculationModel.

The result will be a json string with the following information;
{
    "x":23.5,                   x coordinate of the slope circle
    "z":3.8,                    z coordinate of the slope circle
    "r":5.5,                    radius of the slope circle
    "sf":0.7095200041953091     safety factor (Bishop) // sf == 9999 means that errors were found
}
*/
const char *calculate_bishop(const string &json)
{
    // get the model from the string
    BishopModel model;
    model = parse_bishop_model(json);

    // DEBUG
    // model.print();
    // END DEBUG
    const int num_x = model.bishop_search_grid.num_x;
    const int num_z = model.bishop_search_grid.num_z;
    const int num_t = model.bishop_search_grid.num_tangent;
    const size_t N = num_t * num_x * num_z;

    double x = model.bishop_search_grid.left;
    double z = model.bishop_search_grid.bottom;
    double dx = model.bishop_search_grid.width / double(num_x - 1);
    double dz = model.bishop_search_grid.height / double(num_z - 1);
    double dt = (model.bishop_search_grid.tangents_top - model.bishop_search_grid.tangents_bottom) / double(num_t - 1);

    // temporary code to measure performance
    auto start = std::chrono::high_resolution_clock::now();

    // iterate over the possible slope circle locations, multithreaded
    // prepare the input

    vector<BishopResult> sfs;
    int i = 0;
    for (int nx = 0; nx < num_x; ++nx)
    {
        for (int nz = 0; nz < num_z; ++nz)
        {
            for (int nt = 0; nt < num_t; ++nt)
            {
                double cx = model.bishop_search_grid.left + nx * dx;
                double cz = model.bishop_search_grid.bottom + nz * dz;
                double t = model.bishop_search_grid.tangents_bottom + nt * dt;
                double r = cz - t;
                sfs.push_back(BishopResult{0.0, cx, cz, r});
                // threads.push_back(thread(sf_bishop, i, model, x, z, t, &sfs[i].sf)); this will create mem errors!
                ++i;
            }
        }
    }

    vector<thread> threads; // looks funny but avoids mem errors
    i = 0;
    for (int nx = 0; nx < num_x; ++nx)
    {
        for (int nz = 0; nz < num_z; ++nz)
        {
            for (int nt = 0; nt < num_t; ++nt)
            {
                double cx = model.bishop_search_grid.left + nx * dx;
                double cz = model.bishop_search_grid.bottom + nz * dz;
                double t = model.bishop_search_grid.tangents_bottom + nt * dt;
                double r = cz - t;
                threads.push_back(thread(sf_bishop, i, model, x, z, t, &sfs[i].sf));
                ++i;
            }
        }
    }

    // array<thread, N> threads;
    // array<BishopResult, N> sfs;
    // int i = 0;
    // for (int nx = 0; nx < num_x; ++nx)
    // {
    //     for (int nz = 0; nz < num_z; ++nz)
    //     {
    //         for (int nt = 0; nt < num_t; ++nt)
    //         {
    //             sfs[i].x = model.bishop_search_grid.left + nx * dx;
    //             sfs[i].z = model.bishop_search_grid.bottom + nz * dz;
    //             double t = model.bishop_search_grid.tangents_bottom + nt * dt;
    //             sfs[i].r = sfs[i].z - t;
    //             // todo error handling in thread
    //             threads[i] = thread(sf_bishop, i, model, x, z, t, &sfs[i].sf);
    //             ++i;
    //         }
    //     }
    // }

    for (auto &t : threads) // uncomment this to add threading again
    {
        t.join();
    }

    // temporary code to measure performance
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> elapsed = end - start;
    // std::cout << "Elapsed time " << elapsed.count() << " ms\n";

    BishopResult final_result = {9999, 0, 0, 0};
    for (auto &r : sfs)
    {
        if (r.sf > 0.0 && r.sf < final_result.sf)
        {
            final_result = r;
        }
    }
    std::cout << "result =" << final_result.sf << endl;

    Document d;
    d.SetObject();
    rapidjson::Document::AllocatorType &allocator = d.GetAllocator();
    d.AddMember("x", final_result.x, allocator);
    d.AddMember("z", final_result.z, allocator);
    d.AddMember("r", final_result.r, allocator);
    d.AddMember("sf", final_result.sf, allocator);

    StringBuffer buffer;
    Writer<StringBuffer> writer(buffer);
    d.Accept(writer);

    return buffer.GetString();
}

int main(int argc, char **argv)
{
    string json = "";

    if (argc == 2)
    {
        json = argv[1];
    }
    else
    {
        Document d;
        d.SetObject();
        rapidjson::Document::AllocatorType &allocator = d.GetAllocator();
        d.AddMember("error", "PStab requires a json string as the first and only argument", allocator);
        StringBuffer buffer;
        Writer<StringBuffer> writer(buffer);
        d.Accept(writer);
        cout << buffer.GetString();
        return 1;
    }

    string response;
    try
    {
        response = calculate_bishop(json);
    }
    catch (const char *error)
    {
        Document d;
        d.SetObject();
        rapidjson::Document::AllocatorType &allocator = d.GetAllocator();
        Value s;
        s.SetString(StringRef(error));
        d.AddMember("error", s, allocator);
        StringBuffer buffer;
        Writer<StringBuffer> writer(buffer);
        d.Accept(writer);
        cout << buffer.GetString();
        return 1;
    }
    cout << response;
    return 0;
}
