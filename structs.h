#include <vector>
#include <string>
#include <stdexcept>

#include "clipper2/clipper.h"
#include "homog2d/homog2d.h"

using namespace std;
using namespace Clipper2Lib;

// Structures to save the data
// Soil information
struct Soil
{
    string code;
    string color;
    double y_dry;
    double y_sat;
    double cohesion;
    double friction_angle;
};

// The soilpolygons from the model
struct SoilPolygon
{
    vector<h2d::Point2d> points;
    string soilcode;
};

// The definition of the Bishop search grid
struct BishopSearchGrid
{
    double left;
    double bottom;
    double width;
    double height;
    int num_x;
    int num_z;
    double tangents_top;
    double tangents_bottom;
    int num_tangent;
    double minimum_slip_plane_length;
};

// One model to rule them all
class BishopModel
{
public:
    vector<Soil> soils;
    vector<SoilPolygon> soilpolygons;
    vector<SoilPolygon> soilpolygons_above_pl;
    vector<SoilPolygon> soilpolygons_below_pl;
    BishopSearchGrid bishop_search_grid;
    h2d::OPolyline phreatic_line;
    h2d::OPolyline surface;

    void print()
    {
        // DEBUG INFO
        cout << "BISHOP MODEL INFO" << endl;
        cout << "-----------------" << endl;
        cout << "Number of soilpolygons above phreatic line: " << soilpolygons_above_pl.size() << endl;
        for (SoilPolygon spg : soilpolygons_above_pl)
        {
            cout << "SOILPOLYGON " << spg.soilcode << endl;
            for (h2d::Point2d p : spg.points)
            {
                cout << "(" << p.getX() << ", " << p.getY() << ")" << endl;
            }
        }

        cout << "Number of soilpolygons below phreatic line: " << soilpolygons_below_pl.size() << endl;
        for (SoilPolygon spg : soilpolygons_below_pl)
        {
            cout << "SOILPOLYGON " << spg.soilcode << endl;
            for (h2d::Point2d p : spg.points)
            {
                cout << "(" << p.getX() << ", " << p.getY() << ")" << endl;
            }
        }
        cout << "-----------------" << endl;
        // END DEBUG INFO
    }
};

struct BishopResult
{
    double sf; // safety factor NOTE -9999 means no safety factor found / error in input or calculation
    double x;  // x of circle centre point
    double z;  // z of circle centre point
    double r;  // radius of the circle
};