## CStab

Developed by Rob van Putten | breinbaasnl@gmail.com | leveelogic@gmail.com

## What it is

CStab is a calculation kernel for slope stabiility assessments. By kernel I mean that there will never be user interface 
code in this repo but only a way to send a model to the programme which then returns the slope stability safety factor. 

![Slopecircle](https://github.com/breinbaas/cstab/blob/master/img/bishop.png?raw=true)

Since this programme is meant to be a calculation kernel for my [LeveeLogic web application](https://app.leveelogic.com) 
I expect a certain input which consists of a json file with soil and geometry information. 
You can find an example of the input in this repo under the test path. JSON can easily be edited and -even better- is
perfectly suitable to be used in web communication. A little more info on the used JSON file; the geometry is 
simply a collection of polygons with a certain soilcode. You can find examples in the ```soilpolygons``` section in
the json file. 

```
"soilpolygons": [
    {
        "points": [
            [
                15.0,
                0.0
            ],
            [
                17.0,
                2.0
            ],
            [
                20.0,
                2.0
            ],
            [
                20.0,
                0.0
            ]
        ],
        "soilcode": "clay"
    },
    ...
]
```

The soilcode refers to soil definitions in the ```soilcollection``` section. 

```
"soilcollection": {
    "soils": [
        {
            "code": "preexcavated",
            "color": "#54575c",
            "y_dry": 14.0,
            "y_sat": 14.0,
            "cohesion": 2.0,
            "friction_angle": 20.0
        },
        ...
    ]
}
```

I also implemented a ```phreatic_line``` and ```bishop_search_grid``` section for additional required information. Note that for now
you also have to add the ```surface_line``` in the json input file. In my case this is automatically generated out
of the polygons in my web application but I will probably move that calculation to this calculation kernel. 

You will also find some other input in the json file which simply is additional information that has been generated in
the LeveeLogic web application for other processes.

**Note** that the input will probably change over time since I have ideas and changing requirements for my web 
application but the essence (JSON and polygons) will stay the same.

**Note**
This program is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it
and/or modify it under the terms of the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.nl.html).

## Todo

**Currently not ready for release**

**Before release 0.1:**
* [ ] main should accept argument of type string
* [ ] we now save an array of sf but we want the min(sf) with the xm, zm and radius of the circle as a result
* [ ] benchmark with dstability
* [ ] cleanup code
* [ ] improve code (check referencing etc.)
* [ ] num x, z and tangents lines should be part of the input, not constants
* [ ] maybe? num slices should depend on the geometry but check first if this really makes a difference else keep constant
* [ ] error handling

**Later:**
* [ ] implement SHANSEP for undrained behaviour
* [ ] auto search grid
* [ ] optimize coordinate conversions if necessary (at initialization we use h2d and clipper, find out a way to avoid redundant conversions but only if this is a bottleneck which I doubt)
* [ ] surface calculation in c++ instead of demanding this is the json input (shouldn't be too hard with clipper2)
* [ ] Spencer model
* [ ] LiftVan model


## Dependencies

CStab makes use of the following libraries;

* [homog2d](https://github.com/skramm/homog2d) - homog2d, A single-file header-only C++ library dedicated to handling 2D lines, points and homographies (2D planar transformations), using (internally) homogeneous coordinates.
* [clipper2](http://www.angusj.com/clipper2/Docs/Overview.htm) - clipper2, an open source freeware library (written in C++, C# and Delphi Pascal) that performs line and polygon clipping, and offsetting.
* [rapidjson](https://rapidjson.org/) - rapidjson, A fast JSON parser/generator for C++ with both SAX/DOM style API

## License

GPLv3 so publish the code if you use / enhance it.

## Notes for myself

### CMake using mingw
Windows cmake keeps trying to use VStudio, use ```cmake .. -G "MinGW Makefiles"``` to use MinGW instead

### Testing as commandline programme
Argument to test as a commandline service;

* Windows


```
.\pstab.exe '{\"soilcollection\": {\"soils\": [{\"code\": \"preexcavated\", \"color\": \"#54575c\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 2.0, \"friction_angle\": 20.0}, {\"code\": \"unknown\", \"color\": \"#696969\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 2.0, \"friction_angle\": 20.0}, {\"code\": \"top_material\", \"color\": \"#696969\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 2.0, \"friction_angle\": 22.0}, {\"code\": \"bottom_material\", \"color\": \"#808080\", \"y_dry\": 17.0, \"y_sat\": 19.0, \"cohesion\": 0.0, \"friction_angle\": 30.0}, {\"code\": \"nl_veen\", \"color\": \"#786926\", \"y_dry\": 10.0, \"y_sat\": 10.0, \"cohesion\": 1.5, \"friction_angle\": 17.5}, {\"code\": \"nl_grof_zand\", \"color\": \"#faff00\", \"y_dry\": 19.0, \"y_sat\": 21.0, \"cohesion\": 0.0, \"friction_angle\": 35.0}, {\"code\": \"nl_middelgrof_zand\", \"color\": \"#e0e342\", \"y_dry\": 18.0, \"y_sat\": 20.0, \"cohesion\": 0.0, \"friction_angle\": 32.5}, {\"code\": \"nl_fijn_zand\", \"color\": \"#e6e876\", \"y_dry\": 17.0, \"y_sat\": 19.0, \"cohesion\": 0.0, \"friction_angle\": 30.0}, {\"code\": \"nl_kleiig_zand\", \"color\": \"#9fa12d\", \"y_dry\": 16.0, \"y_sat\": 18.0, \"cohesion\": 1.0, \"friction_angle\": 30.0}, {\"code\": \"nl_siltig_zand\", \"color\": \"#596b15\", \"y_dry\": 16.0, \"y_sat\": 17.0, \"cohesion\": 2.0, \"friction_angle\": 27.5}, {\"code\": \"nl_zandige_klei\", \"color\": \"#596b15\", \"y_dry\": 16.0, \"y_sat\": 16.0, \"cohesion\": 3.0, \"friction_angle\": 27.5}, {\"code\": \"nl_siltige_klei\", \"color\": \"#3c6318\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 4.0, \"friction_angle\": 25.0}, {\"code\": \"nl_klei\", \"color\": \"#39bf1b\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 5.0, \"friction_angle\": 25.0}, {\"code\": \"nl_venige_klei\", \"color\": \"#7d6b0f\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 3.0, \"friction_angle\": 20.0}, {\"code\": \"nl_humeuze_klei\", \"color\": \"#7d6b0f\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 3.0, \"friction_angle\": 20.0}, {\"code\": \"peat\", \"color\": \"#786926\", \"y_dry\": 10.0, \"y_sat\": 10.0, \"cohesion\": 1.5, \"friction_angle\": 17.5}, {\"code\": \"organic_clay\", \"color\": \"#a3de2f\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 2.0, \"friction_angle\": 20.0}, {\"code\": \"clay\", \"color\": \"#3c6318\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 5.0, \"friction_angle\": 25.0}, {\"code\": \"silty_clay\", \"color\": \"#596b15\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 3.0, \"friction_angle\": 25.0}, {\"code\": \"silty_sand\", \"color\": \"#9fa12d\", \"y_dry\": 16.0, \"y_sat\": 18.0, \"cohesion\": 1.0, \"friction_angle\": 27.5}, {\"code\": \"sand\", \"color\": \"#e6e876\", \"y_dry\": 17.0, \"y_sat\": 19.0, \"cohesion\": 0.0, \"friction_angle\": 30.0}, {\"code\": \"dense_sand\", \"color\": \"#fcf403\", \"y_dry\": 19.0, \"y_sat\": 21.0, \"cohesion\": 0.0, \"friction_angle\": 35.0}], \"aliases\": {}}, \"soilprofile2\": {\"soilpolygons\": [{\"points\": [[15.0, 0.0], [17.0, 2.0], [20.0, 2.0], [20.0, 0.0]], \"soilcode\": \"clay\"}, {\"points\": [[0.0, -2.0], [0.0, 0.0], [15.0, 0.0], [20.0, 0.0], [20.0, -0.5], [20.0, -1.8], [20.0, -2.0]], \"soilcode\": \"peat\"}, {\"points\": [[0.0, -5.0], [0.0, -2.0], [20.0, -2.0], [20.0, -4.5], [20.0, -5.0]], \"soilcode\": \"clay\"}, {\"points\": [[0.0, -6.2], [0.0, -5.0], [20.0, -5.0], [20.0, -6.2]], \"soilcode\": \"sand\"}, {\"points\": [[0.0, -8.5], [0.0, -6.2], [20.0, -6.2], [20.0, -6.4], [20.0, -8.5]], \"soilcode\": \"clay\"}, {\"points\": [[20.0, -15.0], [0.0, -15.0], [0.0, -8.5], [20.0, -8.5], [20.0, -8.9]], \"soilcode\": \"sand\"}, {\"points\": [[20.0, -0.5], [20.0, 0.0], [20.0, 2.0], [22.5, -0.5]], \"soilcode\": \"clay\"}, {\"points\": [[40.0, -1.8], [27.2, -1.8], [28.0, -1.0], [40.0, -1.0]], \"soilcode\": \"peat\"}, {\"points\": [[20.0, -1.8], [20.0, -0.5], [22.5, -0.5], [23.0, -1.0], [25.0, -1.0], [25.8, -1.8]], \"soilcode\": \"peat\"}, {\"points\": [[40.0, -4.5], [20.0, -4.5], [20.0, -2.0], [20.0, -1.8], [25.8, -1.8], [26.0, -2.0], [27.0, -2.0], [27.2, -1.8], [40.0, -1.8]], \"soilcode\": \"clay\"}, {\"points\": [[40.0, -6.4], [20.0, -6.4], [20.0, -6.2], [20.0, -5.0], [20.0, -4.5], [40.0, -4.5]], \"soilcode\": \"sand\"}, {\"points\": [[40.0, -8.9], [20.0, -8.9], [20.0, -8.5], [20.0, -6.4], [40.0, -6.4]], \"soilcode\": \"clay\"}, {\"points\": [[40.0, -8.9], [40.0, -15.0], [20.0, -15.0], [20.0, -8.9]], \"soilcode\": \"sand\"}], \"debug\": false}, \"characteristic_points\": [{\"l\": 17.0, \"point_type\": 20}, {\"l\": 17.5, \"point_type\": 30}, {\"l\": 19.5, \"point_type\": 31}, {\"l\": 25.0, \"point_type\": 50}, {\"l\": 28.0, \"point_type\": 51}, {\"l\": 23.0, \"point_type\": 40}], \"phreatic_line\": [[0.0, 1.6], [17.0, 1.6], [18.0, 1.6], [19.0, 0.6], [20.0, 0.5], [22.5, -0.6], [23.0, -1.1], [25.0, -1.3], [25.8, -1.3], [26.0, -1.3], [27.0, -1.3], [27.2, -1.3], [28.0, -1.3], [40.0, -1.3]], \"bishop_search_grid\": {\"left\": 23.5, \"bottom\": 3.8, \"width\": 10.0, \"height\": 5.0, \"num_x\": 10, \"num_z\": 10, \"num_tangent\": 5, \"tangents_top\": -1.7, \"tangents_bottom\": -3.7, \"minimum_slip_plane_length\": 3.0}, \"surface_line\": [[0.0, 0.0], [15.0, 0.0], [17.0, 2.0], [20.0, 2.0], [22.5, -0.5], [23.0, -1.0], [25.0, -1.0], [25.8, -1.8], [26.0, -2.0], [27.0, -2.0], [27.2, -1.8], [28.0, -1.0], [40.0, -1.0]]}'
```

* Faulty command Windows
```
.\pstab.exe '{\"soilcollection\": {\"0]]}'
```

* Missing soiltype
```
.\pstab.exe '{\"soilcollection\": {\"soils\": [{\"code\": \"preexcavated\", \"color\": \"#54575c\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 2.0, \"friction_angle\": 20.0}, {\"code\": \"unknown\", \"color\": \"#696969\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 2.0, \"friction_angle\": 20.0}, {\"code\": \"top_material\", \"color\": \"#696969\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 2.0, \"friction_angle\": 22.0}, {\"code\": \"bottom_material\", \"color\": \"#808080\", \"y_dry\": 17.0, \"y_sat\": 19.0, \"cohesion\": 0.0, \"friction_angle\": 30.0}, {\"code\": \"nl_veen\", \"color\": \"#786926\", \"y_dry\": 10.0, \"y_sat\": 10.0, \"cohesion\": 1.5, \"friction_angle\": 17.5}, {\"code\": \"nl_grof_zand\", \"color\": \"#faff00\", \"y_dry\": 19.0, \"y_sat\": 21.0, \"cohesion\": 0.0, \"friction_angle\": 35.0}, {\"code\": \"nl_middelgrof_zand\", \"color\": \"#e0e342\", \"y_dry\": 18.0, \"y_sat\": 20.0, \"cohesion\": 0.0, \"friction_angle\": 32.5}, {\"code\": \"nl_fijn_zand\", \"color\": \"#e6e876\", \"y_dry\": 17.0, \"y_sat\": 19.0, \"cohesion\": 0.0, \"friction_angle\": 30.0}, {\"code\": \"nl_kleiig_zand\", \"color\": \"#9fa12d\", \"y_dry\": 16.0, \"y_sat\": 18.0, \"cohesion\": 1.0, \"friction_angle\": 30.0}, {\"code\": \"nl_siltig_zand\", \"color\": \"#596b15\", \"y_dry\": 16.0, \"y_sat\": 17.0, \"cohesion\": 2.0, \"friction_angle\": 27.5}, {\"code\": \"nl_zandige_klei\", \"color\": \"#596b15\", \"y_dry\": 16.0, \"y_sat\": 16.0, \"cohesion\": 3.0, \"friction_angle\": 27.5}, {\"code\": \"nl_siltige_klei\", \"color\": \"#3c6318\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 4.0, \"friction_angle\": 25.0}, {\"code\": \"nl_klei\", \"color\": \"#39bf1b\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 5.0, \"friction_angle\": 25.0}, {\"code\": \"nl_venige_klei\", \"color\": \"#7d6b0f\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 3.0, \"friction_angle\": 20.0}, {\"code\": \"nl_humeuze_klei\", \"color\": \"#7d6b0f\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 3.0, \"friction_angle\": 20.0}, {\"code\": \"peat\", \"color\": \"#786926\", \"y_dry\": 10.0, \"y_sat\": 10.0, \"cohesion\": 1.5, \"friction_angle\": 17.5}, {\"code\": \"organic_clay\", \"color\": \"#a3de2f\", \"y_dry\": 14.0, \"y_sat\": 14.0, \"cohesion\": 2.0, \"friction_angle\": 20.0}, {\"code\": \"silty_clay\", \"color\": \"#596b15\", \"y_dry\": 15.0, \"y_sat\": 15.0, \"cohesion\": 3.0, \"friction_angle\": 25.0}, {\"code\": \"silty_sand\", \"color\": \"#9fa12d\", \"y_dry\": 16.0, \"y_sat\": 18.0, \"cohesion\": 1.0, \"friction_angle\": 27.5}, {\"code\": \"sand\", \"color\": \"#e6e876\", \"y_dry\": 17.0, \"y_sat\": 19.0, \"cohesion\": 0.0, \"friction_angle\": 30.0}, {\"code\": \"dense_sand\", \"color\": \"#fcf403\", \"y_dry\": 19.0, \"y_sat\": 21.0, \"cohesion\": 0.0, \"friction_angle\": 35.0}], \"aliases\": {}}, \"soilprofile2\": {\"soilpolygons\": [{\"points\": [[15.0, 0.0], [17.0, 2.0], [20.0, 2.0], [20.0, 0.0]], \"soilcode\": \"clay\"}, {\"points\": [[0.0, -2.0], [0.0, 0.0], [15.0, 0.0], [20.0, 0.0], [20.0, -0.5], [20.0, -1.8], [20.0, -2.0]], \"soilcode\": \"peat\"}, {\"points\": [[0.0, -5.0], [0.0, -2.0], [20.0, -2.0], [20.0, -4.5], [20.0, -5.0]], \"soilcode\": \"clay\"}, {\"points\": [[0.0, -6.2], [0.0, -5.0], [20.0, -5.0], [20.0, -6.2]], \"soilcode\": \"sand\"}, {\"points\": [[0.0, -8.5], [0.0, -6.2], [20.0, -6.2], [20.0, -6.4], [20.0, -8.5]], \"soilcode\": \"clay\"}, {\"points\": [[20.0, -15.0], [0.0, -15.0], [0.0, -8.5], [20.0, -8.5], [20.0, -8.9]], \"soilcode\": \"sand\"}, {\"points\": [[20.0, -0.5], [20.0, 0.0], [20.0, 2.0], [22.5, -0.5]], \"soilcode\": \"clay\"}, {\"points\": [[40.0, -1.8], [27.2, -1.8], [28.0, -1.0], [40.0, -1.0]], \"soilcode\": \"peat\"}, {\"points\": [[20.0, -1.8], [20.0, -0.5], [22.5, -0.5], [23.0, -1.0], [25.0, -1.0], [25.8, -1.8]], \"soilcode\": \"peat\"}, {\"points\": [[40.0, -4.5], [20.0, -4.5], [20.0, -2.0], [20.0, -1.8], [25.8, -1.8], [26.0, -2.0], [27.0, -2.0], [27.2, -1.8], [40.0, -1.8]], \"soilcode\": \"clay\"}, {\"points\": [[40.0, -6.4], [20.0, -6.4], [20.0, -6.2], [20.0, -5.0], [20.0, -4.5], [40.0, -4.5]], \"soilcode\": \"sand\"}, {\"points\": [[40.0, -8.9], [20.0, -8.9], [20.0, -8.5], [20.0, -6.4], [40.0, -6.4]], \"soilcode\": \"clay\"}, {\"points\": [[40.0, -8.9], [40.0, -15.0], [20.0, -15.0], [20.0, -8.9]], \"soilcode\": \"sand\"}], \"debug\": false}, \"characteristic_points\": [{\"l\": 17.0, \"point_type\": 20}, {\"l\": 17.5, \"point_type\": 30}, {\"l\": 19.5, \"point_type\": 31}, {\"l\": 25.0, \"point_type\": 50}, {\"l\": 28.0, \"point_type\": 51}, {\"l\": 23.0, \"point_type\": 40}], \"phreatic_line\": [[0.0, 1.6], [17.0, 1.6], [18.0, 1.6], [19.0, 0.6], [20.0, 0.5], [22.5, -0.6], [23.0, -1.1], [25.0, -1.3], [25.8, -1.3], [26.0, -1.3], [27.0, -1.3], [27.2, -1.3], [28.0, -1.3], [40.0, -1.3]], \"bishop_search_grid\": {\"left\": 23.5, \"bottom\": 3.8, \"width\": 10.0, \"height\": 5.0, \"num_x\": 10, \"num_z\": 10, \"num_tangent\": 5, \"tangents_top\": -1.7, \"tangents_bottom\": -3.7, \"minimum_slip_plane_length\": 3.0}, \"surface_line\": [[0.0, 0.0], [15.0, 0.0], [17.0, 2.0], [20.0, 2.0], [22.5, -0.5], [23.0, -1.0], [25.0, -1.0], [25.8, -1.8], [26.0, -2.0], [27.0, -2.0], [27.2, -1.8], [28.0, -1.0], [40.0, -1.0]]}'
```

## Developer

Rob van Putten, breinbaasnl@gmail.com 