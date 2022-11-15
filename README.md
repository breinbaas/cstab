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

Windows keeps trying to use VStudio, use ```cmake .. -G "MinGW Makefiles"``` to use MinGW instead

## Developer

Rob van Putten, breinbaasnl@gmail.com 