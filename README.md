## CStab

This is an temporary repo that will serve for the c++ code of the pstab Python package for slope stability assessments. 
The final c++ code will be moved to the pstab repository once it is finished.

The pstab Python package will enable developers to calculate slope stability safety factors. The c++ code is needed
to enable much faster calculation times. 

## Dependencies

CStab makes use of the following libraries;

* [homog2d](https://github.com/skramm/homog2d) - homog2d, A single-file header-only C++ library dedicated to handling 2D lines, points and homographies (2D planar transformations), using (internally) homogeneous coordinates.
* [clipper2](http://www.angusj.com/clipper2/Docs/Overview.htm) - clipper2, an open source freeware library (written in C++, C# and Delphi Pascal) that performs line and polygon clipping, and offsetting.

## License

GPLv3 so publish the code if you use / enhance it.

## Notes

Windows keeps trying to use VStudio, use ```cmake .. -G "MinGW Makefiles"``` to use MinGW instead

## Developer

Rob van Putten, breinbaasnl@gmail.com 