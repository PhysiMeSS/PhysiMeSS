# PhysiMeSS unit tests

Tests of the geometry of the cell-fibre interaction (`addons/PhysiMeSS/PhysiMeSS_geometry.h`):
nearest point on a fibre, and rotation of a fibre about its centre or its crosslink point, in 2D and 3D.

The expected values come from the mathematics only (closed forms and properties such as
"the result is a rigid rotation" or "the result does not depend on the frame"), not from
simulation outputs. The tests only need the header and a C++11 compiler.

From the PhysiCell root directory:
```
make physimess-tests
```
or from this directory:
```
make run
```
The exit code is the number of failed tests.
