GeoBIM application
==================

This application is a reimplemetation of the work in:

> Donkers, S., Ledoux, H., Zhao, J., & Stoter, J. (2016). Automatic conversion of IFC datasets to
> geometrically and semantically correct CityGML LOD3 buildings. Transactions in GIS, 20(4), 547-569.

In short, interpreted geometries in IFC building models are enlarged
(dilated) as the Minkowski sum of the IFC volumes and a constant size
polyhedron (small cube in our case for now). The resulting geometries are
all accumulated as the Boolean sum of the enlarged element volumes.

A specific branch of IfcOpenShell is used that can use CGAL to process the
element geometries. The application uses the Nef Polyhedron 3 and Minkowski
sum packages from CGAL.

Compilation
-----------

The application relies on a subset of the dependencies of IfcOpenShell.
Therefore easiest is to install IfcOpenShell first. The CMake file in this
repository can identify the dependencies as installed by the IfcOpenShell
build script. This makes compilation of this executable a fairly trivial,
but lengthy process (the compilation of all IfcOpenShell dependencies takes
about an hour on average hardware, five minutes on 32 core machine).

    git clone --recursive https://github.com/tudelft3d/esri_geobim
    cd esri_geobim/IfcOpenShell/nix
    apt-get install git gcc g++ autoconf bison make \
        libfreetype6-dev mesa-common-dev libffi-dev cmake
    python build-all.py IfcConvert
    cd ../../
    mkdir build
    cd build
    cmake .. -DIFCOPENSHELL_ROOT=`pwd`/../IfcOpenShell
    make
    ./tudelft_esri_geobim       

Instructions for Windows and for  compilation using preexisting dependencies
will follow.
