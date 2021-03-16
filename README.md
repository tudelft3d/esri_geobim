GeoBIM application
==================

Building Information Modelling (BIM) is a fairly novel paradigm in the 
construction sector in which a multi-disciplinary information in 
construction is exchanged between stakeholders as a decomposition of 
individual elements. These elements are defined as a solid volume 
geometry that is procedurally defined, semantic information associated 
to elements as key-value pairs (property sets) and relationships between 
these elements such as decomposition relationships and wall and space 
boundary connectivity. 

Conversely in the field of Geographic Information Systems (GIS) a higher 
level view is required and there is a stronger focus on functional 
characteristics (for the purpose of e.g Facility Management) and spatial 
analysis. Due to these discrepancies in information management a direct 
conversion of an IFC building model to, for example, CityGML (the 
predominant standard for semantic city models) is suboptimal: the 
individual building elements do not form a complete manifold shell that 
separates interior and exterior and there is geometric informtion in the 
model that might be too detailed for this purpose. 

This application offers a transformation from IFC to CityGML (CityJSON 
to be precise) that takes into account these considerations. It is a 
reimplemetation of the work in: 

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
    #   either:
    python build-all.py IfcConvert
    #   or 
    CC=/usr/bin/clang-10 CXX=/usr/bin/clang++-10 python2 build-all.py IfcConvert
    cd ../../
    mkdir build
    cd build
    cmake .. -DIFCOPENSHELL_ROOT=`pwd`/../IfcOpenShell
    make
    ./tudelft_esri_geobim       

Instructions for Windows and for compilation using preexisting dependencies
will follow.

Features
--------

### Closing gaps in the facade

The formation of a complete manifold shell that separates interior from exterior depends on the individual building elements providing this separation. (Small) gaps in this facade are closed by means of the dilation radius that makes the invidual elements slightly larger so that they overlap. This dilation radius is (currently) a global value, so increasing this radius to close gaps might mean loss of geometric detail elsewhere. For this purpose there is a specific feature to overlay results for two radii and create IFC geometries that fill the gaps so that by incorporating these newly constructed elements a smaller global radius can be used preserving geometric detail.


### Gap radius finding using binary search


### Validating IsExternal property


Performance
-----------

CGAL number types

The robustness of the implementation in CGAL stems from what they call exact number types: trees of operands that keep their history and more precision bits than available in hardware. Interval types to encode the uncertainty as a pair of efficient machine native floating points and a fallback to the precise numbers when the interval arithmetic is uncertain.

Minkowski



Padding volume


2D

Openings posthoc

Vertex moving

Multi-threading

Reference counted mutable shared number types (trees of operands). Makes multi-threading difficult. Not only writing, but also reading, as a read can update the "lazy exact" representation. Extra care not to use static data fields such as a global Z axis.

Voxel preselection

Minkowski is O(V), voxel is O(m), so by using a coarse grid size the complexity can be greatly reduced. IFC interpretation could be re-used but is currently not implemented as the voxelization code base is based on the traditional IfcOpenShell version that relies on OpenCASCADE.


Results
-------

Duplex

Schependomlaan

Riverside LOD 300



