cmake_minimum_required (VERSION 2.8.11)
project (tudelft_esri_geobim)

set(OCC_LIBRARY_NAMES TKernel TKMath TKBRep TKGeomBase TKGeomAlgo TKG3d TKG2d TKShHealing TKTopAlgo TKMesh TKOffset
                      TKPrim TKBool TKBO TKFillet TKSTEP TKSTEPBase TKSTEPAttr TKXSBase TKSTEP209 TKIGES TKHLR TKFeat)
                      
# IFCOPENSHELL_ROOT
# or
# BOOST_ROOT
# 
                      
OPTION(USE_STATIC_MSVC_RUNTIME "Link to the static runtime on MSVC." ON)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    add_definitions(-std=c++14)
else()
    message(FATAL_ERROR "requires a compiler with C++14 support")
endif()

if(UNIX)
    set(IFCOPENSHELL_INSTALL ${IFCOPENSHELL_ROOT}/build/Linux/x86_64/install)
    
    foreach(occt_dir occt-7.3.0 occt-7.3.0p3)
        if(EXISTS "${IFCOPENSHELL_INSTALL}/${occt_dir}/")
            set(OCCT_PATH "${IFCOPENSHELL_INSTALL}/${occt_dir}/")
        endif()
    endforeach()

    set(OCC_LIBRARY_DIR  ${OCCT_PATH}/lib)
    set(OCC_INCLUDE_DIR  ${OCCT_PATH}/include/opencascade)
    
    foreach(libname ${OCC_LIBRARY_NAMES})
        unset(lib CACHE)
    	find_library(lib NAMES "${libname}" PATHS "${OCC_LIBRARY_DIR}" NO_DEFAULT_PATH)
    	list(APPEND OCC_LIBRARIES "${lib}")
    endforeach()

    if(NOT APPLE)
        set(LIB_RT rt)
    endif()

    find_package(Threads)
    find_library(dl NAMES dl)

    set(OCC_LIBRARIES -Wl,--start-group ${OCC_LIBRARIES} -Wl,--end-group ${CMAKE_THREAD_LIBS_INIT} ${dl})
    
    add_definitions(-Wno-deprecated-declarations)
else()
    set(MSVC_YEAR 2017)
    if ("${MSVC_VERSION}" STREQUAL "1900")
        set(MSVC_YEAR 2015)
    endif()
    
    file(TO_CMAKE_PATH "${IFCOPENSHELL_ROOT}" IFCOPENSHELL_ROOT)
    if (EXISTS ${IFCOPENSHELL_ROOT}/deps-vs${MSVC_YEAR}-x86-installed/oce/Win32/lib)
        set(OCC_LIBRARY_DIR  ${IFCOPENSHELL_ROOT}/deps-vs${MSVC_YEAR}-x86-installed/oce/Win32/lib)    
    else()
        set(OCC_LIBRARY_DIR  ${IFCOPENSHELL_ROOT}/deps-vs${MSVC_YEAR}-x86-installed/opencascade-7.3.0/win32/lib)
    endif()
                         
    set(OCC_INCLUDE_DIR  ${IFCOPENSHELL_ROOT}/deps-vs${MSVC_YEAR}-x86-installed/oce/include/oce
                         ${IFCOPENSHELL_ROOT}/deps-vs${MSVC_YEAR}-x86-installed/opencascade-7.3.0/inc)
    
	add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS -D_UNICODE)
    
    string(REGEX REPLACE "([^;]+)" "${OCC_LIBRARY_DIR}/\\1.lib" OCC_LIBRARIES "${OCC_LIBRARY_NAMES}")
    
	set(Boost_USE_STATIC_LIBS   ON)
	set(Boost_USE_MULTITHREADED ON)

	if(USE_STATIC_MSVC_RUNTIME)
		set(Boost_USE_STATIC_RUNTIME ON)
    
		# Copied from IfcOpenShell to match compilation flags, todo detect whether necessary
		FOREACH(flag CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL
				CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
				CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
			IF(${flag} MATCHES "/MD")
				STRING(REGEX REPLACE "/MD" "/MT" ${flag} "${${flag}}")
			ENDIF()
			IF(${flag} MATCHES "/MDd")
				STRING(REGEX REPLACE "/MDd" "/MTd" ${flag} "${${flag}}")
			ENDIF()
		ENDFOREACH()
    endif()

    set(WS2_LIBRARIES ws2_32.lib)
endif()

set(BOOST_COMPONENTS regex program_options iostreams system)

# Find IfcOpenShell version
file(READ "${IFCOPENSHELL_ROOT}/src/ifcparse/IfcParse.h" header)
string(REGEX MATCH "#define IFCOPENSHELL_VERSION \"[0-9a-zA-Z\\.\\-]+\"" ifcopenshell_version "${header}")
string(REGEX MATCH "[0-9]\\.[0-9]" ifcopenshell_major_version "${ifcopenshell_version}")
message(STATUS "IfcOpenShell version ${ifcopenshell_major_version}")

set(IFC_LIBRARY_NAMES IfcGeom IfcParse geometry_kernel_cgal geometry_kernel_opencascade geometry_kernels geometry_mapping_ifc2x3 geometry_mapping_ifc4 geometry_mapping_ifc4x1 geometry_mapping_ifc4x2 geometry_mappings serializers serializers_ifc2x3 serializers_ifc4 serializers_ifc4x1 serializers_ifc4x2)

if (UNIX)
set(IFC_LIBRARY_NAMES ${IFC_LIBRARY_NAMES} ${IFC_LIBRARY_NAMES})
endif()

if(UNIX)
    set(BOOST_ROOT       ${IFCOPENSHELL_INSTALL}/boost-1.69.0)
    
    set(EIGEN_DIR        ${IFCOPENSHELL_INSTALL}/Eigen-3.3.7)
    
    set(CGAL_INCLUDE_DIR ${IFCOPENSHELL_INSTALL}/cgal-4.13/include)
    set(CGAL_LIBRARY_DIR ${IFCOPENSHELL_INSTALL}/cgal-4.13/lib)

    file(GLOB CGAL_LIBRARIES "${CGAL_LIBRARY_DIR}/*.a")
    
    set(GMP_INCLUDE_DIR  ${IFCOPENSHELL_INSTALL}/gmp-6.1.2/include)
    set(GMP_LIBRARY_DIR  ${IFCOPENSHELL_INSTALL}/gmp-6.1.2/lib)
    
    file(GLOB GMP_LIBRARIES "${GMP_LIBRARY_DIR}/*.a")
    
    set(MPFR_INCLUDE_DIR  ${IFCOPENSHELL_INSTALL}/mpfr-3.1.5/include)
    set(MPFR_LIBRARY_DIR  ${IFCOPENSHELL_INSTALL}/mpfr-3.1.5/lib)
    
    file(GLOB MPFR_LIBRARIES "${MPFR_LIBRARY_DIR}/*.a")

    set(IFC_INCLUDE_DIR  ${IFCOPENSHELL_ROOT}/src)
    set(IFC_LIBRARY_DIR  ${IFCOPENSHELL_INSTALL}/ifcopenshell/lib)
    
    string(REGEX REPLACE "([^;]+)" "${IFC_LIBRARY_DIR}/lib\\1.a" IFC_LIBRARIES "${IFC_LIBRARY_NAMES}")
else()
    set(BOOST_ROOT       ${IFCOPENSHELL_ROOT}/deps/boost_1_67_0)
    set(BOOST_LIBRARYDIR ${IFCOPENSHELL_ROOT}/deps/boost_1_67_0/stage/vs${MSVC_YEAR}-Win32/lib)
	
    set(IFC_INCLUDE_DIR  ${IFCOPENSHELL_ROOT}/src)
    set(IFC_LIBRARY_DIR  ${IFCOPENSHELL_ROOT}/build-vs${MSVC_YEAR}-x86/${CMAKE_BUILD_TYPE})
    
    string(REGEX REPLACE "([^;]+)" "${IFC_LIBRARY_DIR}/\\1.lib" IFC_LIBRARIES "${IFC_LIBRARY_NAMES}")
endif()

find_package(Boost REQUIRED COMPONENTS ${BOOST_COMPONENTS})

message(STATUS "Boost include files found in ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries found in ${Boost_LIBRARY_DIRS}")

include_directories(${OCC_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${IFC_INCLUDE_DIR}
${EIGEN_DIR} ${CGAL_INCLUDE_DIR} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(tudelft_esri_geobim geobim.cpp)
target_link_libraries(tudelft_esri_geobim ${IFC_LIBRARIES} ${OCC_LIBRARIES} ${CGAL_LIBRARIES}
    ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

install(TARGETS tudelft_esri_geobim RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
