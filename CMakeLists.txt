cmake_minimum_required(VERSION 3.6)
project(maya-pointSetMatch)
set(PLUGIN_NAME pointSetMatch)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELEASE)

# Maya SDK
set(MAYA_ROOT "/usr/autodesk/maya2016")
set(MAYA_INCLUDE ${MAYA_ROOT}/include)
set(MAYA_LIB ${MAYA_ROOT}/lib)

# Lev-Mar
set(LEVMAR_ROOT "/home/davidc/bin/levmar/levmar-2.6")
set(LEVMAR_INCLUDE ${LEVMAR_ROOT})
set(LEVMAR_LIB ${LEVMAR_ROOT})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fPIC -std=c++0x")  # debug mode
## set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fPIC -std=c++0x")  # performance mode
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-Bsymbolic -shared")

set(UTILITIES_FILES
        include/utilities/debugUtils.h
        )

#
set(SOURCE_FILES
        ${UTILITIES_FILES}
        # include/testData1.h
        # include/testData2.h
        # include/testData3.h
        # include/testDataOrigin.h
        include/pointSetMatchUtils.h
        include/pointSetMatchNode.h
        # include/pointSetMatchCmd.h
        src/pointSetMatchNode.cpp
        # src/pointSetMatchCmd.cpp
        src/pointSetMatchMain.cpp
        )

include_directories(
        include
        ${LEVMAR_INCLUDE}
        ${MAYA_INCLUDE}
)

link_directories(
        ${LIBRARY_OUTPUT_PATH}
        ${LEVMAR_LIB}
        ${MAYA_LIB}
)

# 'PointSetMatch' maya plugin library
add_library(${PLUGIN_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${PLUGIN_NAME}
        OpenMaya
        OpenMayaAnim
        Foundation
        levmar
        lapack
        blas
        m)
set_target_properties(${PLUGIN_NAME} PROPERTIES
        PREFIX "" # no 'lib' prefix to .so files
        )