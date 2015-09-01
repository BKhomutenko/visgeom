cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)

project( calibration )

find_package( OpenCV REQUIRED )

find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

include_directories(include)

add_executable( calibration
    src/calibration.cpp
)

target_link_libraries( calibration ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( calibration ${CERES_LIBRARIES})

add_executable( calibration_stereo
    src/calibration_stereo.cpp
)

target_link_libraries( calibration_stereo ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( calibration_stereo ${CERES_LIBRARIES})

add_executable( rectify
    src/rectify.cpp
)

target_link_libraries( rectify ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( rectify ${CERES_LIBRARIES})

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
