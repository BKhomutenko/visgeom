cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)

project( calibration )

find_package( OpenCV REQUIRED )
FIND_PACKAGE(Ceres REQUIRED)

#find_package( Eigen3 REQUIRED )

include_directories(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

include_directories(include)

add_executable( calibration
    src/calibration/calibration.cpp
)

target_link_libraries( calibration ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( calibration ${CERES_LIBRARIES})

add_executable( calibration_mei
    src/calibration/calibration_mei.cpp
)

target_link_libraries( calibration_mei ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( calibration_mei ${CERES_LIBRARIES})

add_executable( calibration_stereo
    src/calibration/calibration_stereo.cpp
)

target_link_libraries( calibration_stereo ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( calibration_stereo ${CERES_LIBRARIES})

add_executable( calibration_robot
    src/calibration/calibration_robot.cpp
)

target_link_libraries( calibration_robot ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( calibration_robot ${CERES_LIBRARIES})

add_executable( rectify
    src/calibration/rectify.cpp
)

target_link_libraries( rectify ${OpenCV_LIBRARIES} )
TARGET_LINK_LIBRARIES( rectify ${CERES_LIBRARIES})

add_executable( stereo 
   src/reconstruction/stereo_test.cpp
   src/reconstruction/eucm_stereo.cpp
   src/reconstruction/depth_map.cpp
)

target_link_libraries( stereo ${OpenCV_LIBS} )

add_executable( multistereo 
   src/reconstruction/stereo_multibase.cpp
   src/reconstruction/eucm_stereo.cpp
   src/reconstruction/depth_map.cpp
)

target_link_libraries( multistereo ${OpenCV_LIBS} )

add_executable( generate_plane 
   src/reconstruction/generate_plane.cpp
   src/reconstruction/eucm_stereo.cpp
   src/reconstruction/depth_map.cpp
)

target_link_libraries( generate_plane ${OpenCV_LIBS} )

add_executable( photometric
   src/localization/photometric_test.cpp
   src/localization/photometric.cpp
   src/reconstruction/eucm_stereo.cpp
   src/reconstruction/depth_map.cpp
)

target_link_libraries( photometric ${OpenCV_LIBS} )
TARGET_LINK_LIBRARIES( photometric ${CERES_LIBRARIES})

add_executable( reproject_plane
   src/reconstruction/reproject_test.cpp
   src/reconstruction/reproject_depth.cpp
   src/reconstruction/depth_map.cpp
)

target_link_libraries( reproject_plane ${OpenCV_LIBS} )


if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
