cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)

project( calibration )

find_package( OpenCV REQUIRED )
FIND_PACKAGE(Ceres REQUIRED)

#find_package( Eigen3 REQUIRED )

include_directories(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})

include_directories(include)

### CALIBRATION ###
#add_executable( calibration
#    test/calibration/calibration.cpp
#)

#target_link_libraries( calibration ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( calibration ${CERES_LIBRARIES})

#add_executable( calibration_mei
#    test/calibration/calibration_mei.cpp
#)

#target_link_libraries( calibration_mei ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( calibration_mei ${CERES_LIBRARIES})

#add_executable( calibration_stereo
#    test/calibration/calibration_stereo.cpp
#)

#target_link_libraries( calibration_stereo ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( calibration_stereo ${CERES_LIBRARIES})

#add_executable( calibration_robot
#    test/calibration/calibration_robot.cpp
#)

#target_link_libraries( calibration_robot ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( calibration_robot ${CERES_LIBRARIES})

#add_executable( rectify
#    test/calibration/rectify.cpp
#)

#target_link_libraries( rectify ${OpenCV_LIBRARIES} )
#TARGET_LINK_LIBRARIES( rectify ${CERES_LIBRARIES})

### RECONSTRUCTION ###

#add_executable( stereo 
#    test/reconstruction/stereo_test.cpp
#    src/reconstruction/eucm_stereo.cpp
#    src/reconstruction/eucm_epipolar.cpp
#    src/reconstruction/depth_map.cpp
#)

#target_link_libraries( stereo ${OpenCV_LIBS} )

#add_executable( multistereo 
#    test/reconstruction/stereo_multibase.cpp
#    src/reconstruction/eucm_stereo.cpp
#    src/reconstruction/eucm_epipolar.cpp
#    src/reconstruction/depth_map.cpp
#)

#target_link_libraries( multistereo ${OpenCV_LIBS} )

#add_executable( generate_plane 
#    test/reconstruction/generate_plane.cpp
#    src/reconstruction/eucm_stereo.cpp
#    src/reconstruction/eucm_epipolar.cpp
#    src/reconstruction/depth_map.cpp
#)

#target_link_libraries( generate_plane ${OpenCV_LIBS} )


add_executable( reproject_plane
    test/reconstruction/reproject_test.cpp
    src/reconstruction/depth_map.cpp
)

target_link_libraries( reproject_plane ${OpenCV_LIBS} )

add_executable( depth_reconstruct_test
    test/reconstruction/depth_reconstruct_test.cpp
    src/reconstruction/depth_map.cpp
)

target_link_libraries( depth_reconstruct_test ${OpenCV_LIBS} )

#add_executable( motion_test
#    test/reconstruction/motion_test.cpp
#    src/reconstruction/eucm_stereo.cpp
#    src/reconstruction/eucm_epipolar.cpp
#    src/reconstruction/depth_map.cpp
#)

#target_link_libraries( motion_test ${OpenCV_LIBS} )

add_executable( merge_test
    test/reconstruction/merge_test.cpp
    src/reconstruction/eucm_stereo.cpp
    src/reconstruction/eucm_epipolar.cpp
    src/reconstruction/depth_map.cpp
)

target_link_libraries( merge_test ${OpenCV_LIBS} )


### LOCALIZATION ###

#add_executable( photometric
#    test/localization/photometric_test.cpp
#    src/localization/photometric.cpp
#    src/localization/cost_functions.cpp
#    src/reconstruction/eucm_stereo.cpp
#    src/reconstruction/eucm_epipolar.cpp
#    src/reconstruction/depth_map.cpp
#)

#target_link_libraries( photometric ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( photometric ${CERES_LIBRARIES})

#add_executable( mi_test
#    test/localization/mi_test.cpp
#    src/localization/photometric.cpp
#    src/localization/cost_functions.cpp
#    src/reconstruction/eucm_stereo.cpp
#    src/reconstruction/eucm_epipolar.cpp
#    src/reconstruction/depth_map.cpp
#)

#target_link_libraries( mi_test ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( mi_test ${CERES_LIBRARIES})

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2 -Wno-deprecated")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
