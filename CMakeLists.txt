cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)

project( calibration )

find_package( OpenCV REQUIRED )
FIND_PACKAGE( Ceres REQUIRED )
find_package( Boost COMPONENTS program_options REQUIRED )

#find_package( Eigen3 REQUIRED )

include_directories(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

include_directories(include)

#Lib compilation

add_library( calibration STATIC     
    src/calibration/unified_calibration.cpp
    src/calibration/cost_functors.cpp
)

target_link_libraries( calibration 
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_library( reconstruction STATIC     
    src/reconstruction/eucm_stereo.cpp
    src/reconstruction/eucm_sgm.cpp
    src/reconstruction/eucm_motion_stereo.cpp
    src/reconstruction/eucm_epipolar.cpp
    src/reconstruction/depth_map.cpp
)

target_link_libraries( reconstruction ${OpenCV_LIBS} )

add_library( localization STATIC 
    src/localization/photometric.cpp
    src/localization/cost_functions.cpp
    src/localization/mono_odom.cpp
)

TARGET_LINK_LIBRARIES( localization reconstruction  ${OpenCV_LIBS} ${CERES_LIBRARIES})

### CALIBRATION ###
add_executable( calib
    test/calibration/generic_calibration.cpp
)

target_link_libraries( calib
    calibration
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)



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

#add_executable( calibration_base
#    test/calibration/calibration_base.cpp
#)

#target_link_libraries( calibration_base ${OpenCV_LIBS} )
#TARGET_LINK_LIBRARIES( calibration_base ${CERES_LIBRARIES})

#add_executable( rectify
#    test/calibration/rectify.cpp
#)

#target_link_libraries( rectify ${OpenCV_LIBRARIES} )
#TARGET_LINK_LIBRARIES( rectify ${CERES_LIBRARIES})

### RECONSTRUCTION ###

add_executable( stereo 
    test/reconstruction/stereo_test.cpp
)

target_link_libraries( stereo 
    reconstruction
    ${OpenCV_LIBS} 
)

add_executable( epipolar 
    test/reconstruction/epipolar_trace.cpp
)

target_link_libraries( epipolar 
    reconstruction
    ${OpenCV_LIBS} 
)


add_executable( descriptor_step 
    test/reconstruction/descriptor_step.cpp
)

target_link_libraries( descriptor_step 
    reconstruction
    ${OpenCV_LIBS} 
)


add_executable( multistereo 
    test/reconstruction/stereo_multibase.cpp
)

target_link_libraries( multistereo 
    reconstruction
    ${OpenCV_LIBS} 
)

add_executable( generate_plane 
    test/reconstruction/generate_plane.cpp
)

target_link_libraries( generate_plane 
    reconstruction
    ${OpenCV_LIBS} 
)


add_executable( reproject_plane
    test/reconstruction/reproject_test.cpp
)

target_link_libraries( reproject_plane 
    reconstruction
    ${OpenCV_LIBS} 
)

add_executable( depth_reconstruct_test
    test/reconstruction/depth_reconstruct_test.cpp
)

target_link_libraries( depth_reconstruct_test 
    reconstruction
    ${OpenCV_LIBS} 
)

add_executable( motion_test
    test/reconstruction/motion_test.cpp
)

target_link_libraries( motion_test 
    reconstruction
    ${OpenCV_LIBS} 
)

add_executable( merge_test
    test/reconstruction/merge_test.cpp
)

target_link_libraries( merge_test 
    reconstruction
    ${OpenCV_LIBS} 
)

add_executable( motion_gradual
    test/reconstruction/motion_gradual.cpp
)

target_link_libraries( motion_gradual 
    reconstruction
    ${OpenCV_LIBS} 
)

## LOCALIZATION ###
add_executable( odometry_test
    test/localization/odometry_test.cpp
)

target_link_libraries( odometry_test  
    reconstruction
    localization 
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_executable( photometric
    test/localization/photometric_test.cpp
)

target_link_libraries( photometric  
    reconstruction
    localization 
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
)

add_executable( mi_test
    test/localization/mi_test.cpp
)

target_link_libraries( mi_test  
    reconstruction
    localization 
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2 -Wno-deprecated")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()
