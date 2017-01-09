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
    src/calibration/calib_cost_functions.cpp
    src/calibration/trajectory_generation.cpp
)

target_link_libraries( calibration 
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_library( reconstruction STATIC     
    src/reconstruction/eucm_sgm.cpp  
    src/reconstruction/eucm_stereo.cpp
    src/reconstruction/eucm_motion_stereo.cpp
    src/reconstruction/eucm_epipolar.cpp
    src/reconstruction/depth_map.cpp
    src/reconstruction/triangulator.cpp
)

target_link_libraries( reconstruction ${OpenCV_LIBS} )

add_library( localization STATIC 
    src/localization/photometric.cpp
    src/localization/local_cost_functions.cpp
    src/localization/mono_odom.cpp
    src/localization/sparse_odom.cpp
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

add_executable( optim_trajectory
    test/calibration/trajectory.cpp
)

target_link_libraries( optim_trajectory
    calibration
    ${CERES_LIBRARIES}
)
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
    test/localization/mono_odom_test.cpp
)

target_link_libraries( odometry_test  
    reconstruction
    localization 
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_executable( sparse_odom_test
    test/localization/sparse_odom_test.cpp
)

target_link_libraries( sparse_odom_test 
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
    set(CMAKE_CXX_FLAGS "-Wno-deprecated -O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary

#    set(CMAKE_CXX_FLAGS "-Wno-deprecated -ggdb")        ## DEBUG    
endif()
