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

### LIB COMPILATION ###

#fisheye simulation
#TODO add multi-camera system
add_library( render STATIC   
    src/render/render.cpp
    src/render/plane.cpp
    src/render/background.cpp
    src/render/texture.cpp
    src/render/aa_filter_lt.cpp
)

target_link_libraries( render 
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_library( calibration STATIC   
    src/calibration/corner_detector.cpp  
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
    src/reconstruction/scale_parameters.cpp
    src/reconstruction/epipoles.cpp
)

target_link_libraries( reconstruction ${OpenCV_LIBS} )

add_library( localization STATIC 
    src/localization/photometric.cpp
    src/localization/local_cost_functions.cpp
    src/localization/mono_odom.cpp
    src/localization/sparse_odom.cpp
    src/localization/mapping.cpp
)

TARGET_LINK_LIBRARIES( localization
    reconstruction
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
)

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


add_executable( rectify
    test/calibration/rectify.cpp
)

target_link_libraries( rectify
    calibration
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_executable( odom_interpolation
    test/calibration/odometry_interpolation.cpp
)

target_link_libraries( odom_interpolation
    ${Boost_LIBRARIES}
)

add_executable( stereo_sync
    test/calibration/stereo_sync.cpp
)

target_link_libraries( stereo_sync
    ${Boost_LIBRARIES}
)

add_executable( optim_trajectory
    test/calibration/trajectory.cpp
)

target_link_libraries( optim_trajectory
    calibration
    ${CERES_LIBRARIES}
)

add_executable(board_detector
    test/calibration/board_detector.cpp
)

target_link_libraries( board_detector
    calibration
    ${OpenCV_LIBS}  
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

### RECONSTRUCTION ###

add_executable( equalizer_test 
    test/reconstruction/equalizer_test.cpp
)

target_link_libraries( equalizer_test
    ${OpenCV_LIBS} 
)

add_executable( renderer_test 
    test/reconstruction/renderer_test.cpp
)

target_link_libraries( renderer_test 
    render
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

add_executable( stereo_test
    test/reconstruction/stereo_test.cpp
)

target_link_libraries( stereo_test
    reconstruction
    render
    ${OpenCV_LIBS} 
)

## LOCALIZATION ###
add_executable( mapping
    test/localization/map_test.cpp
)

target_link_libraries( mapping  
    reconstruction
    localization
    render
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_executable( map_real_data
    test/localization/map_real_data.cpp
)

target_link_libraries( map_real_data  
    reconstruction
    localization
    render
    ${OpenCV_LIBS}
    ${CERES_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_executable( odometry_test
    test/localization/mono_odom_test.cpp
)

target_link_libraries( odometry_test  
    reconstruction
    localization
    render
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

#    set(CMAKE_CXX_FLAGS "-Wno-deprecated -ggdb")        # DEBUG    
endif()
