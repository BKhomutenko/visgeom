# Visual Geometry Project --- Fisheye Calibration
A C++ implementation of different aspects of geometric vision.
This project includes following features:
* Geometric framework including transformations, rotation parametrization and quaternions
* Fast subpixel detector of calibration boards (checkerboards)
* Calibration toolbox for intrinsic, stereo and extrinsic calibration of fisheye cameras
* Image synthesis for simulated-data-based tests
* Algorithms of direct stereo correspondence
* Sparse and dense localization algorithms

Some of these features are quite mature and can be used for practical purposes. For example, calibration toolbox is automated and easy to use. The part on 3D reconstruction requires a certain degree of refactoring. The localization part is quite experimental.

## Installation

Requred packages:
* ceres-solver 1.10.0+
* Eigen3
* OpenCV 2.4.9+

After you have installed the needed libraries, go to the directory visgeom in the terminal and run:
```
$ mkdir build
$ cd build
$ cmake ..
$ make 
```
## Monocular Calibration


To run the calibration process you need to collect calibration data and describe the problem in a .json file. An example of such a file for monocular calibration is given in data/calib_example.json

### Calibration Variables

First all the calibration variables have to be defined:

```
"transformations" : [
    {
        "name" : "xiCamBoard",
        "global" : false,
        "constant" : false,
        "prior" : false
    }
],
```

There is a single block of transformations **xiCamBoard**. The name is arbitrary.
- **global** is **false** because this transformation is unique for every image.
- **constant** is **false** because we have to compute it.
- **prior** is **false** because we don't have any prior estimation of these transformations.


```
"cameras": [
    {
       "name" : "camera1",
       "type" : "eucm",
       "constant" : false,
       "value" : [0.5, 1, 300, 300, 600, 500]
    }
],
```

In this case, we have one single camera called **camera1**. The following types are supported:
* **ucm** --- Unified Camera Model [J. Courbon et. al., A generic fisheye camera model for robotic applications. IROS 2007]
intrinsics: [xi, fu, fv, u0, v0]
* **eucm** --- Enhanced Unified Camera Model [B. Khomutenko et. al., An enhanced unified camera model. IEEE RA-L 2016]
intrinsics: [alpha, beta, fu, fv, u0, v0]
* **mei** --- Unified Camera Model with distortion [C. Mei and P. Rives, Single view point omnidirectional camera calibration from planar grids. ICRA 2007]
intrinsics: [xi, k1, k2, k3, k4, k5, fu, fv, u0, v0]

**value** contains the initial guess for the parameters, or if **constant** is **true**, then it is the intrinsic parameters which will be used during the whole optimization process. In the case of monocular calibration **constant** must be **false**.

### Dataset Description
A quick overview of different values:
#### type
Type of the data used in this block. For monocular calibration is always **images**
#### camera 
Camera name, should be declared before
### transform_chain
Describes the sequence of transformations between the camera and the calibration board. Each transformation is represented by a tuple:
*  **name** --- transformation declared before
*  **direct** --- **true** if the points must be transformed as X_cam = T(X_board). 

#### init
Defines which transformation from the transformation chain has to be initialized. In the case of monocular calibration it should be always the only transformation.
#### object
Calibration object description.
*  **type** --- the only supported type so far is **checkboard**,
*  **rows** and **cols** --- size of the calibration patterns. Only internal corners are counted
*  **size**  --- size of one square in meters

####  parameters  
Various flags to tweak the calibration process. If a particular flag is not in the list, then it is considered as false.
Hence undersore in front of a flag disables it
* **check_extraction** --- show the detected grid for visual examination
* **improve_detection** --- use subpixel detector
* **show_outliers** --- show the images with reprojected pattern if some of the errors are ubnormally high
* **save_outlire_images** --- save these images 
* **do_not_solve_global** --- don't solve the calibration problem, is used to see the initialization 
* **do_not_solve** --- don't solve the optimization for transformation initialization, is used to see the computed initial pose
#### images
contains the image names
* **prefix** --- concatenated to every element of **names**, can be simply "" (an empty line)
* **names** --- a list of file names


## Stereo Calibration

In the case of stereo calibration we will use three datasets. It combines two monocular calibration datasets 
and one stereo calibration dataset.

### Declaration of Transformations

In this case we have four transformations.

```
"transformations" : [
    {
        "name" : "xiCamBoard1",
        "global" : false,
        "constant" : false,
        "prior" : false
    },
    {
        "name" : "xiCamBoard2",
        "global" : false,
        "constant" : false,
        "prior" : false
    },
    {
        "name" : "xiCamBoardStereo",
        "global" : false,
        "constant" : false,
        "prior" : false
    },
    {
        "name" : "xiCam12",
        "global" : true,
        "constant" : false,
        "prior" : true,
        "value" : [0, -0.3, 0, 0, 0, 0]
    }
],
```

In this example (which is a typical example of stereo calibration definition) 
* **xiCamBoard1** and **xiCamBoard2** are used in monocular calibration subproblems. They are used exactly in the same manner as in monocular calibration.
* **xiCamBoardStereo** is used for either of stereo images.
* **xiCam12** is the same for all the stereo images and defines the transformation between the two cameras. It is **global** and has a **prior**.

### Transformation Prior Format

The following formats are supported to define a transformation 

* 3 values  [x, y, theta] --- 2D posture. z-coordinate, x- and y-rotation are 0
* 6 values  [x, y, z, rx, ry, rz] --- full 3D parametrized by translation and rotation vectors (see Rodriguez' rotation formula)
* 7 values  [x, y, z, qx, qy, qz, qw] ---  3D parametrized by a translation vector and a normalized quaternion
* 12 values --- [r11, r12, r13, t1, r21, r22, r23, t2, r31, r32, r33, t3] homogeneous transformation matrix, stored row-wise:
```
    [   r11     r12     r13     t1  ]
    [   r21     r22     r23     t2  ]
    [   r31     r32     r33     t3  ]
    [   0       0       0       1   ]
```

### Declaration of Cameras

In the case of stereo calibration, we'll have two cameras.
```
"cameras": [
    {
       "name" : "camera1",
       "type" : "eucm",
       "constant" : false,
       "value" : [0.55, 1.2, 200, 200,  320,  240]
    },
    {
       "name" : "camera2",
       "type" : "eucm",
       "constant" : false,
       "value" : [0.55, 1.2, 200, 200,  320,  240]
    }
],
```
### Data Definition

We will have four data blocks. The two monocular calibration blocks as well as the stereo calibration block for the first camera are exactly the same as for the monocular calibration.
The difference appears in the definition of the stereo calibration block for the second camera.
**The first thing to note** is that the images in the name list should correspond to the ones declared in the first stereo calibration block (they must be the same number and in the same order).
**The second thing** is that now the transformation chain includes one more transformation:
```
"transform_chain" : [
            {"name" : "xiCam12", "direct" : false},
            {"name" : "xiCamBoardStereo", "direct" : true}
        ]
        
```
The transformations are combined in the same order before being applied to the points of the calibration board.
In this case the points in the second camera frame will be defined as (pseudo-code):
```
X_2[i, j] = xiCam12.inverse().compose(xiCamBoardStereo[j]).transform(X_b[i])
```
* X_2[i, j] is the i-th point in the camera frame for j-th image
* X_b[i] it the i-th point in the board frame

Notice that for the first camera points will be defined as (in the same pseudo-code):
```
X_1[i, j] = xiCamBoardStereo[j].transform(X_b[i])
```


