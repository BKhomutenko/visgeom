# Visual Geometry Projects
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
## Calibration


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
* **ucm** --- the Unified Camera Model [J. Courbon et. al., A generic fisheye camera model for robotic applications. IROS 2007]
* **eucm** --- the Enhanced Unified Camera Model [B. Khomutenko et. al., An enhanced unified camera model. IEEE RA-L 2016]
* **mei** --- the Unified Camera Model with distortion [C. Mei and P. Rives, Single view point omnidirectional camera calibration from planar grids. ICRA 2007]

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

## Technical Details

The transformations in the software are represented by translation vector and uTheta rotation vector (normalized vector represents the axis of rotation, the norm is the angle in radians)


