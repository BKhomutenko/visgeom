This project includes following features:
* Fast subpixel detector of calibration boards (checkerboards)
* Calibration toolbox for intrinsic, stereo and extrinsic calibration of fisheye cameras
* Algorithms of direct stereo correspondence
* Sparse and dense localization algorithms

Some of these features are quite mature and can be used for practical purposes. 
For example, calibration toolbox is automated and easy to use.
The part on 3D reconstruction requires a certain degree of refactoring.
The localization part is quite experimental.

# INSTALLATION

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
# CALIBRATION


To run the calibration process you need to collect calibration data and describe the problem in a .json file.
An example of such a file for monocular calibration is given in data/calib_example.json

## Calibration Variables

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

There is a single block of transformations **xiCamBoard**. 
**global** is __false__ because this transformation is unique for every image.
**constant** is __false__ because we have to compute it.
**prior** is __false__ because we don't have any prior estimation of these transformations.


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

In this case, we have one single camera called **camera1**. 
The following types are supported:
* ucm -- the Unified Camera Model 
* eucm -- the Enhanced Unified Camera Model
* mei -- the Unified Camera Model with 

The structure of the file is the following:
```
%NX %NY %W %EMAX %CheckExtraction
%Directory/
%File1
%File2
...
```
%NX and %NY define the size of the calibration pattern 
(number of internal corners along horizontal and vertical direction)

%W is the size of one square in the calibration pattern in meters

%EMAX defines the error threshold. It is used during the error analysis.
Program will show all the images with points that have reprojection error greater than this threshold
the points would be marked by dark circles; also it will output the image name. 
If %EMAX is 0, nothing will be shown

%CheckExtraction defines if the extracted patterns are shown to the user. If it is 0, the extraction is
automatic (pattern may be extracted incorrectly).
It it is 1, each extracted pattern is shown and the user can discard it pressing 'N' or accept 
it pressing any other key (e.g. 'space').

%Directory/ should contain all the files %File1 .. and all these files must be valid image files

Example of calling the calibration function:
```
$ ./calibration ex_calibrate.txt

Intrinsics :
    0.630112     1.02342     266.294     265.871     639.902     487.559
Ex = 0.169697; Ey = 0.161184; Emax = 1.79888
```
STEREO CALIBRATION
To run stereo calibration three files are required. Example:
```
$ ./calibration_stereo calibInfoLeft.txt calibInfoRight.txt calibInfoStereo.txt
```
calibInfoLeft.txt and calibInfoRight.txt are just the same as for the monocular calibration
The former provides images just for the left camera, the latter just for the right.
The special file that contain stereo information (calibInfoStereo.txt) has the following structure:
```
%NX %NY %W %EMAX %CheckExtraction
%Directory/
%LeftPrefix
%RightPrefix
%File1
%File2
...
```
The only difference is that the Directory has to contain a pair of images LefPrefixFile1 and RightPrefixFile1
that are images of the same scene taken from left and right cameras respectively.

In case of stereo calibration %W, %EMAX, and %CheckExtraction from the first file are used 
(the one that has stereo information).


RECOMMENDATIONS

First try to perform monocular calibration for each camera separately.
start always with %CheckExtraction = 1 and then just delete from the file names 
of images for which the extraction fails or is incorrect.
After you are sure that the extraction is correct put %CheckExtraction = 0
and future extractions will be completely automatic
Also start with %EMAX = 0 and see the 
standard deviation (Ex and Ey in the output). Then you can take for example 3*Ex as the threshold


RECTIFICATION AND UNDISTORTION

This program is used to check how the model cope with undistortion part.
Example:
$ ./rectify ex_undistort_rectify.txt

The structure of ex_undistort_rectify.txt:
```
%Intrinsic parameters
%u0 %v0 %f -- virtual pinhole camera parameters
%Rotation Vector
%Directory
%File1
%File2
...
```

TECHNICAL DETAILS

The transformations in the software are represented by translation vector and 
uTheta rotation vector (normalized vector represents the axis of rotation, the norm is the angle in radians)


