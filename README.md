INSTALLATION
Requred packages:
* ceres-solver 1.10.0+
* Eigen3
* OpenCV 2.4.9+

After you have installed the needed libraries, go to the directory visgeom and run in bash:
```
$ cmake .
$ make 
```
CALIBRATION
(The file formatting will be changed to xml or yaml in the nearest future, the given version is temporary)
To run the calibration process you need to write calibration info files.
ex_calibrate.txt is an example of such file
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
$ ./calibration_stereo calibInfoStereo.txt calibInfoLeft.txt calibInfoRight.txt
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

In case of stereo calibration %W %EMAX %CheckExtraction from the first file are used 
(the one that has stereo information)


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


