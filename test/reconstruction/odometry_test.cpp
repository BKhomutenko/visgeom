// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"

#include "reconstruction/eucm_stereo.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "localization/photometric.h"

using namespace std;


// Helper function to read values from a file
template <int N>
array<double, N> readValues(std::ifstream& file)
{
	array<double, N> values;
	for (double & p: values)
	{
		file >> p;
	}
	file.ignore(65536, '\n'); //skip line
	return values;
}


int main(int argc, char** argv)
{
	//Open parameter file from input
	std::ifstream paramFile(argv[1]);
	if (not paramFile.is_open())
	{
		std::cerr << argv[1] << " : ERROR - file not found" << std::endl;
		return 1;
	}

	//Accept filename for images for left and right camera
	std::string datasetFoldername;
	std::string img1Filename, img2Filename; // Files that contain the image foldername as well as image+groundTruth data
	getline(paramFile, datasetFoldername);
	img1Filename = datasetFoldername + "/leftImageGroundTruth.txt";
	img2Filename = datasetFoldername + "/rightImageGroundTruth.txt";
	std::ifstream img1File(img1Filename);
	std::ifstream img2File(img2Filename);
	if( not( img1File.is_open() and img2File.is_open() ) )
	{
		std::cerr << "ERROR : Unable to open image list files" << std::endl;
		return 1;
	}
	std::string foldername1, foldername2; // Folders where the images are stored
	getline(img1File, foldername1);
	getline(img2File, foldername2);
	foldername1 = datasetFoldername + "/" + foldername1;
	foldername2 = datasetFoldername + "/" + foldername2;

	//Get first stereo pair (Mat8u) and store as keyframe
	std::string img1name, img2name;
	img1File >> img1name;
	img2File >> img2name;
	Mat8u keyframe1 = imread(datasetFoldername + "/" + img1name, 0);
	Mat8u keyframe2 = imread(datasetFoldername + "/" + img2name, 0);
	array<double, 6> prevPose1 = readValues<6>(img1File);
	array<double, 6> prevPose2 = readValues<6>(img1File);
	Transformation<double> T0key (prevPose1.data());
	std::cerr << "Attempting to open file:" << datasetFoldername << "/" << img1name << std::endl;
	if ( (keyframe1.data == NULL) or (keyframe2.data == NULL) )
	{
		std::cerr << "ERROR : Could not load first image pair" << std::endl;
		std::cerr << "Attempted to open file:" << datasetFoldername << "/" << img1name << std::endl;
		return 1;
	}

	//Accept EUCM params for first (left) camera, and create EnhancedCamera
	array<double, 6> cam1Params = readValues<6>(paramFile);
	EnhancedCamera camera1( cam1Params.data() );

	//Accept EUCM params for second (right) camera, and create Enhanced Camera
	array<double, 6> cam2Params = readValues<6>(paramFile);
	EnhancedCamera camera2( cam2Params.data() );

	//Accept pose transformation between stereo cameras
	array<double, 6> stereoPose = readValues<6>(paramFile);
	Transformation<double> T12( stereoPose.data() );

	//Accept parameters for stereo, and create StereoParameters, and then EnhancedStereo and MotionStereo
	StereoParameters stereoParams;
	stereoParams.verbosity = 0;
	stereoParams.salientPoints = false;
	paramFile >> stereoParams.u0;
	paramFile >> stereoParams.v0;
	paramFile >> stereoParams.dispMax;
	paramFile >> stereoParams.scale;
	paramFile.ignore();
	stereoParams.uMax = keyframe1.cols;
	stereoParams.vMax = keyframe1.rows;
	stereoParams.setEqualMargin();
	EnhancedStereo stereoSGM(T12, &camera1, &camera2, stereoParams);

	MotionStereoParameters mstereoParams;
	mstereoParams.verbosity = 0;
	mstereoParams.scale = stereoParams.scale;
	MotionStereo stereoMotion(&camera1, &camera2, mstereoParams);

	//Create ScaleParameters object for depthmap
	ScaleParameters scaleParams;
	scaleParams.scale = stereoParams.scale;
	scaleParams.u0 = 25;
	scaleParams.v0 = 25;
	scaleParams.uMax = keyframe1.cols;
	scaleParams.vMax = keyframe1.rows;
	scaleParams.setEqualMargin();

	//Accept camera pose wrt the robot, and create transformation from base to camera frame
	array<double, 6> cameraPose = readValues<6>(paramFile);
	Transformation<double> Tbase1(cameraPose.data());

	//Create ScalePhotometric object
	ScalePhotometric photometricLocalizer(5, &camera1);
	photometricLocalizer.setVerbosity(0);
	photometricLocalizer.computeBaseScaleSpace(keyframe1);

	//Call setBaseImage of MotionStereo to set keyframe as base image
	stereoMotion.setBaseImage(keyframe1);

	//Get initial SGM depthmap
	DepthMap keyDepth(&camera1, scaleParams);
	stereoSGM.computeStereo(keyframe1, keyframe2, keyDepth);

	int refinement = 0;

	while(1)
	{
		//Get next image pair to evaluate
		img1File >> img1name;
		img2File >> img2name;
		Mat8u newframe1 = imread(datasetFoldername + "/" + img1name, 0);
		Mat8u newframe2 = imread(datasetFoldername + "/" + img2name, 0);
		array<double, 6> newPose1 = readValues<6>(img1File);
		array<double, 6> newPose2 = readValues<6>(img2File);
		Transformation<double> T0new (newPose1.data());

		std::cout << "Attempting to open image: " << datasetFoldername << "/" << img1name << std::endl;
		if( (newframe1.data == NULL) or (newframe2.data == NULL) )
		{
			std::cout << "Next image could not be loaded. Terminating program." << std::endl;
			std::cout << "Attempted to open image: " << datasetFoldername << "/" << img1name << std::endl;
			return 0;
		}

		//Transformation between camera locations from camera at keyframe to camera right now
		Transformation<double> Tmotion = T0key.compose(Tbase1).inverseCompose(T0new.compose(Tbase1));


		//Localization ~ after 5 refinements
		if (refinement > 5)
		{
			//Set depthmap for localizer
			photometricLocalizer.depth() = keyDepth;

			//Use photometric localization to refine pose estimation
			photometricLocalizer.computePose(newframe1, Tmotion);
		}

		
		//Step: Use motion stereo to refine depthmap
		//Seed new depth estimation as equal to keyframe depthmap
		DepthMap newDepth = keyDepth.wrapDepth(Tmotion, scaleParams);

		//Compute new depth using computeDepth() of MotionStereo
		stereoMotion.reprojectDepth(Tmotion, newframe1, newDepth);

		//Merge the new depthmap into the keyframe depthmap
		keyDepth.merge(newDepth);

		//Output region
		Mat32f newDepthMat, keyDepthMat;
		newDepth.toMat(newDepthMat);
		keyDepth.toMat(keyDepthMat);
		cv::imshow("Current image", newframe1);
		cv::imshow("Motion stereo depthmap", newDepthMat);
		cv::imshow("Merged keyframe depthmap", keyDepthMat);
		std::cout << "Transformation: " << Tmotion << std::endl;
		cv::waitKey(1);

		refinement++; //increase the level of refinement of the keyframe depthmap
	}
	return 0;
}