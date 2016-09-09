// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"

#include "reconstruction/eucm_sgm.h"
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
		std::cout << " " << p;
	}
	file.ignore(5, '\n'); //skip line
	std::cout << std::endl;
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
	int skiprows;
	paramFile >> skiprows; paramFile.ignore();
	for(int i=1; i<=skiprows; i++)
	{
		std::string dump;
		getline(img1File, dump);
		getline(img2File, dump);
	}
	std::string keyimg1name, keyimg2name;
	img1File >> keyimg1name;
	img2File >> keyimg2name;
	Mat8u keyframe1 = imread(datasetFoldername + "/" + keyimg1name, 0);
	Mat8u keyframe2 = imread(datasetFoldername + "/" + keyimg2name, 0);
	std::cout << "Keyframe1 robot pose:" << std::endl;
	array<double, 6> prevPose1 = readValues<6>(img1File);
	std::cout << "Keyframe2 robot pose:" << std::endl;
	array<double, 6> prevPose2 = readValues<6>(img2File); //Unused
	Transformation<double> T0key (prevPose1.data());
	std::cerr << "Attempting to open image:" << datasetFoldername << "/" << keyimg1name << std::endl;
	std::cerr << "Attempting to open image:" << datasetFoldername << "/" << keyimg2name << std::endl;
	if ( (keyframe1.data == NULL) or (keyframe2.data == NULL) )
	{
		std::cerr << "ERROR : Could not load first image pair" << std::endl;
		std::cerr << "Attempted to open image:" << datasetFoldername << "/" << keyimg1name << std::endl;
		std::cerr << "Attempted to open image:" << datasetFoldername << "/" << keyimg2name << std::endl;
		return 1;
	}

	//Accept EUCM params for first (left) camera, and create EnhancedCamera
	std::cout << "EUCM Camera 1 Intrinsic parameters:" << std::endl;
	array<double, 6> cam1Params = readValues<6>(paramFile);
	EnhancedCamera camera1( cam1Params.data() );

	//Accept EUCM params for second (right) camera, and create Enhanced Camera
	std::cout << "EUCM Camera 2 Intrinsic parameters:" << std::endl;
	array<double, 6> cam2Params = readValues<6>(paramFile);
	EnhancedCamera camera2( cam2Params.data() );

	//Accept pose transformation between stereo cameras
	std::cout << "EUCM Stereo Camera Extrinsic parameters:" << std::endl;
	array<double, 6> stereoPose = readValues<6>(paramFile);
	Transformation<double> T12( stereoPose.data() );

	//Accept parameters for stereo, and create StereoParameters, and then EnhancedStereo and MotionStereo
	SGMParameters stereoSgmParams;
	stereoSgmParams.salientPoints = true;
	stereoSgmParams.verbosity = 0;
	paramFile >> stereoSgmParams.u0;
	paramFile >> stereoSgmParams.v0;
	paramFile >> stereoSgmParams.dispMax;
	paramFile >> stereoSgmParams.scale;
	paramFile.ignore();
	stereoSgmParams.uMax = keyframe1.cols;
	stereoSgmParams.vMax = keyframe1.rows;
	stereoSgmParams.setEqualMargin();
	EnhancedSGM stereoSGM(T12, &camera1, &camera2, stereoSgmParams);

	MotionStereoParameters mstereoParams;
	mstereoParams.verbosity = 0;
	mstereoParams.descLength = (stereoSgmParams.scale / 2 + 1) * 2 + 1;
	// mstereoParams.scale = stereoSgmParams.scale;
	mstereoParams.dispMax = 32;
	MotionStereo stereoMotion(&camera1, &camera2, mstereoParams);

	//Create ScaleParameters object for depthmap
	ScaleParameters scaleParams;
	scaleParams.scale = stereoSgmParams.scale;
	scaleParams.u0 = 25;
	scaleParams.v0 = 25;
	scaleParams.uMax = keyframe1.cols;
	scaleParams.vMax = keyframe1.rows;
	scaleParams.setEqualMargin();

	//Accept camera pose wrt the robot, and create transformation from base to camera frame
	std::cout << "Pose between robot and camera1:" << std::endl;
	array<double, 6> cameraPose = readValues<6>(paramFile);
	Transformation<double> Tbase1(cameraPose.data());

	//Create ScalePhotometric object
	ScalePhotometric photometricLocalizer(5, &camera1);
	photometricLocalizer.setVerbosity(0);
	photometricLocalizer.computeBaseScaleSpace(keyframe1);

	//Call setBaseImage of MotionStereo to set keyframe as base image
	stereoMotion.setBaseImage(keyframe1);

	//Get initial SGM depthmap
	DepthMap keyDepth(&camera1, scaleParams, 3); // Use 3-hypothesis depthmap
	stereoSGM.computeStereo(keyframe1, keyframe2, keyDepth);

	Mat32f keyDepthMap;
	keyDepth.toMat(keyDepthMap);
	cv::imshow("SGM Depthmap", keyDepthMap/20);
	cv::waitKey(0);

	int refinement = 0;

	while(1)
	{
		//Get next image pair to evaluate
		std::string img1name, img2name;
		img1File >> img1name;
		img2File >> img2name;
		Mat8u newframe1 = imread(datasetFoldername + "/" + img1name, 0);
		Mat8u newframe2 = imread(datasetFoldername + "/" + img2name, 0);
		std::cout << "New frame pose (cam1):" << std::endl;
		array<double, 6> newPose1 = readValues<6>(img1File);
		std::cout << "New frame pose (cam2):" << std::endl;
		array<double, 6> newPose2 = readValues<6>(img2File);
		Transformation<double> T0new (newPose1.data());

		std::cout << "Attempting to open image:" << datasetFoldername << "/" << img1name << std::endl;
		std::cout << "Attempting to open image:" << datasetFoldername << "/" << img2name << std::endl;
		if( (newframe1.data == NULL) or (newframe2.data == NULL) )
		{
			std::cout << "Next image could not be loaded. Terminating program." << std::endl;
			std::cout << "Attempted to open image:" << datasetFoldername << "/" << img1name << std::endl;
			std::cout << "Attempted to open image:" << datasetFoldername << "/" << img2name << std::endl;
			return 0;
		}

		//Transformation between camera locations from camera at keyframe to camera right now
		Transformation<double> Tmotion = T0key.compose(Tbase1).inverseCompose(T0new.compose(Tbase1));
		if(Tmotion.trans().squaredNorm() < 0.0001) continue;
		cout << "Motion:" << Tmotion << endl;


		//Localization ~ after 5 refinements
		if (refinement > 5)
		{
			cout << "Calculating photometric localization..." << endl;

			//Set depthmap for localizer
			photometricLocalizer.depth() = keyDepth;

			//Use photometric localization to refine pose estimation
			photometricLocalizer.computePose(newframe1, Tmotion);
		}

		
		//Step: Use motion stereo to refine depthmap
		//Seed new depth estimation as equal to keyframe depthmap
		cout << "Wrapping depthmap ..." << endl;
		DepthMap newDepth = keyDepth;

		//Compute new depth using computeDepth() of MotionStereo
		cout << "Calculating stereo from motion ..." << endl;
		stereoMotion.computeDepth(Tmotion, newframe1, newDepth);

		//Merge the new depthmap into the keyframe depthmap
		cout << "Merging depthmaps ..." << endl;
		keyDepth.merge(newDepth);

		//Output region
		Mat32f newDepthMat, keyDepthMat;
		newDepth.toMat(newDepthMat);
		keyDepth.toMat(keyDepthMat);
		cv::imshow("Current image", newframe1);
		cv::imshow("Motion stereo depthmap", newDepthMat/20);
		cv::imshow("Merged keyframe depthmap", keyDepthMat/20);
		if(cv::waitKey(0)==1048603) break; // Break on ESC key

		refinement++; //increase the level of refinement of the keyframe depthmap
	}
	return 0;
}