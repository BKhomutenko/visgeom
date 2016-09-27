// A test to see how we can perform combined localization and reconstruction
// Uses the data recorded from the Fluence

#include "io.h"
#include "ocv.h"
#include "timer.h"

#include <boost/program_options.hpp>

#include "reconstruction/eucm_sgm.h"
#include "reconstruction/eucm_motion_stereo.h"
#include "localization/photometric.h"


using namespace std;
namespace po = boost::program_options;

const float imgIntensityScale = 0.2;

// Helper function to read values from a file
template <int N>
array<double, N> readValues(std::ifstream& file, bool output = false)
{
	array<double, N> values;
	for (double & p: values)
	{
		file >> p;
		if(output) std::cout << " " << p;
	}
	file.ignore(5, '\n'); //skip line
	if(output) std::cout << std::endl;
	return values;
}

template <int N>
array<double, N> readValues(const po::variables_map & vmap, const std::string & key, bool output = false)
{
	vector<double> values = vmap[key].as< vector<double> >();
	if ( values.size()!=N )
	{
		std::cerr << "ERROR : Option " << key << " does not have " << N << "arguments" << endl;
		std::exit(1);
	}
	array<double, N> valArray;
	for (int i = 0; i < N; ++i)
	{
		valArray[i] = values[i];
		if(output) std::cout << " " << valArray[i];
	}
	if(output) std::cout << std::endl;
	return valArray;
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

	Timer timer;

	po::options_description options("Odometry Test options");
	options.add_options() ("help", "Produce help message");
	options.add_options() ("foldername", po::value<string>()->required(), "Folder to search for the dataset");
	options.add_options() ("skipRows", po::value<int>()->required(), "Number of rows of images to skip in the Image files");
	options.add_options() ("refinement", po::value<int>()->required(), "Refinement level threshold before photometric localization");

	options.add_options() ("Intrinsics.Camera1", po::value< vector<double> >()->required(), "EUCM Intrinsics for camera 1 (6 doubles)");
	options.add_options() ("Intrinsics.Camera2", po::value< vector<double> >()->required(), "EUCM Intrinsics for camera 2 (6 doubles)");
	options.add_options() ("Extrinsics.Stereo", po::value< vector<double> >()->required(), "EUCM Stereo extrinsic calibration (6 doubles)");
	
	options.add_options() ("StereoSGMParams.u0", po::value<double>()->required(), "Stereo SGM Parameters - u0");
	options.add_options() ("StereoSGMParams.v0", po::value<double>()->required(), "Stereo SGM Parameters - v0");
	options.add_options() ("StereoSGMParams.dispMax", po::value<double>()->required(), "Stereo SGM Parameters - dispMax");
	options.add_options() ("StereoSGMParams.scale", po::value<double>()->required(), "Stereo SGM Parameters - scale");
	options.add_options() ("StereoSGMParams.verbosity", po::value<int>()->required(), "Stereo SGM Parameters - verbosity level");
	options.add_options() ("StereoSGMParams.salientPoints", po::value<bool>()->required(), "Stereo SGM Parameters - use only salient points");
	
	options.add_options() ("MotionStereoParams.descLength", po::value<int>()->required(), "Motion Stereo descriptor length in pixels (needs to be odd)");
	options.add_options() ("MotionStereoParams.dispMax", po::value<int>()->required(), "Motion Stereo maximum disparity");
	options.add_options() ("MotionStereoParams.scale", po::value<int>()->required(), "Motion Stereo scale - unused?");
	options.add_options() ("MotionStereoParams.verbosity", po::value<int>()->required(), "Motion Stereo verbosity level");

	options.add_options() ("Pose.Base2Camera", po::value< vector<double> >()->required(), "6dof pose between robot base and camera 1 (6 doubles)");
	//TODO

	po::variables_map vmap;
	po::store( po::parse_config_file(paramFile, options, true), vmap);
	
	if( vmap.count("help") )
	{
		cout << options << endl;
		return 1;
	}
	
	po::notify( vmap ); // Checks if all required values have been defined


	//Accept filename for images for left and right camera
	std::string datasetFoldername;
	std::string img1Filename, img2Filename; // Files that contain the image foldername as well as image+groundTruth data
	datasetFoldername = vmap["foldername"].as<string>();
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
	std::cout << "Folder for camera1 : " << foldername1 << std::endl;
	std::cout << "Folder for camera2 : " << foldername2 << std::endl;

	//Skip the first few images of the dataset
	int skiprows = vmap["skipRows"].as<int>();
	for(int i=1; i<=skiprows; i++)
	{
		std::string dump;
		getline(img1File, dump);
		getline(img2File, dump);
	}

	//Get first stereo pair (Mat8u) and store as keyframe
	std::string keyimg1name, keyimg2name;
	img1File >> keyimg1name;
	img2File >> keyimg2name;
	Mat8u keyframe1 = imread(datasetFoldername + "/" + keyimg1name, 0);
	Mat8u keyframe2 = imread(datasetFoldername + "/" + keyimg2name, 0);
	std::cout << "Keyframe1 robot pose:" << std::endl;
	array<double, 6> prevPose1 = readValues<6>(img1File, true);
	std::cout << "Keyframe2 robot pose:" << std::endl;
	array<double, 6> prevPose2 = readValues<6>(img2File, true); //Unused
	Transformation<double> T0key (prevPose1.data());
	std::cout << "Attempting to open image:" << datasetFoldername << "/" << keyimg1name << std::endl;
	std::cout << "Attempting to open image:" << datasetFoldername << "/" << keyimg2name << std::endl;
	if ( (keyframe1.data == NULL) or (keyframe2.data == NULL) )
	{
		std::cerr << "ERROR : Could not load first image pair" << std::endl;
		std::cerr << "Attempted to open image:" << datasetFoldername << "/" << keyimg1name << std::endl;
		std::cerr << "Attempted to open image:" << datasetFoldername << "/" << keyimg2name << std::endl;
		return 1;
	}
	cv::imshow("Keyframe", keyframe1);

	//Accept EUCM params for first (left) camera, and create EnhancedCamera
	std::cout << "EUCM Camera 1 Intrinsic parameters:" << std::endl;
	array<double, 6> cam1Params = readValues<6>(vmap, "Intrinsics.Camera1", true);
	EnhancedCamera camera1( cam1Params.data() );

	//Accept EUCM params for second (right) camera, and create Enhanced Camera
	std::cout << "EUCM Camera 2 Intrinsic parameters:" << std::endl;
	array<double, 6> cam2Params = readValues<6>(vmap, "Intrinsics.Camera2", true);
	EnhancedCamera camera2( cam2Params.data() );

	//Accept pose transformation between stereo cameras
	std::cout << "EUCM Stereo Camera Extrinsic parameters:" << std::endl;
	array<double, 6> stereoPose = readValues<6>(vmap, "Extrinsics.Stereo", true);
	Transformation<double> T12( stereoPose.data() );

	//Accept parameters for stereo, and create StereoParameters, and then EnhancedStereo and MotionStereo
	SGMParameters stereoSgmParams;
	stereoSgmParams.salientPoints = vmap["StereoSGMParams.salientPoints"].as<bool>();
	stereoSgmParams.verbosity = vmap["StereoSGMParams.verbosity"].as<int>();
	stereoSgmParams.u0 = vmap["StereoSGMParams.u0"].as<double>();
	stereoSgmParams.v0 = vmap["StereoSGMParams.v0"].as<double>();
	stereoSgmParams.dispMax = vmap["StereoSGMParams.dispMax"].as<double>();
	stereoSgmParams.scale = vmap["StereoSGMParams.scale"].as<double>();
	std::cout << "Accepted parameters: u0 = " << stereoSgmParams.u0 << "; v0 = " << stereoSgmParams.v0 << "; dispMax = " << stereoSgmParams.dispMax << "; scale = " << stereoSgmParams.scale << std::endl;
	stereoSgmParams.uMax = keyframe1.cols;
	stereoSgmParams.vMax = keyframe1.rows;
	stereoSgmParams.setEqualMargin();
	// stereoSgmParams.setXMargin(stereoSgmParams.u0);
	stereoSgmParams.setYMargin(330);
	stereoSgmParams.hypMax = 1;
	EnhancedSGM stereoSGM(T12, &camera1, &camera2, stereoSgmParams);

	MotionStereoParameters mstereoParams;
	mstereoParams.verbosity = vmap["MotionStereoParams.verbosity"].as<int>();
	mstereoParams.descLength =  vmap["MotionStereoParams.descLength"].as<int>(); //(stereoSgmParams.scale / 2 + 1) * 2 + 1;
	mstereoParams.scale = vmap["MotionStereoParams.scale"].as<int>();
	mstereoParams.dispMax = vmap["MotionStereoParams.dispMax"].as<int>();
	MotionStereo stereoMotion(&camera1, &camera1, mstereoParams);

	//Create ScaleParameters object for depthmap
	// ScaleParameters scaleParams;
	// scaleParams.scale = stereoSgmParams.scale;
	// scaleParams.u0 = 25;
	// scaleParams.v0 = 25;
	// scaleParams.uMax = keyframe1.cols;
	// scaleParams.vMax = keyframe1.rows;
	// // scaleParams.setEqualMargin();
	// scaleParams.setXMargin(25);
	// scaleParams.setYMargin(300);

	//Accept camera pose wrt the robot, and create transformation from base to camera frame
	std::cout << "Pose between robot and camera1:" << std::endl;
	array<double, 6> cameraPose = readValues<6>(vmap, "Pose.Base2Camera", true);
	Transformation<double> Tbase1(cameraPose.data());

	//Create ScalePhotometric object
	ScalePhotometric photometricLocalizer(5, &camera1);
	photometricLocalizer.setVerbosity(0);
	photometricLocalizer.computeBaseScaleSpace(keyframe1);

	//Call setBaseImage of MotionStereo to set keyframe as base image
	stereoMotion.setBaseImage(keyframe1);

	//Get initial SGM depthmap
	DepthMap keyDepth; // Use 3-hypothesis depthmap
	stereoSGM.computeStereo(keyframe1, keyframe2, keyDepth);

	Mat32f keyDepthMap;
	keyDepth.toInverseMat(keyDepthMap);
	cv::namedWindow("SGM Depthmap", 1);
	cv::imshow("SGM Depthmap", keyDepthMap / imgIntensityScale);
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
		// std::cout << "New frame pose (cam1):" << std::endl;
		array<double, 6> newPose1 = readValues<6>(img1File, false);
		// std::cout << "New frame pose (cam2):" << std::endl;
		array<double, 6> newPose2 = readValues<6>(img2File, false);
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
		cout << "Motion (input):" << Tmotion << endl;


		//Images for photometric residual
		Mat32f wrappedImg;
		//Localization ~ after 5 refinements
		if (refinement >= vmap["refinement"].as<int>())
		{
			// cout << "Calculating photometric localization..." << endl;

			//Set depthmap for localizer
			photometricLocalizer.depth() = keyDepth;

			//Calculate photometric residual
			Mat32f keyImg_2, residualImg2;
			keyframe1.copyTo(keyImg_2);
			photometricLocalizer.wrapImage(newframe1, wrappedImg, Tmotion);
			residualImg2 = keyImg_2 - wrappedImg + 128;
			cv::threshold(wrappedImg, wrappedImg, 0, 1/255.0, cv::THRESH_BINARY);
			cv::imshow("Photometric Residual prev", residualImg2.mul(wrappedImg) );

			//Use photometric localization to refine pose estimation
			photometricLocalizer.computePose(newframe1, Tmotion);
			cout << "Motion (photo):" << Tmotion << endl;

			//Calculate photometric residual
			Mat32f keyImg_1, residualImg;
			keyframe1.copyTo(keyImg_1);
			photometricLocalizer.wrapImage(newframe1, wrappedImg, Tmotion);

			residualImg = keyImg_1 - wrappedImg + 128;
			cv::threshold(wrappedImg, wrappedImg, 0, 1/255.0, cv::THRESH_BINARY);
			cv::imshow("Photometric Residual", residualImg.mul(wrappedImg) );
		}

		
		//Step: Use motion stereo to refine depthmap
		//Seed new depth estimation as equal to keyframe depthmap
		DepthMap newDepth = keyDepth; 
//		newDepth.setTo(252, 100);
//		newDepth.setDefault();
		//Compute new depth using computeDepth() of MotionStereo
		// cout << "Calculating stereo from motion ..." << endl;
		stereoMotion.computeDepth(Tmotion, newframe1, newDepth);

		//Merge the new depthmap into the keyframe depthmap
		// cout << "Merging depthmaps ..." << endl;
		keyDepth.merge(newDepth);
		// keyDepth.filterNoise();

		//Output region
		Mat32f newDepthMat, keyDepthMat;
		newDepth.toInverseMat(newDepthMat);
		keyDepth.toInverseMat(keyDepthMat);
		cv::imshow("Current image", newframe1);
		cv::imshow("Motion stereo depthmap", newDepthMat / imgIntensityScale);
		cv::imshow("Merged keyframe depthmap", keyDepthMat / imgIntensityScale);
		cv::imshow("diff", (newDepthMat - keyDepthMap)*25 + 0.5);
		int pressedKey = cv::waitKey(0);
		cout << "Pressed key: " << pressedKey << endl;
		if( pressedKey == 1048603 or pressedKey == 27 ) break; // Break on ESC key

		refinement++; //increase the level of refinement of the keyframe depthmap
	}
	while( cv::waitKey(500)!=27 )
	{
		// Do nothing
	}
	return 0;
}
