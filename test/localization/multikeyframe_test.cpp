// Testing out the multi-keyframe SLAM

#include "io.h"
#include "ocv.h"
#include "timer.h"

#include <boost/program_options.hpp>

#include "reconstruction/key_map.h"

using namespace std;
namespace po = boost::program_options;

const float imgIntensityScale = 0.5;

// Helper function to load values into variable
template <int N>
array<double, N> readValues(const po::variables_map & vmap, const std::string & key, bool output = false)
{
	const vector<double> values = vmap[key].as< vector<double> >();
	if ( values.size() != N )
	{
		std::cerr << "ERROR : Option " << key << " does not have " << N << " arguments" << endl;
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
		std::cerr << "ERROR : Parameter file not found at " << argv[1] << std::endl;
		return 1;
	}

	po::options_description options("Multi-Keyframe Test options");
	options.add_options() ("help", "Produce help message");
	options.add_options() ("foldername", po::value<string>()->required(), "Folder to search for the dataset, should contain a \"camera1\" and \"camera2\" folder");
	options.add_options() ("skipRows", po::value<string>()->default_value(0), "Number of rows of images to skip in the Image files");
	options.add_options() ("refinement", po::value<int>()->default_value(0), "Refinement level threshold before photometric localization");

	options.add_options() ("Intrinsics.Camera1", po::value< vector<double> >()->required(), "EUCM Intrinsics for camera 1 (6 doubles)");
	options.add_options() ("Intrinsics.Camera2", po::value< vector<double> >()->required(), "EUCM Intrinsics for camera 2 (6 doubles)");
	options.add_options() ("Extrinsics.Stereo", po::value< vector<double> >()->required(), "EUCM Stereo extrinsic calibration (6 doubles)");

	options.add_options() ("StereoSGMParams.u0", po::value<double>()->required(), "Stereo SGM Parameters - u0");
	options.add_options() ("StereoSGMParams.v0", po::value<double>()->required(), "Stereo SGM Parameters - v0");
	options.add_options() ("StereoSGMParams.dispMax", po::value<int>()->required(), "Stereo SGM Parameters - dispMax (int)");
	options.add_options() ("StereoSGMParams.scale", po::value<int>()->required(), "Stereo SGM Parameters - scale (int)");
	options.add_options() ("StereoSGMParams.verbosity", po::value<int>()->default_value(0), "Stereo SGM Parameters - verbosity level (int)");
	options.add_options() ("StereoSGMParams.salientPoints", po::value<bool>()->required(), "Stereo SGM Parameters - use only salient points? (bool)");
	options.add_options() ("StereoSGMParams.hypMax", po::value<int>()->required(), "Stereo SGM Parameters - number of hypotheses layers (int)");

	options.add_options() ("MotionStereoParams.descLength", po::value<int>()->required(), "Motion Stereo descriptor length in pixels (needs to be an odd int)");
	options.add_options() ("MotionStereoParams.dispMax", po::value<int>()->required(), "Motion Stereo maximum disparity (int)");
	options.add_options() ("MotionStereoParams.scale", po::value<int>()->required(), "Motion Stereo scale (int)");
	options.add_options() ("MotionStereoParams.verbosity", po::value<int>()->default_value(0), "Motion Stereo verbosity level (int)");

	options.add_options() ("Pose.Base2Camera", po::value< vector<double> >()->required(), "6dof pose between robot base and camera 1 (6 doubles)");

	po::variables_map vmap;
	po::store( po::parse_config_file(paramFile, options, true), vmap);

	// Print the help message and quit
	if( vmap.count("help") )
	{
		cout << options << endl;
		return 1;
	}

	po::notify( vmap ); // Checks if all required values have been defined

	// Accept the filename for images for left and right cameras
	const std::string datasetFoldername = vmap["foldername"].as<string>();
	const std::string img1Filename = datasetFoldername + "/leftImageGroundTruth.txt";
	const std::string img2Filename = datasetFoldername + "/rightImageGroundTruth.txt";
	std::ifstream img1File(img1Filename);
	std::ifstream img2File(img2Filename);
	if( not( img1File.is_open() and img2File.is_open() ) )
	{
		std::cerr << "ERROR : Unable to open image list files" << std::endl;
		std::exit(1);
	}
	std::string foldername1, foldername2; // Folders where the images are stored
	getline(img1File, foldername1);
	getline(img2File, foldername2);
	foldername1 = datasetFoldername + "/" + foldername1;
	foldername2 = datasetFoldername + "/" + foldername2;
	std::cout << "Folder for camera1 : " << foldername1 << std::endl;
	std::cout << "Folder for camera2 : " << foldername2 << std::endl;

	// Accept the EUCM params for intrinsic and extrinsic calibratin of the cameras
	array<double, 6> cam1Params = readValues<6>(vmap, "Intrinsics.Camera1", true);
	array<double, 6> cam2Params = readValues<6>(vmap, "Intrinsics.Camera2", true);
	array<double, 6> stereoPose = readValues<6>(vmap, "Extrinsics.Stereo", true);
	EnhancedCamera camera1( cam1Params.data() );
	EnhancedCamera camera2( cam2Params.data() );
	Transformation<double> T12( stereoPose.data() );

	// Skip the first few images of the dataset
	int skiprows = vmap["skipRows"].as<int>();
	for( int i = 1; i < skiprows; ++i )
	{
		std::string dump;
		getline(img1File, dump);
		getline(img2File, dump);
	}

	// Initialize empty keyframe

	// start loop
		// If getkeyframe state set, collect next keyframe
		// Calculate 
		// ???
		// Profit

	return 0;
}