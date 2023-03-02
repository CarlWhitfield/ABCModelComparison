#include"read_mbw_data.h"
#include"mbw_processing_params.h"
#include<file_manip.h>
#include<boost/filesystem.hpp>

using namespace std;

int main(int argc, char *argv[])
{
	string filepath, subject;
	vector<std::shared_ptr<MBWRawFile>> mbw_files;
	if(argc > 1)  //extract subject name first (folder name above current folder)
	{
		boost::filesystem::path path_ss(argv[1]);
		boost::filesystem::path full_path = canonical(path_ss);
		filepath = full_path.parent_path().parent_path().string();
		subject = full_path.parent_path().parent_path().filename().string();
	}
	else
	{
		cout << "Error: not enough arguments!" << endl;
		return 1;
	}

	std::shared_ptr<LCIOptionList> options = std::make_shared<LCIOptionList>();         //start by assigning default option and paramter lists
	std::shared_ptr<LCIParameterList> params = get_raw_mbw_data(argc, argv, mbw_files, options.get());

	std::shared_ptr<MBWData> MBWprocessed = std::make_shared<MBWData>(mbw_files, 
		options.get(), params.get());
	
	MBWprocessed->subject_name = subject;
	//these can all be moved to ABC input
	MBWprocessed->machine_ds = 0.001*params->get_param<double>(MACHINE_DS_PARAM_NAME)->get_value();   //convert to L
	MBWprocessed->extra_rebreathe_vol = 0.001*params->get_param<double>(REBREATHE_VOL_PARAM_NAME)->get_value();  //convert to L
	MBWprocessed->Vbag = params->get_param<double>(MRI_VBAG_PARAM_NAME)->get_value();  //L

	write_processed_washout_data(filepath,MBWprocessed.get());
	write_mbw_summary(filepath,MBWprocessed.get());

	return 0;
}