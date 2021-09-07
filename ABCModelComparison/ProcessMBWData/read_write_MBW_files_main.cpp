#include"read_mbw_data.h"
#include<file_manip.h>
#include<boost/filesystem.hpp>

using namespace std;

int main(int argc, char *argv[])
{
	std::vector<std::string>argv_copy(argc);
	for (int i=0;i<argc;i++) 
	{
		argv_copy[i]=argv[i];
	}
	double Vbag, machine_ds;
	string filepath, subject;
	vector<std::shared_ptr<MBWRawFile>> mbw_files;
	if(argc > 2)
	{
		boost::filesystem::path path_ss(argv[2]);
		filepath = path_ss.parent_path().string();
		subject = path_ss.parent_path().filename().string();
	}
	else
	{
		cout << "Error: not enough arguments!" << endl;
		return 1;
	}

	if(get_raw_mbw_data(argv_copy, mbw_files, Vbag, machine_ds)) return 1;

	std::shared_ptr<MBWData> MBWprocessed = std::make_shared<MBWData>(mbw_files);
	MBWprocessed->subject_name = subject;
	MBWprocessed->Vbag = Vbag;
	MBWprocessed->machine_ds = machine_ds;
	
	write_processed_washout_data(filepath,MBWprocessed.get());
	write_mbw_summary(filepath,MBWprocessed.get());

	return 0;
}