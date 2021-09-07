#include"read_mbw_data.h"
#include"read_write_codes.h"
#include"read_write_settings.h"
#include<file_manip.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<globals.h>
#include<unordered_map>
#include<boost/filesystem.hpp>

using namespace std;

void linear_fit(const Eigen::VectorXd & x, const Eigen::VectorXd & y, double & m, double & c)
{
	if(x.size() != y.size())
	{
		cerr << "Error! Cannot perform linear fit as x and y vectors have different sizes\n";
	}
	else
	{
		int n = x.size();
		Eigen::VectorXd x2m = x.array()*x.array()/double(n), xym = x.array()*y.array()/double(n);
		Eigen::VectorXd xm = x/double(n), ym = y/double(n);
		double x2mean = x2m.sum();
		double xymean = xym.sum();
		double xmean = xm.sum(), ymean = ym.sum();
		m = (xymean - xmean*ymean)/(x2mean - xmean*xmean);
		c = (ymean*x2mean - xmean*xymean)/(x2mean - xmean*xmean);
	}
}

void MBWRawFile::read_raw_data(const std::string & fname)
{
	if(!check_infile(fname))
	{
		ifstream infile;
		string s;
		stringstream ss;
		//flags indicating what to do next
		bool read_conc = false;
		unsigned nlines = count_lines(fname);
		infile.open(fname);
		unsigned count = 0;
		this->conc_data.reserve(nlines);
		this->flow_data.reserve(nlines);
		this->CO2_data.reserve(nlines);
		this->O2_data.reserve(nlines);
		this->time.reserve(nlines);
		while(getline(infile, s))
		{
			ss.str(s);
			string nc;
			getline(ss, nc, '\t');
			if(read_conc)
			{
				this->time.push_back(StringToNumber<double>(nc));   //first column is time
				int n = 1;
				while(getline(ss, nc, '\t'))
				{
					switch(n)    //assign data by column
					{
					case 1:
						{
							this->flow_data.push_back(StringToNumber<double>(nc));
						} break;

					case 2:
						{
							this->conc_data.push_back(StringToNumber<double>(nc));
						} break;

					case 3:
						{
							this->O2_data.push_back(StringToNumber<double>(nc));
						} break;

					case 4:
						{
							this->CO2_data.push_back(StringToNumber<double>(nc));
						} break;

					default:
						break;
					}
					n++;
				}
			}

			if(nc == "Time[ms]") 
			{
				read_conc = true;  //once this appears, start reading conc
			}
			ss.clear();
			count++;
		}
		infile.close();
		this->trim_and_apply_delay();
		//define breaths (in and out)
		this->measure_breath_start_end_points();

		//nothing stored in washin vectors here yet, define them here
		this->define_washin_and_washout();

		//print washout data
		//{
		//	vector<vector<double>> test;
		//	vector<string> test_head;
		//	test_head.push_back("Time");
		//	test.push_back(this->time);
		//	test_head.push_back("Flow");
		//	test.push_back(this->flow_data);
		//	test_head.push_back("SF6");
		//	test.push_back(this->conc_data);
		//	test_head.push_back("CO2");
		//	test.push_back(this->CO2_data);
		//	test_head.push_back("O2");
		//	test.push_back(this->O2_data);
		//	write_csv_file_cols("washout.csv",test_head,test);
		//}

		//print washin data
		//{
		//	vector<vector<double>> test;
		//	vector<string> test_head;
		//	test_head.push_back("Time");
		//	test.push_back(this->time_washin);
		//	test_head.push_back("Flow");
		//	test.push_back(this->flow_data_washin);
		//	test_head.push_back("SF6");
		//	test.push_back(this->conc_data_washin);
		//	test_head.push_back("CO2");
		//	test.push_back(this->CO2_data_washin);
		//	test_head.push_back("O2");
		//	test.push_back(this->O2_data_washin);
		//	write_csv_file_cols("washin.csv",test_head,test);
		//}
	}
	else
	{
		abort_on_failure();
	}
}

//tested and works
void MBWRawFile::trim_and_apply_delay()  //unfinished -- need to think about raw data format
{
	//file in 10 ms steps
	vector<double> time_new, SF6_new, O2_new, CO2_new, flow_new;

	//delay in time steps (converted from secs)
	unsigned CO2d_i = ((unsigned) (100*this->CO2_delay));
	unsigned O2d_i = ((unsigned) (100*this->O2_delay));

	//first identify regions withh 100 flow, data is between these
	size_t max_valid_length = 0, current_valid_length = 0,
		   start_point = 0, current_start_point = 0;
	for(size_t n = 0; n < this->time.size(); n++)   //find bounds of signal
	{
		if(this->flow_data[n] > 99) //means no flow signal
		{
			if(current_valid_length > max_valid_length) //check if this is longest run of valid points
			{
				//if so, store it
				max_valid_length = current_valid_length;   
				start_point = current_start_point;
			}
			current_valid_length = 0;   //end run of valid points
		}
		else  //valid point
		{
			if(current_valid_length == 0)  //if start of run, mark it
			{
				current_start_point = n;
			}
			current_valid_length++;
		}
	}
	if(current_valid_length > max_valid_length)   //check if last run was longest
	{
		max_valid_length = current_valid_length;
		start_point = current_start_point;
	}




	//trim down data and adjust for delays
	//last_invalid_point = std::max(std::max(CO2d_i, O2d_i), last_invalid_point);
	//maximum size is max_valid length
	size_t vect_size = min(max_valid_length - CO2d_i, max_valid_length - O2d_i);
	unsigned end_cutoff = std::max(CO2d_i, O2d_i);
	this->time = vector<double>(this->time.begin() + start_point, 
		                        this->time.begin() + start_point + vect_size);
	//take flow measurement as same as time
	this->flow_data = vector<double>(this->flow_data.begin() + start_point,
		                             this->flow_data.begin() + start_point + vect_size);
	//correct SF6 and C02 for delay
	this->conc_data= vector<double>(this->conc_data.begin() + start_point + CO2d_i,
		                            this->conc_data.begin() + start_point + CO2d_i + vect_size);
	this->CO2_data = vector<double>(this->CO2_data.begin() + start_point + CO2d_i,
		                            this->CO2_data.begin() + start_point + CO2d_i + vect_size);
	this->O2_data = vector<double>(this->O2_data.begin() + start_point + O2d_i, 
		                           this->O2_data.begin() + start_point + O2d_i + vect_size);
}

void MBWRawFile::measure_breath_start_end_points()
{
	bool inhaling;
	this->inhalation_start_pts.reserve(size_t(this->time.size()/(500)));
	this->exhalation_start_pts.reserve(size_t(this->time.size()/(500)));
	//determine if inhale or exhale for first 0.5s
	double sum = 0;
	for(size_t np = 0; np < std::min(size_t(50), this->time.size()); np++)
	{
		sum += this->flow_data[np];
	}
	if(sum >= 0)   //inhalation first
	{
		inhaling = true;
		this->inhalation_start_pts.push_back(0);
	}
	else    //exhalation first
	{
		inhaling = false;
		this->exhalation_start_pts.push_back(0);
	}
	for(int n = 0; n < int(this->time.size()); n++) //loop over whole washout
	{
		if(inhaling)
		{
			if(this->flow_data[n] < 0)   //exhalation
			{
				double exhaled_pts = 1;
				int np = 1;
				while(np < 20 && n+np < int(this->time.size())) //if next .2 of a second mostly exhaling
				{
					if( this->flow_data[n+np] < 0 ) exhaled_pts += 1;
					np++;
				}
				if(exhaled_pts > 15)   //need 75% for change of direction
				{
					inhaling = false;
					this->exhalation_start_pts.push_back(n);
				}
			}
		}
		else
		{
			if(this->flow_data[n] > 0)   //inhalation
			{
				int inhaled_pts = 1;
				int np = 1;
				while(np < 20 && n+np < int(this->time.size())) //if next .2 of a second mostly inhaling
				{
					if(this->flow_data[n+np] > 0) inhaled_pts += 1;
					np++;
				}
				if(inhaled_pts > 15)
				{
					inhaling = true;
					this->inhalation_start_pts.push_back(n);
				}
			}
		}
	}
	//check breaths have the correct sign
	//find exh offset
	int exh_offset = 0;
	if(this->inhalation_start_pts[0] > this->exhalation_start_pts[0])
	{
		exh_offset = 1;
	}

	vector<size_t> wrong_inhalations;
	for(size_t nb = 0; nb < this->inhalation_start_pts.size(); nb++)
	{
		size_t nb_end = nb + exh_offset;
		size_t ib = this->inhalation_start_pts[nb];
		size_t ib_end = 0;
		if(nb_end >= this->exhalation_start_pts.size()) ib_end = this->time.size();  //no exh to follow
		else ib_end = this->exhalation_start_pts[nb_end];
		double sum = 0;
		for(size_t i = ib; i < ib_end; i++)  //loop over current inhalation
		{
			sum += this->flow_data[i];
		}		
		if(sum < 0) wrong_inhalations.push_back(nb);
	}
	//remove incorrect inhalations and following exhalation
	for(int iw = wrong_inhalations.size()-1; iw >= 0 ; iw--) //go backwards
	{
		inhalation_start_pts.erase(inhalation_start_pts.begin() + wrong_inhalations[iw]);
		int exh_index = wrong_inhalations[iw] + exh_offset;  //check following exhalation exists
		if(exh_index >= 0 && exh_index < int(exhalation_start_pts.size()))
		{
			exhalation_start_pts.erase(exhalation_start_pts.begin() + wrong_inhalations[iw] + exh_offset);
		}
	}
	vector<size_t> wrong_exhalations;
	for(size_t nb = 0; nb < this->exhalation_start_pts.size(); nb++)
	{
		size_t nb_end = nb + 1 - exh_offset;
		size_t ib = this->exhalation_start_pts[nb];
		size_t ib_end = 0;
		if(nb_end >= this->inhalation_start_pts.size()) ib_end = this->time.size();  //no inh to follow
		else ib_end = this->inhalation_start_pts[nb_end];
		double sum = 0;
		for(size_t i = ib; i < ib_end; i++)  //loop over current inhalation
		{
			sum += this->flow_data[i];
		}		
		if(sum > 0) wrong_exhalations.push_back(nb);
	}
	for(int iw = wrong_exhalations.size()-1; iw >= 0 ; iw--) //go backwards
	{
		exhalation_start_pts.erase(exhalation_start_pts.begin() + wrong_exhalations[iw]);
		int inh_index =  wrong_exhalations[iw] + 1 - exh_offset; 
		//check following inhalation exists
		if(inh_index >= 0 && inh_index < int(inhalation_start_pts.size()))
		{
			inhalation_start_pts.erase(inhalation_start_pts.begin() + wrong_exhalations[iw]
		                            + 1 - exh_offset);
		}
	}
}

void MBWRawFile::define_washin_and_washout()
{
	bool go = true;
	size_t nb = 0;
	size_t washin_start = 0, washin_breath_start = 0, washin_exh_start = 0;
	//find washin start point
	int exh_offset = 0;
	if(this->inhalation_start_pts[0] > this->exhalation_start_pts[0]) 
	//if starts with an exhalation add a dummy inhalation
	{
		exh_offset = 1;
	}
	while(go && nb < this->inhalation_start_pts.size())
	{
		size_t ib = this->inhalation_start_pts[nb];   //inhalation start point
		size_t ib_end = 0, nb_end = nb+exh_offset;   //get next exh start point
		if(nb_end >= this->exhalation_start_pts.size()) ib_end = this->time.size();  //no exh to follow
		else ib_end = this->exhalation_start_pts[nb_end];
		size_t i = ib;
		while(go && i < ib_end)
		{
			if(this->conc_data[i] > washin_min_conc)  //if inhaled conc exceeds 0.15% there must be a source
			{
				go = false;  //stop searching
				washin_start = ib;   //record start of inhalation as first point of washin
				washin_breath_start = nb;   //record breath number too
				washin_exh_start = nb_end;
			}
			i++;
		}
		nb++;
	}
	
	//work backwards from end to avoid false starts
	size_t washout_start = 0, washout_breath_start = 0, washout_exh_start = 0;
	nb = this->inhalation_start_pts.size()-1;
	go = true;
	//find washout start point
	while(go && nb > washin_breath_start)
	{
		size_t ib = this->inhalation_start_pts[nb];   //inhalation start point
		size_t ib_end = 0, nb_end = nb + exh_offset;   //get next exh start point
		if(nb_end >= this->exhalation_start_pts.size()) ib_end = this->time.size();  //no exh to follow
		else ib_end = this->exhalation_start_pts[nb_end];	
		double inhaled_conc_av = 0;
		for(size_t i = ib; i < ib_end; i++)
		{
			inhaled_conc_av += this->conc_data[i];
		}
		if(ib_end > ib) inhaled_conc_av /= (ib_end - ib);

		if(inhaled_conc_av > washin_min_conc)  //washin inhalation
		{
			go = false;
			washout_breath_start = nb+1;   //record brath number too
			washout_exh_start = nb_end+1;
			washout_start = this->inhalation_start_pts[washout_breath_start];   //record start of inhalation as first point of washin
		}
		nb--;
	}

	//now need to trim down to store washin and washout separately
	//washin data
	this->time_washin = vector<double>(this->time.begin() + washin_start, 
		                          this->time.begin() + washout_start);
	this->flow_data_washin = vector<double>(this->flow_data.begin() + washin_start, 
		                          this->flow_data.begin() + washout_start);
	this->conc_data_washin = vector<double>(this->conc_data.begin() + washin_start, 
		                          this->conc_data.begin() + washout_start);
	this->CO2_data_washin = vector<double>(this->CO2_data.begin() + washin_start, 
		                          this->CO2_data.begin() + washout_start);
	this->O2_data_washin = vector<double>(this->O2_data.begin() + washin_start, 
		                          this->O2_data.begin() + washout_start);

	//washin breath data
	this->inhalation_start_pts_washin = vector<size_t>(this->inhalation_start_pts.begin() + washin_breath_start,
		                                               this->inhalation_start_pts.begin() + washout_breath_start);
	this->exhalation_start_pts_washin = vector<size_t>(this->exhalation_start_pts.begin() + washin_exh_start,
		                                               this->exhalation_start_pts.begin() + washout_exh_start);
	for(size_t nb = 0; nb < this->inhalation_start_pts_washin.size(); nb++)  //correct start points
	{
		this->inhalation_start_pts_washin[nb] -= washin_start;
	}
	for(size_t nb = 0; nb < this->exhalation_start_pts_washin.size(); nb++)  //correct start points
	{
		this->exhalation_start_pts_washin[nb] -= washin_start;
	}

	//washout data
	this->time = vector<double>(this->time.begin() + washout_start, this->time.end());
	this->flow_data = vector<double>(this->flow_data.begin() + washout_start, this->flow_data.end());
	this->conc_data = vector<double>(this->conc_data.begin() + washout_start, this->conc_data.end());
	this->CO2_data = vector<double>(this->CO2_data.begin() + washout_start, this->CO2_data.end());
	this->O2_data = vector<double>(this->O2_data.begin() + washout_start, this->O2_data.end());

	//washout breath data
	this->inhalation_start_pts = vector<size_t>(this->inhalation_start_pts.begin() + washout_breath_start,
		                                        this->inhalation_start_pts.end());
	this->exhalation_start_pts = vector<size_t>(this->exhalation_start_pts.begin() + washout_exh_start,
		                                        this->exhalation_start_pts.end());
	for(size_t nb = 0; nb < this->inhalation_start_pts.size(); nb++)  //correct start points
	{
		this->inhalation_start_pts[nb] -= washout_start;
	}
	for(size_t nb = 0; nb < this->exhalation_start_pts.size(); nb++)  //correct start points
	{
		this->exhalation_start_pts[nb] -= washout_start;
	}
}

int get_raw_mbw_data(const std::vector<std::string> & argv, 
					 std::vector<std::shared_ptr<MBWRawFile>> & mbw_files,
					 double & Vbag, double & machine_ds)
{
	//command line inputs
	//nfiles, file, file, ..., file,  CO2delay, VMRIbag, machine_ds, O2delay
	int argc = int(argv.size());
	if(argc > 2)
	{
		//fist argument is number of tests to do
		int Nfiles = StringToNumber<int>(argv[1]);
		mbw_files.resize(Nfiles);
		if(argc > 1 + Nfiles)
		{
			double CO2_delay, O2_delay;
			//C02 delay will be first argument after filenames
			if(argc > 2 + Nfiles)
			{
				CO2_delay = StringToNumber<double>(argv[2+Nfiles]);
			}
			else 
			{
				cout << "Please give Flow Gas Delay (in seconds):";
				cin >> CO2_delay;
			}
			//O2 delay will be same (it isn't, but not used currently anyway)
			
			if(argc > 3 + Nfiles) 
			{
				Vbag = StringToNumber<double>(argv[3+Nfiles]);
			}
			else 
			{
				cout << "Enter MRI bag volume used (in L):";
				cin >> Vbag;
			}

			if(argc > 4 + Nfiles)
			{
				machine_ds = StringToNumber<double>(argv[4+Nfiles]);
			}
			else 
			{
				cout << "Enter MBW machine deadspace used (in L):";
				cin >> machine_ds;
			}

			if(argc > 5 + Nfiles)
			{
				O2_delay = StringToNumber<double>(argv[5+Nfiles]);
			}
			else 
			{
				cout << "Assuming O2 delay same as CO2 delay." << endl;
				O2_delay = CO2_delay;
			}

			for(int n = 0; n < Nfiles; n++)
			{
				mbw_files[n] = std::make_shared<MBWRawFile>(argv[2+n], CO2_delay, O2_delay);
			}
		}
		else
		{
			cout << "Error: not enough arguments!" << endl;
			return 1;
		}
	}
	else
	{
		cout << "Error: not enough arguments!" << endl;
		return 1;
	}

	return 0;
}

void MBWTest::process_MBW_inputs(std::shared_ptr<MBWRawFile> MBWFile)
{
	//adjust fluxes for BTPS
	Eigen::VectorXd adjusted_washin_flux, adjusted_washout_flux;
	this->adjust_flow_for_BTPS_and_bias(MBWFile, adjusted_washin_flux, adjusted_washout_flux);

	//first get end washin conc
	this->measure_Cinit(MBWFile, adjusted_washin_flux);

	//measure washout breath volumes
	this->measure_washout_breath_volumes(MBWFile, adjusted_washout_flux);

	this->measure_LCI_FRC_FDS(MBWFile, adjusted_washout_flux);
	
	//trim dataset down
	size_t BreathsToKeep;
	if(breath_cutoff_option == 'l')
	{
		BreathsToKeep = min(this->Nbreaths, this->LCI_point+3);
	}
	else
	{
		BreathsToKeep = min(this->Nbreaths, Nbreaths_keep);
	}
	this->FRC0 = (this->cumulative_exhaled_SF6vols[BreathsToKeep-1] 
	             - this->cumulative_inhaled_SF6vols[BreathsToKeep-1])
	           / (this->Cinit - this->Cet[BreathsToKeep-1]);
	this->LCI = this->cumulative_exhaled_vols[this->LCI_point]/this->FRC0;

	this->cumulative_exhaled_vols = vector<double>(this->cumulative_exhaled_vols.begin(),
		                            this->cumulative_exhaled_vols.begin() + BreathsToKeep);
	this->cumulative_inhaled_vols = vector<double>(this->cumulative_inhaled_vols.begin(),
		                            this->cumulative_inhaled_vols.begin() + BreathsToKeep);
	this->cumulative_exhaled_SF6vols = vector<double>(this->cumulative_exhaled_SF6vols.begin(),
		                            this->cumulative_exhaled_SF6vols.begin() + BreathsToKeep);
	this->cumulative_inhaled_SF6vols = vector<double>(this->cumulative_inhaled_SF6vols.begin(),
		                            this->cumulative_inhaled_SF6vols.begin() + BreathsToKeep);
	this->exhaled_breath_vols = vector<double>(this->exhaled_breath_vols.begin(),
		                            this->exhaled_breath_vols.begin() + BreathsToKeep);
	this->inhaled_breath_vols = vector<double>(this->inhaled_breath_vols.begin(),
		                            this->inhaled_breath_vols.begin() + BreathsToKeep);
	this->exhaled_SF6vols = vector<double>(this->exhaled_SF6vols.begin(),
		                            this->exhaled_SF6vols.begin() + BreathsToKeep);
	this->inhaled_SF6vols = vector<double>(this->inhaled_breath_vols.begin(),
		                            this->inhaled_breath_vols.begin() + BreathsToKeep);
	this->inhaled_breath_pts = vector<size_t>(this->inhaled_breath_pts.begin(),
		                       this->inhaled_breath_pts.begin() + BreathsToKeep);
	this->exhaled_breath_pts = vector<size_t>(this->exhaled_breath_pts.begin(),
		                       this->exhaled_breath_pts.begin() + BreathsToKeep);
	this->FDS = vector<double>(this->FDS.begin(), this->FDS.begin() + BreathsToKeep);
	this->Cet = vector<double>(this->Cet.begin(), this->Cet.begin() + BreathsToKeep);
	this->CO2et = vector<double>(this->CO2et.begin(), this->CO2et.begin() + BreathsToKeep);
	double median_VT = median<double>(this->exhaled_breath_vols);
	if(median<double>(this->CO2et) == 0 && BreathsToKeep > 3)
	{
		vector<double> FDSh(this->FDS);
		FDSh.erase(FDSh.begin()+3, FDSh.end());
		this->median_FDS = median<double>(FDSh);
	}
	else
	{
		this->median_FDS = median<double>(this->FDS);
	}
	this->Nbreaths = BreathsToKeep;  //update Nbreaths
	//do exhalation measurements
	this->define_breath_measurements(MBWFile, adjusted_washout_flux);
}

void MBWTest::adjust_flow_for_BTPS_and_bias(std::shared_ptr<MBWRawFile> MBWfile, 
			                       Eigen::VectorXd & adjusted_washin_flow, 
						           Eigen::VectorXd & adjusted_washout_flow)
{
	Eigen::VectorXd volume, v_in;
	size_t total_timesteps = MBWfile->count_washin_timesteps() + MBWfile->count_washout_timesteps();
	size_t Nwashin = MBWfile->count_washin_timesteps();
	size_t Nwashout = MBWfile->count_washout_timesteps();
	adjusted_washin_flow = Eigen::VectorXd::Zero(Nwashin);
	adjusted_washout_flow = Eigen::VectorXd::Zero(Nwashout);

	//define window for counting up tidal volume in and out (trim off first and last breaths)
	double cumul_exh = 0;
	double cumul_inh = 0;
	size_t start_count, end_count;
	if(MBWfile->count_washin_inhalations() > 1) start_count = MBWfile->get_washin_inhalation_start(1);
	else start_count = MBWfile->get_washout_inhalation_start(0);
	size_t Nbwashout = MBWfile->count_washout_inhalations();
	if(Nbwashout > 1) end_count = MBWfile->get_washout_inhalation_start(Nbwashout-2);
	else end_count = MBWfile->get_washout_inhalation_start(0);
	//volume = Eigen::VectorXd::Zero(total_timesteps + 1);
	//v_in = Eigen::VectorXd::Zero(total_timesteps + 1);
	/*double m=1,c,m_last=1;
	int iteration = 0;
	while(abs(m_last) > 0.001)
	{*/
	for(size_t i = 0; i < total_timesteps; i++)
	{
		double adjusted_flow_h;
		if(MBWfile->get_flow(i) < 0) //expiration
		{
			if(i < Nwashin)
			{
				/*if(iteration==0)*/ adjusted_washin_flow[i] = MBWfile->get_flow(i)*BTPSout;
				adjusted_flow_h = adjusted_washin_flow[i];
			}
			else
			{
				/*if(iteration==0)*/ adjusted_washout_flow[i - Nwashin] = MBWfile->get_flow(i)*BTPSout;
				adjusted_flow_h = adjusted_washout_flow[i - Nwashin];
			}
			if(i >= start_count && i < end_count) cumul_exh += -adjusted_flow_h*0.01;
			/*v_in[i+1] = v_in[i];*/
		}
		else  //inspiration
		{
			if(i < Nwashin)
			{
				/*if(iteration==0)*/ adjusted_washin_flow[i] = MBWfile->get_flow(i);
				/*else adjusted_washin_flow[i] *= (1-m);*/
				adjusted_flow_h = adjusted_washin_flow[i];
				
			}
			else
			{
				/*if(iteration==0)*/ adjusted_washout_flow[i - Nwashin] = MBWfile->get_flow(i);
				/*else adjusted_washout_flow[i - Nwashin] *= (1-m);*/
				adjusted_flow_h = adjusted_washout_flow[i - Nwashin];
			}
			/*v_in[i+1] = v_in[i] + adjusted_flow_h*0.01;*/
			if(i >= start_count && i < end_count) cumul_inh += adjusted_flow_h*0.01;
		}
		/*volume[i+1] = volume[i] + adjusted_flow_h*0.01;*/
	}
	//m_last=m;
	//linear_fit(v_in, volume, m, c);
	//iteration++;
	//}

	//normalise inhalation volumes by exhalation
	double inh_sf = cumul_exh/cumul_inh;
	for(size_t i = 0; i < total_timesteps; i++)
	{
		if(MBWfile->get_flow(i) > 0) //inspiration
		{
			if(i < Nwashin) adjusted_washin_flow[i] *= inh_sf;
			else adjusted_washout_flow[i - Nwashin] *= inh_sf;
		}
	}
}

void MBWTest::measure_Cinit(std::shared_ptr<MBWRawFile> MBWFile, 
							const Eigen::VectorXd & washin_flux)
{
	size_t last_washin_exh = MBWFile->count_washin_exhalations() - 1;
	double cumulative_vol = 0;
	for(size_t i = MBWFile->get_washin_exhalation_start(last_washin_exh);
		       i < MBWFile->count_washin_timesteps(); i++)   //loop over last exhalation of washout
	{
		cumulative_vol -= washin_flux[i]*0.01;
	}
	//measure Cinit
	double m1 = 0, m2 = 0;
	size_t Cinit_pts = 0;
	double sum_exhaled = 0;
	for(size_t i = MBWFile->get_washin_exhalation_start(last_washin_exh);
		       i < MBWFile->count_washin_timesteps(); i++)   //loop over last exhalation of washout
	{
		sum_exhaled -= washin_flux[i]*0.01;
		if(sum_exhaled >= 0.6*cumulative_vol && sum_exhaled < 0.9*cumulative_vol)  //measure only in these limits
		{
			m1 += MBWFile->get_washin_conc(i);
			m2 += MBWFile->get_washin_conc(i)*MBWFile->get_washin_conc(i);
			Cinit_pts++;
		}
	}
	m1 /= double(Cinit_pts);
	m2 /= double(Cinit_pts);
	this->Cinit = m1;
	double var = m2 - m1*m1;  //machine error can cause issues if var very small
	this->Cinit_std = sqrt(max(var,0.0));
}

void MBWTest::measure_washout_breath_volumes(std::shared_ptr<MBWRawFile> MBWFile, 
											 const Eigen::VectorXd & washout_flux)
{
	this->Nbreaths = min(MBWFile->count_washout_inhalations(), MBWFile->count_washout_exhalations());
	this->inhaled_breath_vols.resize(this->Nbreaths);
	this->cumulative_inhaled_vols.resize(this->Nbreaths);
	this->inhaled_SF6vols.resize(this->Nbreaths);
	this->cumulative_inhaled_SF6vols.resize(this->Nbreaths);
	this->inhaled_breath_pts.resize(this->Nbreaths);

	this->exhaled_breath_vols.resize(this->Nbreaths);
	this->exhaled_breath_pts.resize(this->Nbreaths);
	this->cumulative_exhaled_vols.resize(this->Nbreaths);
	this->exhaled_SF6vols.resize(this->Nbreaths);
	this->cumulative_exhaled_SF6vols.resize(this->Nbreaths);

	//store all relevant inhalation and exhalation volumes in sequence
	for(int nb = 0; nb < this->Nbreaths; nb++)  //loop over breaths
	{
		size_t i = MBWFile->get_washout_inhalation_start(nb);   //start at beginning of inhalation
		this->inhaled_breath_vols[nb] = 0;  //initialise
		this->inhaled_SF6vols[nb] = 0;
		if(size_t(nb) < MBWFile->count_washout_exhalations())   //if an exhalation is to follow
		{
			while(i < MBWFile->get_washout_exhalation_start(nb))   //keep going until next exhalation
			{
				this->inhaled_breath_vols[nb] += washout_flux[i]*0.01; //count up breath vol
				this->inhaled_SF6vols[nb] += washout_flux[i]*0.01*MBWFile->get_washout_conc(i);
				i++;
			}
			//count number of pts in breath
			this->inhaled_breath_pts[nb] = MBWFile->get_washout_exhalation_start(nb) 
			                             - MBWFile->get_washout_inhalation_start(nb);
		}
		if(nb > 0)
		{
			this->cumulative_inhaled_vols[nb] = this->cumulative_inhaled_vols[nb-1]
			                                 + this->inhaled_breath_vols[nb];
			this->cumulative_inhaled_SF6vols[nb] = this->cumulative_inhaled_SF6vols[nb-1] 
		                                          + this->inhaled_SF6vols[nb];
		}
		else
		{
			this->cumulative_inhaled_vols[nb] = this->inhaled_breath_vols[nb];
			this->cumulative_inhaled_SF6vols[nb] = this->inhaled_SF6vols[nb];
		}

		i = MBWFile->get_washout_exhalation_start(nb); //start at beginning of exhalation
		this->exhaled_breath_vols[nb] = 0;  //initialise
		this->exhaled_SF6vols[nb] = 0;
		if(size_t(nb+1) < MBWFile->count_washout_inhalations())   //if an inhalation is to follow
		{
			while(i < MBWFile->get_washout_inhalation_start(nb+1))    //keep going until next inhalation
			{
				this->exhaled_breath_vols[nb] -= washout_flux[i]*0.01;
				this->exhaled_SF6vols[nb] -= washout_flux[i]*0.01*MBWFile->get_washout_conc(i);
				i++;
			}
			//count number of pts in breath
			this->exhaled_breath_pts[nb] = MBWFile->get_washout_inhalation_start(nb+1)
				                         - MBWFile->get_washout_exhalation_start(nb);
		}
		//record cumulative expired volume
		if(nb > 0)
		{
			this->cumulative_exhaled_vols[nb] = this->cumulative_exhaled_vols[nb-1]
		                                             + this->exhaled_breath_vols[nb];
			this->cumulative_exhaled_SF6vols[nb] = this->cumulative_exhaled_SF6vols[nb-1]
		                                             + this->exhaled_SF6vols[nb];
		}
		else
		{
			this->cumulative_exhaled_vols[nb] = this->exhaled_breath_vols[nb];
			this->cumulative_exhaled_SF6vols[nb] = this->exhaled_SF6vols[nb];
		}
	}
}

void MBWTest::measure_LCI_FRC_FDS(std::shared_ptr<MBWRawFile> MBWFile, 
								  const Eigen::VectorXd & washout_flux)
{
	//approximate LCI and FDS
	this->Cet.resize(this->Nbreaths);   //end tidal conc container
	this->CO2et.resize(this->Nbreaths);   //end tidal conc container
	this->FDS.resize(this->Nbreaths);   //FDS container
	for(int ne = 0; ne < this->Nbreaths; ne++)
	{
		size_t n = MBWFile->get_washout_exhalation_start(ne);
		double cev = 0;
		this->Cet[ne] = 0;
		this->CO2et[ne] = 0;
		int samples = 0;
		while (cev < 0.95*this->exhaled_breath_vols[ne])
		{
			cev -= washout_flux[n]*0.01;
			if(cev >= 0.9*this->exhaled_breath_vols[ne])
			{
				this->Cet[ne] += MBWFile->get_washout_conc(n);
				//do same for CO2
				this->CO2et[ne] += MBWFile->get_washout_CO2(n);
				samples = 0;
			}
			n++;
		}
		if(samples > 0)
		{
			this->Cet[ne] /= double(samples);
			this->CO2et[ne] /= double(samples);
		}
		else 
		{
			this->Cet[ne] = MBWFile->get_washout_conc(n-1);
			this->CO2et[ne] = MBWFile->get_washout_CO2(n-1);
		}

		bool CO2_exists = false;
		if(this->CO2et[ne] > 0) CO2_exists = true;

		vector<double> AUC, AOC;
		AUC.resize(exhaled_breath_pts[ne]+1);
		AOC.resize(exhaled_breath_pts[ne]+1);
		AUC[0] = 0;
		AOC[this->exhaled_breath_pts[ne]] = 0;
		for(size_t nt = 0; nt < this->exhaled_breath_pts[ne]; nt++) //integrate under conc-vol curve
		{
			n = MBWFile->get_washout_exhalation_start(ne) + nt;
			size_t nt_back = this->exhaled_breath_pts[ne] - nt;
			size_t nback = MBWFile->get_washout_exhalation_start(ne) 
				         +  this->exhaled_breath_pts[ne] - 1 - nt;
			double ch = MBWFile->get_washout_conc(n);
			double ceth = this->Cet[ne];
			double cback = MBWFile->get_washout_conc(nback);
			if(CO2_exists)
			{
				ch = MBWFile->get_washout_CO2(n);
				ceth = this->CO2et[ne];
				cback = MBWFile->get_washout_CO2(nback);
			}
			AUC[nt+1] = AUC[nt] - washout_flux[n]*0.01*ch;
			AOC[nt_back - 1] = AOC[nt_back] - washout_flux[nback]*0.01*(ceth - cback);	
		}
		size_t nt = 0;
		double exp_vol = 0, exp_vol_last = 0;
		while(nt < AUC.size() && AUC[nt] < AOC[nt])
		{
			n = MBWFile->get_washout_exhalation_start(ne) + nt;
			size_t nback = MBWFile->get_washout_exhalation_start(ne) + this->exhaled_breath_pts[ne] - 1 - nt;
			exp_vol_last = exp_vol;
			exp_vol += -washout_flux[nback]*0.01;
			nt++;
		}
		if(nt == AUC.size()) this->FDS[ne] = exp_vol;
		else
		{
			if(nt > 0)
			{
				this->FDS[ne] = exp_vol_last + -((AOC[nt - 1] - AUC[nt - 1]) 
							   / (AOC[nt] - AOC[nt - 1] - AUC[nt] + AUC[nt - 1]))
							   * (exp_vol - exp_vol_last);
			}
			else
			{
				this->FDS[ne] = exp_vol_last;
			}
		}
	}

	//cut down test size based on LCI
	size_t last_valid_breath =  this->Nbreaths - 1;
	while(this->exhaled_breath_vols[last_valid_breath] < 2.0*this->median_FDS)
	{
		last_valid_breath--;
	}

	//find LCI
	bool finished = false;
	this->LCI_point = 0;
	while(!finished)
	{
		if(this->exhaled_breath_vols[this->LCI_point] >= 2.0*this->median_FDS) //only check valid points
		{
			if(this->Cet[this->LCI_point] < LCI_frac*this->Cinit)  //check if LCI cutoff passed
			{
				int consec = 1;
				int i = 1;
				while(size_t(this->LCI_point + i) <= last_valid_breath
				   && this->Cet[this->LCI_point + i] < LCI_frac*this->Cinit && consec < 3)
				{
					if(this->exhaled_breath_vols[this->LCI_point + i] >= 2.0*this->median_FDS) // check if valid
					{
						consec++;
					}
					i++;
				}
				if( consec == 3 || size_t(this->LCI_point + i) == this->Cet.size())  finished = true;
			}
		}
		if(size_t(this->LCI_point) < last_valid_breath && !finished) this->LCI_point++;
		else finished = true;
	}
}

void MBWTest::define_breath_measurements(std::shared_ptr<MBWRawFile> MBWFile, 
										 const Eigen::VectorXd & washout_flux)
{
	size_t BreathsToMeasure;
	if(MEASURE_SUBSET) BreathsToMeasure = min(size_t(this->Nbreaths), Nmeasure);
	else  BreathsToMeasure = size_t(this->Nbreaths);
	double NextToMeasure = 0, step = double(this->Nbreaths) / double(BreathsToMeasure);
	double VTmed = median<double>(this->exhaled_breath_vols);   //median breath vol
	int Npts_per_breath = NphaseII + NphaseIII;
	this->conc_measurements.reserve(Npts_per_breath*BreathsToMeasure);
	this->measurement_times.reserve(Npts_per_breath*BreathsToMeasure);
	this->measurement_steps.reserve(Npts_per_breath*BreathsToMeasure);
	double vol_sim_step = 0;
	if(FIXED_STEP_SIZE)  //run simulation with fixed step size on both inhalation and exhalation
	{
		double totbvol = 0;
		for(size_t i = 0; i < this->exhaled_breath_vols.size(); i++)
		{
			totbvol += this->exhaled_breath_vols[i];
		}
		for(size_t i = 0; i < this->inhaled_breath_vols.size(); i++)
		{
			totbvol += this->inhaled_breath_vols[i];
		}
		this->sim_vol_steps.reserve(size_t(2*this->Nbreaths/vol_sim_step_frac)+1);
		this->sim_step_durations.reserve(size_t(2*this->Nbreaths/vol_sim_step_frac)+1);
		vol_sim_step = vol_sim_step_frac*totbvol/(2*this->Nbreaths);
	}
	else   //run simulation with vol steps defined by measurement points
	{
		this->sim_vol_steps.reserve(Npts_per_breath*BreathsToMeasure + 2*this->Nbreaths);
		this->sim_step_durations.reserve(Npts_per_breath*BreathsToMeasure + 2*this->Nbreaths);
	}

	for(int CurrentBreath = 0; CurrentBreath < this->Nbreaths; CurrentBreath++)
	{
		//add inhalation steps to simulation points
		vector<double> cumul_vol_steps;
		if(FIXED_STEP_SIZE)
		{
			size_t vsteps =  size_t(this->inhaled_breath_vols[CurrentBreath]/vol_sim_step) + 1;
			cumul_vol_steps.resize(vsteps);
			for(size_t j = 0; j < vsteps; j++) 
			{
				if(j < vsteps - 1) //all but last step the same
				{
					this->sim_vol_steps.push_back(vol_sim_step);
				}
				else //final step completes breath
				{
					if(j > 0)
					{
						this->sim_vol_steps.push_back(this->inhaled_breath_vols[CurrentBreath]  
			                                       - cumul_vol_steps[j-1]);
					}
					else
					{
						this->sim_vol_steps.push_back(this->inhaled_breath_vols[CurrentBreath] );
					}
				}
				if(j > 0) cumul_vol_steps[j] = cumul_vol_steps[j-1] + this->sim_vol_steps.back();   //cumulative vol
				else cumul_vol_steps[j] = this->sim_vol_steps.back();
			}
			
			size_t j = 0;
			double civ = 0;   //cumulative inspired vol
			double t_last = 0;
			for(size_t it = 0; it < inhaled_breath_pts[CurrentBreath]; it++)
			{
				size_t i = MBWFile->get_washout_inhalation_start(CurrentBreath) + it;
				double civ_old = civ;
				civ += washout_flux[i]*0.01;
				while(j < cumul_vol_steps.size() && civ >= cumul_vol_steps[j])  //measure once cev exceeds volume
				{
					double th;
					if(it > 0) th = ((civ - cumul_vol_steps[j])*(it - 1) + (cumul_vol_steps[j] - civ_old)*it)*0.01 / 
						              (civ - civ_old);
					else th = 0;
					this->sim_step_durations.push_back(th - t_last);
					t_last = th;
					j++;
				}
			}
			if(j == cumul_vol_steps.size() - 1) this->sim_step_durations.push_back(
				                     (inhaled_breath_pts[CurrentBreath] - 1)*0.01 - t_last);

		}
		else
		{
			//otherwise inhalation is a single step
			this->sim_vol_steps.push_back(this->inhaled_breath_vols[CurrentBreath]);
			this->sim_step_durations.push_back(this->inhaled_breath_pts[CurrentBreath]*0.01);
		}

		//exhalation
		cumul_vol_steps.resize(0);
		if(FIXED_STEP_SIZE)
		{
			size_t vsteps =  size_t(this->exhaled_breath_vols[CurrentBreath]/vol_sim_step) + 1;
			cumul_vol_steps.resize(vsteps);
			for(size_t j = 0; j < vsteps; j++) 
			{
				if(j < vsteps - 1) //all but last step the same
				{
					this->sim_vol_steps.push_back(-vol_sim_step);  //-ve for exhalation
				}
				else //final step completes breath
				{
					if(j > 0) 
					{
						this->sim_vol_steps.push_back(cumul_vol_steps[j-1] 
							                - this->exhaled_breath_vols[CurrentBreath]);
					}
					else
					{
						this->sim_vol_steps.push_back(-this->exhaled_breath_vols[CurrentBreath]);
					}
				}
				if(j > 0) cumul_vol_steps[j] = cumul_vol_steps[j-1] - this->sim_vol_steps.back();   //cumulative vol
				else cumul_vol_steps[j] = -this->sim_vol_steps.back();
			}
			size_t j = 0;
			double cev = 0;   //cumulative inspired vol
			double t_last = 0;
			for(size_t it = 0; it < exhaled_breath_pts[CurrentBreath]; it++)
			{
				size_t i = MBWFile->get_washout_exhalation_start(CurrentBreath) + it;
				double cev_old = cev;
				cev -= washout_flux[i]*0.01;
				while(j < cumul_vol_steps.size() && cev >= cumul_vol_steps[j])  //measure once cev exceeds volume
				{
					double th;
					if(it > 0) th = ((cev - cumul_vol_steps[j])*(it - 1) + (cumul_vol_steps[j] - cev_old)*it)*0.01 / 
						              (cev - cev_old);
					else th = 0;
					this->sim_step_durations.push_back(th - t_last);
					t_last = th;
					j++;
				}
			}
			if(j == cumul_vol_steps.size() - 1) this->sim_step_durations.push_back(
				               (exhaled_breath_pts[CurrentBreath] - 1)*0.01 - t_last);
		}


		//get measurements
		if(CurrentBreath == size_t(NextToMeasure))
		{
			if(this->exhaled_breath_vols[CurrentBreath] >= 2.0*this->median_FDS)  //check breath is valid
			{
				//measure_points
				vector<double> measure_vols(Npts_per_breath);
				for(size_t i = 0; i < NphaseII; i++)
				{
					measure_vols[i] = this->median_FDS*(0.5 + (double(i)/double(NphaseII-1)));
				}
				for(size_t i = 0; i < NphaseIII; i++)
				{
					double vol_left = 0.95*this->exhaled_breath_vols[CurrentBreath] - 1.5*this->median_FDS;
					measure_vols[NphaseII+i] = 1.5*this->median_FDS + (double(i+1)/double(NphaseIII))*vol_left;
				}
				if(FIXED_STEP_SIZE)  //shift measured vols to align with sim steps
				{
					size_t i = 0, j=0;
					size_t start_step = this->sim_vol_steps.size() - cumul_vol_steps.size();
					while(j < cumul_vol_steps.size())
					{
						while(i < measure_vols.size() && cumul_vol_steps[j] >= measure_vols[i])
						{
							if(j > 0) measure_vols[i] = cumul_vol_steps[j-1];
							else measure_vols[i] = cumul_vol_steps[j];
							this->measurement_steps.push_back(start_step + j);
							i++;
						}
						j++;
					}
				} 
				else   //measurement steps are the same as sim steps -> add sim steps
				{
					cumul_vol_steps = measure_vols;
					size_t start_step = this->sim_vol_steps.size();
					for(size_t j = 0; j < cumul_vol_steps.size(); j++)
					{
						if(j > 0) this->sim_vol_steps.push_back(cumul_vol_steps[j-1] - cumul_vol_steps[j]); //-ve
						else this->sim_vol_steps.push_back(- cumul_vol_steps[j]);
						this->measurement_steps.push_back(start_step + j);
					}
					this->sim_vol_steps.push_back(cumul_vol_steps.back() - this->exhaled_breath_vols[CurrentBreath]); //-ve: rest of exhalation
				}


				//measure conc values at these steps from data
				size_t j = 0;
				double cev = 0, cev_old = 0;
				double th = 0, t_last = 0;
				//loop over time points from current breath
				for(size_t it = 0; it < exhaled_breath_pts[CurrentBreath]; it++)
				{
					cev_old = cev;
					size_t i = MBWFile->get_washout_exhalation_start(CurrentBreath) + it;
					cev -= washout_flux[i]*0.01;
					//loop until we exceed next vol step
					while(j < measure_vols.size() && cev > measure_vols[j])  //measure once cev exceeds volume
					{
						double conc = 0;
						//interpolate between neighbouring points if there is an i-1 point
						if(it > 0) 
						{
							conc = ((measure_vols[j] - cev_old)*MBWFile->get_washout_conc(i)
							            +  (cev - measure_vols[j])*MBWFile->get_washout_conc(i-1))/ (cev - cev_old);
							th = 0.01*((measure_vols[j] - cev_old)*it + (cev - measure_vols[j])*(it-1))/ (cev - cev_old);
						}
						else conc = MBWFile->get_washout_conc(i);
						//add to conc measurement
						this->conc_measurements.push_back(conc);
						//record vol step where measurement was taken
						this->measurement_times.push_back(i*0.01);
						if(!FIXED_STEP_SIZE) this->sim_step_durations.push_back(th-t_last);
						t_last = th;
						j++;
					}
				}
				double t_end = 0.01*(exhaled_breath_pts[CurrentBreath] - 1);
				this->sim_step_durations.push_back(t_end - t_last);
				NextToMeasure += step;
			}
			else   //breath is not measured
			{
				NextToMeasure += 1;
				if(!FIXED_STEP_SIZE) 
				{
					this->sim_vol_steps.push_back(-this->exhaled_breath_vols[CurrentBreath]);
					this->sim_step_durations.push_back(this->exhaled_breath_pts[CurrentBreath]*0.01);
				}
			}
		}
		else  //breath is not measured
		{
			if(!FIXED_STEP_SIZE) 
			{
				this->sim_vol_steps.push_back(-this->exhaled_breath_vols[CurrentBreath]);
				this->sim_step_durations.push_back(this->exhaled_breath_pts[CurrentBreath]*0.01);
			}
		}
	}
}

void MBWData::process_data_inputs(const std::vector<std::shared_ptr<MBWRawFile>> & MBWfiles)
{
	this->sim_vol_steps.resize(MBWfiles.size());
	this->sim_step_durations.resize(MBWfiles.size());
	this->measurement_times.resize(MBWfiles.size());
	this->measurement_steps.resize(MBWfiles.size());
	this->conc_measurements.resize(MBWfiles.size());
	this->Cinit.resize(MBWfiles.size());
	this->Cinit_std.resize(MBWfiles.size());
	this->LCI.resize(MBWfiles.size());

	this->Tinsp = 0;
	this->Texp = 0;
	this->VT0 = 0;
	this->FRC0 = 0;
	this->dead_space = 0;
	for(size_t n = 0; n < MBWfiles.size(); n++)
	{
		MBWTest MBWtest(MBWfiles[n]);
		this->sim_vol_steps[n] = MBWtest.sim_vol_steps;
		this->sim_step_durations[n] = MBWtest.sim_step_durations;
		this->measurement_times[n] = MBWtest.measurement_times;
		this->measurement_steps[n] = MBWtest.measurement_steps;
		this->conc_measurements[n] = MBWtest.conc_measurements;
		this->Cinit[n] = MBWtest.Cinit;
		this->Cinit_std[n] = MBWtest.Cinit_std;
		this->LCI[n] = MBWtest.LCI;

		this->Tinsp += 0.01*double(median(MBWtest.inhaled_breath_pts));
		this->Texp += 0.01*double(median(MBWtest.exhaled_breath_pts));
		this->VT0 += median(MBWtest.exhaled_breath_vols);
		this->FRC0 += MBWtest.FRC0;
		this->dead_space += MBWtest.median_FDS;
	}
	this->FRC0 /= double(MBWfiles.size());
	this->dead_space /= double(MBWfiles.size());
	this->Tinsp /= double(MBWfiles.size());
	this->Texp /= double(MBWfiles.size());
	this->VT0 /= double(MBWfiles.size());
}

void write_processed_washout_data(const std::string & filepath, const MBWData* mbw_tests)
{
	for(int it = 0; it < int(mbw_tests->sim_vol_steps.size()); it++)
	{
		boost::filesystem::path path(filepath);
		stringstream ss;
		ss << MBWTestFileHead << "_" << mbw_tests->subject_name 
			<< "_" << it << MBWFileExtension;
		string filename = ss.str().c_str();
		path /= filename;

		ofstream outfile;
		outfile.open(path.string());

		outfile << DURATION_NAME << "," << VOLUME_NAME << "," 
			    << MEASURED_NAME << "," << CONC_NAME << ","
				<< DIST_WEIGHT_NAME << endl;
		int im = 0;
		for(int is = 0; is < int(mbw_tests->sim_vol_steps[it].size()); is++)
		{
			outfile << mbw_tests->sim_step_durations[it][is] << "," 
					<< mbw_tests->sim_vol_steps[it][is] << ",";
			if(mbw_tests->measurement_steps[it][im] == is)
			{
				outfile << "1," << mbw_tests->conc_measurements[it][im];
				if(im < mbw_tests->conc_measurements[it].size() - 1) im++;
			}
			else
			{
				outfile << "0,0,0";
			}
			outfile << endl;
		}	
		outfile.close();
	}
}

void write_mbw_summary(const std::string & filepath, const MBWData* mbw_tests)
{
	boost::filesystem::path path(filepath);
	stringstream ss;
	ss << MBWSummaryFileHead << "_" << mbw_tests->subject_name 
	   << MBWFileExtension;
	string filename = ss.str().c_str();
	path /= filename;

	ofstream outfile;
	outfile.open(path.string());
	outfile << CINIT_NAME << "," << CINIT_STD_NAME << "," << FRC_NAME << "," 
		    << FDS_NAME << "," << MACHINE_DS_NAME << "," << VBAG_NAME << ","
			<< LCI_NAME << endl;
	for(int it = 0; it < int(mbw_tests->Cinit.size()); it++)
	{
		outfile << mbw_tests->Cinit[it] << "," << mbw_tests->Cinit_std[it] << ",";
		if(it == 0)
		{
			outfile << mbw_tests->FRC0 << "," << mbw_tests->dead_space << ","
				    << mbw_tests->machine_ds << "," << mbw_tests->Vbag;
		}
		else
		{
			outfile << ",,,";
		}
		outfile << "," << mbw_tests-> LCI[it] << endl;
	}

	outfile.close();
}