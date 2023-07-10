#ifndef READ_MBW_DATA_H
#define READ_MBW_DATA_H

#include<string>
#include<vector>
#include<unordered_map>
#include<Eigen\Dense>
#include<memory>
#include"mbw_processing_params.h"

#define KELVIN0 273.15
#define MMHG_2_PA 133.32
#define BODY_TEMP_K 310.15

//file extensions
#define RAW_FILE_EXT "txt"
#define OPTIONS_FILE_EXT "options"
#define PARAMS_FILE_EXT "params"

//file headers to find
#define SF6_FHEAD "SF6"
#define MARKER_FHEAD "Markers"
#define BREATHS_FHEAD "Breaths"

//variable names in files
//marker files
#define CET_MARKER_HEADER "r_ExpMarkers[][3]"
#define C0_MARKER_HEADER "CinitMarker[][0]"
#define C0_HEADER "CinitMarker[][1]"
#define FGD_HEADER "FGD"

//breath files
#define SF6_INSP_HEADER "InspSF6"
#define VOL_INSP_HEADER "InspVol"
#define CET_HEADER "Cet"
#define SF6_EXP_HEADER "ExpSF6"
#define VOL_EXP_HEADER "ExpVol"
#define TERM_HEADER "Term"
#define FRC_HEADER "FRC"
#define TO_HEADER "TO"
#define BREATH_START_HEADER "r_ExpTime[][0]"
#define BREATH_END_HEADER "r_ExpTime[][1]"
#define P3_START_HEADER "r_BreathMarkers[][1]"
#define P3_SLOPE_HEADER "SIII" 
#define FOWLER_DS_HEADER "FDS"
#define MACHINE_DS_HEADER "DS"

//flow file
#define CONC_HEADER "SF6"
#define FLOW_HEADER "Flow"

//parameter keys
#define FLOW_DATA_KEY 1
#define CONC_DATA_KEY 2
#define CO2_DATA_KEY 3
#define O2_DATA_KEY 4

void linear_fit(const Eigen::VectorXd & x, const Eigen::VectorXd & y, double & m, double & c);
template<typename T> T median(const std::vector<T> & values, std::vector<T> & bottom_half = std::vector<T>(),
							  std::vector<T> & top_half = std::vector<T>())
{
	vector<T> valsort = values;
	std::sort(valsort.begin(), valsort.begin() + valsort.size(), less<T>());
	if(valsort.size()%2==0)  //even
	{
		size_t halfway = valsort.size()/2;
		bottom_half = vector<T>(valsort.begin(), valsort.begin() + halfway);
		top_half = vector<T>(valsort.begin() + halfway, valsort.begin() + valsort.size());
		return ((valsort[halfway - 1] + valsort[halfway])/2);   
	}
	else
	{
		size_t halfway = (valsort.size() - 1)/2;
		bottom_half = vector<T>(valsort.begin(), valsort.begin() + halfway);
		top_half = vector<T>(valsort.begin() + halfway + 1, valsort.begin() + valsort.size());
		return (valsort[halfway]);
	}
}

class MBWRawFile
{
protected:
	void trim_and_apply_delay();
	void measure_breath_start_end_points();
	void define_washin_and_washout(const double & washin_min_conc);
	std::vector<double> conc_data, CO2_data, O2_data, time, flow_data;  //store washout data
	std::vector<double> conc_data_washin, CO2_data_washin, O2_data_washin, 
		                time_washin, flow_data_washin;  //store washout data
	std::vector<size_t> inhalation_start_pts, exhalation_start_pts;
	std::vector<size_t> inhalation_start_pts_washin, exhalation_start_pts_washin;
	double CO2_delay, O2_delay;  //delay between flow and gas signal in seconds
	inline double get_full_test_data(const size_t & k, const int & param) const
	{
		if(k < this->count_washin_timesteps())
		{
			switch(param)
			{
			case FLOW_DATA_KEY:
				{
					return (this->get_washin_flow(k));
				} break;
			case CONC_DATA_KEY:
				{
					return (this->get_washin_conc(k));
				} break;
			case CO2_DATA_KEY:
				{
					return (this->get_washin_CO2(k));
				} break;
			case O2_DATA_KEY:
				{
					return (this->get_washin_O2(k));
				} break;
			}
		}
		else 
		{
			size_t k_rel = k - this->count_washin_timesteps();
			switch(param)
			{
			case FLOW_DATA_KEY:
				{
					return (this->get_washout_flow(k_rel));
				} break;
			case CONC_DATA_KEY:
				{
					return (this->get_washout_conc(k_rel));
				} break;
			case CO2_DATA_KEY:
				{
					return (this->get_washout_CO2(k_rel));
				} break;
			case O2_DATA_KEY:
				{
					return (this->get_washout_O2(k_rel));
				} break;
			}
		}
		return 0;
	}
public:
	MBWRawFile(const std::string & fname, const double & CO2d, 
		       const double & O2d, const double & washin_min_conc)    //construct from raw files
	{
		this->CO2_delay = CO2d;
		this->O2_delay = O2d;
		read_raw_data(fname, washin_min_conc);
	}

	void read_raw_data(const std::string & fname, const double & washin_min_conc);

	inline size_t count_washin_timesteps() const { return conc_data_washin.size(); }
	inline size_t count_washout_timesteps() const { return conc_data.size(); }
	inline size_t count_washin_inhalations() const { return inhalation_start_pts_washin.size(); }
	inline size_t count_washin_exhalations() const { return exhalation_start_pts_washin.size(); }
	inline size_t count_washout_inhalations() const { return inhalation_start_pts.size(); }
	inline size_t count_washout_exhalations() const { return exhalation_start_pts.size(); }
	inline size_t get_washin_inhalation_start(const size_t & n) const 
	{ 
		return this->inhalation_start_pts_washin[n]; 
	}
	inline size_t get_washin_exhalation_start(const size_t & n) const 
	{ 
		return this->exhalation_start_pts_washin[n]; 
	}
	inline size_t get_washout_inhalation_start(const size_t & n) const 
	{ 
		return this->inhalation_start_pts[n]; 
	}
	inline size_t get_washout_exhalation_start(const size_t & n) const 
	{ 
		return this->exhalation_start_pts[n]; 
	}

	inline double get_washin_flow(const size_t & k) const { return flow_data_washin[k]; }
	inline double get_washout_flow(const size_t & k) const { return flow_data[k]; }
	inline double get_flow(const size_t & k) const {return (this->get_full_test_data(k, FLOW_DATA_KEY)); }
	inline double get_washin_conc(const size_t & k) const { return conc_data_washin[k]; }
	inline double get_washout_conc(const size_t & k) const { return conc_data[k]; }
	inline double get_conc(const size_t & k) const {return (this->get_full_test_data(k, CONC_DATA_KEY)); }
	inline double get_washin_CO2(const size_t & k) const { return CO2_data_washin[k]; }
	inline double get_washout_CO2(const size_t & k) const { return CO2_data[k]; }
	inline double get_CO2(const size_t & k) const {return (this->get_full_test_data(k, CO2_DATA_KEY)); }
	inline double get_washin_O2(const size_t & k) const { return O2_data_washin[k]; }
	inline double get_washout_O2(const size_t & k) const { return O2_data[k]; }
	inline double get_O2(const size_t & k) const {return (this->get_full_test_data(k, O2_DATA_KEY)); }
};

std::shared_ptr<LCIParameterList> get_raw_mbw_data(const int &argc, char** argv, 
					 std::vector<std::shared_ptr<MBWRawFile>> & mbw_files, 
					 LCIOptionList *options);

//class for processing raw data
struct MBWTest
{
private:
	void process_MBW_inputs(std::shared_ptr<MBWRawFile> MBW_raw_file, LCIOptionList *opt,
		LCIParameterList *par, const int & Test_no);
	void adjust_flow_for_BTPS(std::shared_ptr<MBWRawFile> MBWfile, 
		       Eigen::VectorXd & adjusted_washin_flow, Eigen::VectorXd & adjusted_washout_flow, 
			   const double & BTPSin);
	void measure_Cinit(std::shared_ptr<MBWRawFile> MBWFile, const Eigen::VectorXd & washin_flux, const double & Cinit_input);
	void measure_washout_breath_volumes(std::shared_ptr<MBWRawFile> MBWFile, 
		                                const Eigen::VectorXd & washout_flux);
	void measure_LCI_FRC_FDS(std::shared_ptr<MBWRawFile>, const Eigen::VectorXd & washout_flux,
		                     const double & LCI_frac);
	void define_breath_measurements(std::shared_ptr<MBWRawFile> MBWFile, 
									const Eigen::VectorXd & washout_flux,
									const bool & measure_subset, const bool & fixed_step_size,
									const bool & measure_inspired_vols,
									const int & Nmeasure, const int & NphaseII, 
									const int & NphaseIII, const double & vol_sim_step_frac);
public:
	int Nbreaths, LCI_point;
	double Cinit, Cinit_std, median_FDS, FRC0, LCI;  //MBW test measurements
	std::vector<size_t> inhaled_breath_pts, exhaled_breath_pts,  //number of time points in each breath
		                conc_measurement_steps, igvol_measurement_steps;      //stores location sim_steps where measurements take place
	std::vector<double> inhaled_breath_vols, exhaled_breath_vols,    //breath volumes
		                cumulative_exhaled_vols, cumulative_inhaled_vols,   //cumulative volumes
						Cet, CO2et, FDS,    //store standard measures
						conc_measurements, igvol_measurements,    //store measurements for distance function
						sim_vol_steps, sim_step_durations,   //stores info for simulation timesteps
						exhaled_SF6vols, inhaled_SF6vols, cumulative_exhaled_SF6vols,    //SF6 breath volumes
						cumulative_inhaled_SF6vols;   //cumulative SF6 volumes

	MBWTest(std::shared_ptr<MBWRawFile> MBW_raw_file, LCIOptionList *opt,
		LCIParameterList *par, const int & Test_no)
	{
		this->process_MBW_inputs(MBW_raw_file, opt, par, Test_no);
	}
};

//class for storing multiple tests
struct MBWData
{
private:
	void process_data_inputs(const std::vector<std::shared_ptr<MBWRawFile>> & MBWfiles,
		 LCIOptionList * opt, LCIParameterList * par);

public:
	std::vector<std::vector<double>> sim_vol_steps, sim_step_durations,
		conc_measurements, igvol_measurements;
	std::vector<std::vector<size_t>> conc_measurement_steps, igvol_measurement_steps;
	std::vector<double> Cinit, Cinit_std, LCI;
	double dead_space, FRC0, machine_ds, extra_rebreathe_vol, Vbag, VT0, Tinsp, Texp;     //timestep and median fowler dead-space
	//Eigen::MatrixXd weights;
	std::string subject_name;

	MBWData()
	{
	}
	MBWData(const std::vector<std::shared_ptr<MBWRawFile>> & MBWfiles, 
		    LCIOptionList * opt, LCIParameterList * par)
	{
		this->process_data_inputs(MBWfiles, opt, par);
	}
	~MBWData()
	{	
	};
};

void write_processed_washout_data(const std::string & filepath, const MBWData* mbw_data);

void write_mbw_summary(const std::string & filepath, const MBWData* mbw_data);

inline double PH20_Buck(const double & HumidityFrac, const double & TempC)
{
	double P = 0.61121*exp((18.678-(TempC/234.5))*(TempC/(257.14+TempC)));
	return HumidityFrac*P;
}

double calc_BTPS(LCIParameterList *par, const int & Test_no);

#endif