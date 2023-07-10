#include"MBW_models.h"
#include"MBW_model_params.h"
#include"read_write_codes.h"
#include<boost\filesystem.hpp>
#include<boost\regex.hpp>
#include<file_manip.h>
#include<map>

//TODO
//using NCOMPS in parts and opt.Nunits in others, needs to be consistent
// 


extern std::shared_ptr<boost::random::mt19937> rng;  //pseudorandom number generator

void MBWModelInputs::read_inputs(const std::string & filepath)
{
	//filepath in this case should be a directory containing the following info
	//..._MBW_summary_data.csv
	//..._MBW_test1.csv
	//..._MBW_test2.csv
	//etc
	//get all file paths
	std::stringstream ss;
	ss << MBWTestFileHead << ".*\\" << MBWFileExtension;
	boost::regex MBWTestRegex(ss.str().c_str());

	ss.clear();
	ss.str("");
	ss << MBWSummaryFileHead << ".*\\" << MBWFileExtension;
	boost::regex MBWSummaryRegex(ss.str().c_str());

	std::vector<boost::filesystem::path> MBWTestPaths, MBWSummaryPaths;
	boost::filesystem::path p_in(filepath.c_str()); 
	boost::filesystem::directory_iterator it_end; //yields beyond end iterator 
	for(boost::filesystem::directory_iterator it(p_in); it!=it_end; ++it)
	{
		boost::smatch result1, result2;
		if(boost::regex_match(it->path().filename().string(), result1, MBWTestRegex)) //returns true for match
		{
			MBWTestPaths.push_back(it->path());
		}
		if(boost::regex_match(it->path().filename().string(), result2, MBWSummaryRegex)) //returns true for match
		{
			MBWSummaryPaths.push_back(it->path());
		}
	}
	
	//sort MBW Test Paths
	std::vector<int> MBWTestNumbers, Indices;
	MBWTestNumbers.resize(MBWTestPaths.size());
	Indices.resize(MBWTestPaths.size());
	for(size_t n = 0; n < MBWTestPaths.size(); n++)
	{
		std::string filehead;
		size_t lastdot = MBWTestPaths[n].filename().string().find_last_of(".");
		if (lastdot == std::string::npos) filehead = MBWTestPaths[n].filename().string();
		filehead = MBWTestPaths[n].filename().string().substr(0, lastdot); 
		std::vector<std::string> spst = string_split(filehead,"_");
		MBWTestNumbers[n] = StringToNumber<int>(spst.back());
		Indices[n] = n;
	}
	std::sort(Indices.begin(), Indices.end(),
       [&MBWTestNumbers](size_t i1, size_t i2) {return MBWTestNumbers[i1] < MBWTestNumbers[i2];});
	std::vector<boost::filesystem::path> MBWTestPathsNew;
	MBWTestPathsNew.resize(MBWTestPaths.size());
	for(size_t n = 0; n < MBWTestPaths.size(); n++)
	{
		MBWTestPathsNew[n] = MBWTestPaths[Indices[n]];
	}
	MBWTestPaths = MBWTestPathsNew;

	//read test data
	this->Ntests = int(MBWTestPaths.size());
	this->conc_measurements.resize(this->Ntests);
	this->conc_weights.resize(this->Ntests);
	this->igvol_diff_measurements.resize(this->Ntests);
	this->igvol_weights.resize(this->Ntests);
	this->conc_measurement_steps.resize(this->Ntests);
	this->igvol_measurement_steps.resize(this->Ntests);
	this->sim_vol_steps.resize(this->Ntests);
	this->sim_step_durations.resize(this->Ntests);
	this->av_vol_step = 0;
	unsigned int i_count = 0;
	for(int n = 0; n < this->Ntests; n++)
	{
		std::vector<std::string> infile = get_all_lines(MBWTestPaths[n].string());
		std::vector<std::string> headers = string_split(infile[0], ",");
		std::map<std::string,int> col_numbers;
		for(int ih = 0; ih < int(headers.size()); ih++)
		{
			col_numbers[headers[ih]] = ih;
		}
		double old_cumul_igvol = 0, cumul_igvolh = 0;
		double old_cumul_volTO = 0, cumul_volTO = 0;
		for(int iline = 1; iline < int(infile.size()); iline++)
		{
			std::vector<std::string> line = string_split(infile[iline], ",");
			double vol_step = StringToNumber<double>(line[col_numbers.at(VOLUME_NAME)]);
			this->sim_vol_steps[n].push_back(vol_step);
			cumul_volTO += abs(vol_step);
			this->av_vol_step += abs(this->sim_vol_steps[n][iline-1]);
			this->sim_step_durations[n].push_back(StringToNumber<double>(line[col_numbers.at(DURATION_NAME)]));
			//get conc measurements
			this->conc_measurement_steps[n].push_back(iline-1);
			this->conc_measurements[n].push_back(StringToNumber<double>(line[col_numbers.at(CONC_NAME)]));
			if(USE_SF6_CONC) this->conc_weights[n].push_back(StringToNumber<double>(line[col_numbers.at(CONC_MEASURED_NAME)]));
			else this->conc_weights[n].push_back(0.0);
			//get igvol diff measurements
			this->igvol_measurement_steps[n].push_back(iline-1);
			cumul_igvolh = StringToNumber<double>(line[col_numbers.at(IGVOL_NAME)]);
			this->igvol_diff_measurements[n].push_back((cumul_igvolh - old_cumul_igvol));
			old_cumul_igvol = cumul_igvolh;
			old_cumul_volTO = cumul_volTO;
			if(USE_SF6_VOL) this->igvol_weights[n].push_back(StringToNumber<double>(line[col_numbers.at(IGVOL_MEASURED_NAME)]));
			else this->igvol_weights[n].push_back(0.0);
		}
		i_count += this->sim_vol_steps[n].size();
	}
	this->av_vol_step /= double(i_count);
	//read summary data
	if(MBWSummaryPaths.size() > 1) 
	{
		std::cerr << "Warning, more than one summary file, only reading in "
			<< MBWSummaryPaths[0].filename().string() << std::endl;
	}
	if(MBWSummaryPaths.size() == 0) 
	{
		std::cerr << "No summary file found, aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<std::string> infile = get_all_lines(MBWSummaryPaths[0].string());
	std::vector<std::string> headers = string_split(infile[0], ",");
	std::map<std::string,int> col_numbers;
	for(int ih = 0; ih < int(headers.size()); ih++)
	{
		col_numbers[headers[ih]] = ih;
	}
	std::vector<std::string> topline = string_split(infile[1], ",");
	this->FRC0 = StringToNumber<double>(topline[col_numbers.at(FRC_NAME)]);
	this->machine_ds = StringToNumber<double>(topline[col_numbers.at(MACHINE_DS_NAME)]);
	this->dead_space = StringToNumber<double>(topline[col_numbers.at(FDS_NAME)]);
	this->rebreathe_vol = StringToNumber<double>(topline[col_numbers.at(REBREATHE_VOL_NAME)]);
	//our definition of FRC in this model excludes dead space
	this->FRC0 -= this->dead_space;   //IMPORTANT

	//get cinit vals
	this->Cinit.resize(this->Ntests);
	this->Cinit_std.resize(this->Ntests);
	for(int it = 0; it < this->Ntests; it++)
	{
		std::vector<std::string> line = string_split(infile[1+it], ",");
		this->Cinit[it] = StringToNumber<double>(line[col_numbers.at(CINIT_NAME)]);
		this->Cinit_std[it] = StringToNumber<double>(line[col_numbers.at(CINIT_STD_NAME)]);
	}

	//get subject name
	std::vector<std::string> filehead = string_split(MBWSummaryPaths[0].filename().string(),MBWFileExtension);
	std::vector<std::string> fileparts = string_split(filehead[0], "_");
	ss.clear();
	ss.str("");
	for(int ip = 1; ip < int(fileparts.size()); ip++)
	{
		if(ip > 1) ss << "_";
		ss << fileparts[ip];
	}
	this->subject_name = ss.str().c_str();

	//needs to be reformatted into single vector
	for(int it = 0; it < this->Ntests; it++)
	{
		this->measured.insert(this->measured.end(),this->conc_measurements[it].begin(), 
			                                       this->conc_measurements[it].end());
		this->dist_weights.insert(this->dist_weights.end(), this->conc_weights[it].begin(), this->conc_weights[it].end());
		this->measured.insert(this->measured.end(),this->igvol_diff_measurements[it].begin(), 
			                                       this->igvol_diff_measurements[it].end());
		this->dist_weights.insert(this->dist_weights.end(), this->igvol_weights[it].begin(), this->igvol_weights[it].end());
	}
}

void MBWModelInputs::generate_inputs(const std::string & params_file, MBWModelOptions & opts, 
									 std::vector<double> & MBWModelParams,  std::vector<std::string> & MBWModelParamNames)
{
	//write this
	//params: 
	//vent_dist_type: 'l' or 'n'
	//lung_unit_type: 's' or 'a'
	//Nunits: int (default 50)
	//NMR_samples: int (default 1000)
	//mix_vol_frac_step: float (default 0.2)
	//model params
	//FRC, VD, VD_shared_fraction, 

	//this could be a read options & params function
	//this->Ntests=3;
	//std::shared_ptr<MBWOptionList> options = std::make_shared<MBWOptionList>();
	//std::shared_ptr<MBWParameterList> params = std::make_shared<MBWParameterList>();
	//params->add_start_inflation_params(this->Ntests);

	//std::vector<std::string> infile = get_all_lines(params_file);
	//for(unsigned int i = 0; i < infile.size(); i++)
	//{
	//	std::vector<std::string> line = string_split(infile[i]," ");
	//	if(line.size() > 1)
	//	{
	//		if (!(options->parse(line[0],line[1]))) 
	//		{
	//			if(!(params->parse(line[0],line[1])))
	//			{
	//				std::cerr << "Error, did not recognise option " << line[0] << std::endl;
	//			}
	//		}
	//	}
	//}
	//params->check_validity(options.get());
	////now use these to create input file

	////generate flux -- can be moved to another function
	////first generate breath sizes
	//double bp = params->get_param<double>(BPERIOD_PARAM_NAME)->get_value();
	//double vt = params->get_param<double>(VT_PARAM_NAME)->get_value();
	//double dvt = params->get_param<double>(VTFRANGE_PARAM_NAME)->get_value();

	//int Nwashin = MBW_BMAX;
	//int Nwashin_steps = int(Nwashin*bp*100);
	//int Nwashout = MBW_BMAX;
	//int Nwashout_steps = int(Nwashout*bp*100);
	//
	//boost::random::uniform_01<double> u01;
	//for(int n = 0; n < this->Ntests; n++)
	//{
	//	std::vector<double> breath_vols(Nwashin+Nwashout, 0.0);
	//	for(int j = 0; j < Nwashin; j++)
	//	{
	//		breath_vols[j] = 2.0*(vt + 2*dvt*(u01(*(rng.get())) - 0.5)); 
	//		//washin breath size is double
	//	}
	//	for(int j = Nwashin; j < Nwashin + Nwashout; j++)
	//	{
	//		breath_vols[j] = vt + 2*dvt*(u01(*(rng.get())) - 0.5);
	//	}
	//	this->sim_step_durations.push_back(std::vector<double>(Nwashin_steps+Nwashout_steps, 0.01));
	//	this->sim_vol_steps.push_back(std::vector<double>(Nwashin_steps+Nwashout_steps, 0.0));
	//	this->measurement_steps.push_back(std::vector<int>(Nwashin_steps+Nwashout_steps, 0));
	//	for(int i = 0; i < Nwashin_steps+Nwashout_steps; i++)
	//	{
	//		int nb = int(i*0.01/bp);
	//		if(i*0.01/bp - nb < 0.5) //inhale
	//		{
	//			this->sim_vol_steps[n][i] = 0.01*breath_vols[nb]/(0.5*bp);
	//		}
	//		else //exhale
	//		{
	//			this->sim_vol_steps[n][i] = -0.01*breath_vols[nb]/(0.5*bp);
	//		}
	//		this->measurement_steps[n][i] = i;
	//	}
	//	this->Cinit.push_back(0.2);  //in %
	//	this->Cinit_std.push_back(SF6_NOISE);
	//}
	//this->machine_ds = params->get_param<double>(MACHINE_DS_PARAM_NAME)->get_value();
	//this->av_vol_step = 2*0.01*vt/bp;
	////create and return options file
	//opts.lung_unit_type = options->get_option<char>(LUNG_UNIT_KEY)->get_value();
	//opts.vent_dist_type = options->get_option<char>(VDIST_KEY)->get_value();
	//opts.breath_model_type = options->get_option<char>(BREATHING_MODEL_KEY)->get_value();
	//opts.NMR_samples = NMRI_SAMPLES;            //could be user inputted
	//opts.Nunits = NCOMPS;                       //could be user inputted
	//opts.mix_vol_frac_step = NMV_FRAC_STEP;     //could be user inputted
	//opts.simulate_washin = true;
	//opts.washout_start_timepoint = Nwashin_steps;


	////create params for MBW model
	//for(int n = 1; n < this->Ntests; n++)
	//{
	//	std::stringstream ss;
	//	ss << FRC_TEST_MODIFIER_NAME << '_' << n;
	//	std::string pname = ss.str().c_str();
	//	MBWModelParams.push_back(params->get_param<double>(pname)->get_value());
	//	MBWModelParamNames.push_back(pname);
	//}

	//MBWModelParams.push_back(params->get_param<double>(FRC_PARAM_NAME)->get_value());
	//MBWModelParamNames.push_back(FRC_PARAM_NAME);
	//MBWModelParams.push_back(params->get_param<double>(VD_PARAM_NAME)->get_value());
	//MBWModelParamNames.push_back(VD_PARAM_NAME);
	//MBWModelParams.push_back(params->get_param<double>(SIGMA_PARAM_NAME)->get_value());
	//MBWModelParamNames.push_back(SIGMA_PARAM_NAME);
	//bool SDS = (params->get_param<double>(VDSFRAC_PARAM_NAME)->get_value() > 0);
	//if(SDS)
	//{
	//	MBWModelParams.push_back(params->get_param<double>(VDSFRAC_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(VDSFRAC_PARAM_NAME);
	//}
	//if(opts.vent_dist_type == BIMODAL_CODE)
	//{
	//	MBWModelParams.push_back(params->get_param<double>(MURATIO_LS_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(MURATIO_LS_PARAM_NAME);
	//	MBWModelParams.push_back(params->get_param<double>(SIGRATIO_LS_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(SIGRATIO_LS_PARAM_NAME);
	//	MBWModelParams.push_back(params->get_param<double>(VFASTFRAC_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(VFASTFRAC_PARAM_NAME);
	//}
	//if(opts.lung_unit_type == ASYMM_UNIT_CODE)
	//{
	//	MBWModelParams.push_back(params->get_param<double>(ASYMM_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(ASYMM_PARAM_NAME);
	//	MBWModelParams.push_back(params->get_param<double>(DIFFSCALE_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(DIFFSCALE_PARAM_NAME);
	//}
	//if(opts.breath_model_type == ASYNC_MODEL_CODE)
	//{
	//	MBWModelParams.push_back(params->get_param<double>(DELAY_PARAM_NAME)->get_value());
	//	MBWModelParamNames.push_back(DELAY_PARAM_NAME);
	//}
}

bool MBWModelOutputs::extra_outputs() const
{
	return true;
}

void MBWModelOutputs::print_extra_outputs(std::string & line) const
{
	std::stringstream ss;
	for(int iV = 0; iV < NCOMPS; iV++)
	{
		if(iV > 0) ss << ",";
		ss << this->Vrates[iV];
	}
	for(int iV = 0; iV < NCOMPS; iV++)
	{
		ss << ",";
		ss << this->Vdelays[iV];
	}
	for(int isim = 0; isim < int(this->simulated.size()); isim++)
	{
		ss << ",";
		ss << this->simulated[isim];
	}
	line = ss.str().c_str();
}

void MBWModelOutputs::get_headers(std::string & line) const
{
	std::stringstream ss;
	for(int iV = 0; iV < NCOMPS; iV++)
	{
		if(iV > 0) ss << ",";
		ss << "VentRate" << iV;
	}
	for(int iV = 0; iV < NCOMPS; iV++)
	{
		ss << ",";
		ss << "VentDelay" << iV;
	}
	for(int isim = 0; isim < int(this->simulated.size()); isim++)
	{
		ss << ",";
		ss << "MBW_sample_" << isim;
	}
	line = ss.str().c_str();
}

MBWModelOptions::MBWModelOptions()   //default model options
{
	this->Nunits = NCOMPS;
	this->mix_vol_frac_step = NMV_FRAC_STEP;
	this->washout_start_inflation = 0.0;
	this->washout_start_inflation_duration = 0.0;
	this->simulate_washin = false;
	this->washout_start_timepoint = 0;
}

void AsymmLungUnit::inhale(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & fv, 
		                   const double & dt) //absorb and destroy volume element
{
	//should really solve simultaneous step -> should be simple enough
	this->conc_old = this->conc;
	//advection step
	for(size_t n = 0; n < fv.size(); n++)   //add elements in order inhaled
	{
		this->vol1 += 0.5*fv[n]->get_volume();
		this->vol2 += 0.5*fv[n]->get_volume();
		this->IGvol1 += 0.5*fv[n]->get_igvol();
		this->IGvol2 += 0.5*fv[n]->get_igvol();
	}
	this->conc1 = IGvol1 / vol1;
	this->conc2 = IGvol2 / vol2;
	//diffusion step   --no diff on inhale
	double dc = conc1 - conc2;
	double dcend;
	if(DT > 0 && vol1 > 0 && vol2 > 0) dcend = dc*exp(-0.5*FRC*(1.0/vol1 + 1/vol2)*dt/DT);
	else dcend = 0;
	this->conc1 = (IGvol1 + IGvol2)/(vol1 + vol2) + dcend/(1 + vol1/vol2);
	this->conc2 = conc1 - dcend;
	this->IGvol1 = vol1*conc1;
	this->IGvol2 = vol2*conc2;

	//update whole acinar params
	this->volume = vol1 + vol2;
	this->IGvolume = IGvol1 + IGvol2;
	this->conc = 0.5*conc1 + 0.5*conc2;
}
					
void AsymmLungUnit::exhale(const double & dv, std::vector<std::shared_ptr<FlexibleVolumeElement>> & exhaled, 
			               const double & dt)
{
	this->conc_old = this->conc;
	//advection step -- exhale unit with constant conc
	this->conc = 0.5*conc1 + 0.5*conc2;
	double dv1 = std::min(0.5*dv,(1-1E-06)*vol1);
	double dv2 = dv - dv1;
	this->vol1 -= dv1;
	this->vol2 -= dv2;
	this->IGvol1 -= conc1*dv1;
	this->IGvol2 -= conc2*dv2;
	
	exhaled.resize(1);
	exhaled[0] = std::make_shared<FlexibleVolumeElement>(dv, this->conc, this->conc);

	//diffusion step
	double dc = conc1 - conc2;
	double dcend;
	if(DT > 0 && vol1 > 0 && vol2 > 0) dcend = dc*exp(-0.5*FRC*(1.0/vol1 + 1/vol2)*dt/DT);
	else dcend = 0;
	this->conc1 = (IGvol1 + IGvol2)/(vol1 + vol2) + dcend/(1 + vol1/vol2);
	this->conc2 = conc1 - dcend;
	this->IGvol1 = vol1*conc1;
	this->IGvol2 = vol2*conc2;

	//update whole acinar params
	this->volume = vol1 + vol2;
	this->IGvolume = IGvol1 + IGvol2;
	this->conc = 0.5*conc1 + 0.5*conc2;
}

void CompartmentalModelGeneratorBase::generate_from_prior(std::vector<double> & 
													     params_generated) const
{
	boost::random::uniform_01<> udist;
	params_generated.resize(this->param_names.size());
	for(int ip = 0; ip < int(this->param_names.size()); ip++)
	{
		params_generated[ip] = this->param_min[ip] + (this->param_max[ip]
		                        - this->param_min[ip])*udist(*(rng.get()));
	}
}
	
double CompartmentalModelGeneratorBase::prior_density(const std::vector<double> &
												      params) const
{
	for(int ip = 0; ip < int(this->param_names.size()); ip++)
	{
		if(params[ip] < this->param_min[ip] || params[ip] > this->param_max[ip])
		{
			return 0.0;		
		}
	}

	return 1.0;
}

void CompartmentalModelGeneratorBase::initialise_FRC_modifier_params(MBWModelInputs* inputs)
{
	size_t Ntests = inputs->Cinit.size();
	this->param_names.resize(Ntests-1);
	this->param_min.resize(Ntests-1);
	this->param_max.resize(Ntests-1);

	for(size_t n = 1; n < Ntests; n++)
	{
		std::stringstream ss;
		ss << FRC_TEST_MODIFIER_NAME << "_" << n;
		this->param_names[n-1] = ss.str().c_str();
		this->param_min[n-1] = FRC_TEST_MODIFIER_MIN;
		this->param_max[n-1] = FRC_TEST_MODIFIER_MAX;
	}
	this->N_FRC_mod_params = Ntests-1;
}

void BasicLognormalModelGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+3);
	this->param_min.resize(this->N_FRC_mod_params+3);
	this->param_max.resize(this->N_FRC_mod_params+3);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;

	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
}

void BasicLognormalModelGenerator::set_model_name()
{
	model_name = "Lognormal_NoSDS";
}

void CompartmentalModelGeneratorBase::generate_model(const std::vector<double> & params,
		                                             std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = LOGNORMAL_CODE;
	opts.lung_unit_type = BASIC_UNIT_CODE;
	opts.breath_model_type = SYNC_MODEL_CODE;
	m->build_model(opts, params, this->param_names);
}

void LognormalModelSDSGenerator::set_model_name()
{
	model_name = "Lognormal_SDS";
}

void LognormalModelSDSGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+4);
	this->param_min.resize(this->N_FRC_mod_params+4);
	this->param_max.resize(this->N_FRC_mod_params+4);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+3] = VDSFRAC_PARAM_NAME;

	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;
	this->param_min[this->N_FRC_mod_params+3] = 0;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
	this->param_max[this->N_FRC_mod_params+3] = 1;
}

void BimodalModelSDSGenerator::set_model_name()
{
	model_name = "Bimodal_SDS";
}

void BimodalModelSDSGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+7);
	this->param_min.resize(this->N_FRC_mod_params+7);
	this->param_max.resize(this->N_FRC_mod_params+7);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+3] = VDSFRAC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+4] = MURATIO_LS_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+5] = SIGRATIO_LS_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+6] = VFASTFRAC_PARAM_NAME;


	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;
	this->param_min[this->N_FRC_mod_params+3] = 0;
	this->param_min[this->N_FRC_mod_params+4] = 1.0;   //ratio of large to small, cannot be < 1
	this->param_min[this->N_FRC_mod_params+5] = 1.0/SIGRATIO_MAX;
	this->param_min[this->N_FRC_mod_params+6] = 0.0;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
	this->param_max[this->N_FRC_mod_params+3] = 1;
	this->param_max[this->N_FRC_mod_params+4] = MURATIO_MAX;   //ratio of large to small, cannot be < 1
	this->param_max[this->N_FRC_mod_params+5] = SIGRATIO_MAX;
	this->param_max[this->N_FRC_mod_params+6] = 1.0;
}

void BimodalModelSDSGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = BIMODAL_CODE;
	opts.lung_unit_type = BASIC_UNIT_CODE;
	opts.breath_model_type = SYNC_MODEL_CODE;
	m->build_model(opts, params, this->param_names);
}

void LognormalAsymmSDSModelGenerator::set_model_name()
{
	model_name = "Lognormal_SDS_Asymm";
}

void LognormalAsymmSDSModelGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+6);
	this->param_min.resize(this->N_FRC_mod_params+6);
	this->param_max.resize(this->N_FRC_mod_params+6);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+3] = VDSFRAC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+4] = ASYMM_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+5] = DIFFSCALE_PARAM_NAME;

	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;
	this->param_min[this->N_FRC_mod_params+3] = 0;
	this->param_min[this->N_FRC_mod_params+4] = 0.0;
	this->param_min[this->N_FRC_mod_params+5] = DIFFSCALE_MIN;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
	this->param_max[this->N_FRC_mod_params+3] = 1;
	this->param_max[this->N_FRC_mod_params+4] = 1.0;
	this->param_max[this->N_FRC_mod_params+5] = DIFFSCALE_MAX;
}

void LognormalAsymmSDSModelGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = LOGNORMAL_CODE;
	opts.lung_unit_type = ASYMM_UNIT_CODE;
	opts.breath_model_type = SYNC_MODEL_CODE;
	m->build_model(opts, params, this->param_names);
}

void BasicLognormalAsymmModelGenerator::set_model_name()
{
	model_name = "Lognormal_NoSDS_Asymm";
}

void BasicLognormalAsymmModelGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+5);
	this->param_min.resize(this->N_FRC_mod_params+5);
	this->param_max.resize(this->N_FRC_mod_params+5);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+3] = ASYMM_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+4] = DIFFSCALE_PARAM_NAME;

	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;
	this->param_min[this->N_FRC_mod_params+3] = 0.0;
	this->param_min[this->N_FRC_mod_params+4] = DIFFSCALE_MIN;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
	this->param_max[this->N_FRC_mod_params+3] = 1.0;
	this->param_max[this->N_FRC_mod_params+4] = DIFFSCALE_MAX;
}

void BasicLognormalAsymmModelGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = LOGNORMAL_CODE;
	opts.lung_unit_type = ASYMM_UNIT_CODE;
	opts.breath_model_type = SYNC_MODEL_CODE;
	m->build_model(opts, params, this->param_names);
}

void BimodalAsymmSDSModelGenerator::set_model_name()
{
	model_name = "Bimodal_SDS_Asymm";
}

void BimodalAsymmSDSModelGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+9);
	this->param_min.resize(this->N_FRC_mod_params+9);
	this->param_max.resize(this->N_FRC_mod_params+9);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+3] = VDSFRAC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+4] = MURATIO_LS_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+5] = SIGRATIO_LS_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+6] = VFASTFRAC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+7] = ASYMM_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+8] = DIFFSCALE_PARAM_NAME;

	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;
	this->param_min[this->N_FRC_mod_params+3] = 0;
	this->param_min[this->N_FRC_mod_params+4] = 1.0;   //ratio of large to small, cannot be < 1
	this->param_min[this->N_FRC_mod_params+5] = 1.0/SIGRATIO_MAX;
	this->param_min[this->N_FRC_mod_params+6] = 0.0;
	this->param_min[this->N_FRC_mod_params+7] = 0.0;
	this->param_min[this->N_FRC_mod_params+8] = DIFFSCALE_MIN;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
	this->param_max[this->N_FRC_mod_params+3] = 1;
	this->param_max[this->N_FRC_mod_params+4] = MURATIO_MAX;   //ratio of large to small, cannot be < 1
	this->param_max[this->N_FRC_mod_params+5] = SIGRATIO_MAX;
	this->param_max[this->N_FRC_mod_params+6] = 1.0;
	this->param_max[this->N_FRC_mod_params+7] = 1.0;
	this->param_max[this->N_FRC_mod_params+8] = DIFFSCALE_MAX;
}

void BimodalAsymmSDSModelGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = BIMODAL_CODE;
	opts.lung_unit_type = ASYMM_UNIT_CODE;
	opts.breath_model_type = SYNC_MODEL_CODE;
	m->build_model(opts, params, this->param_names);
}

void BasicLognormalAsyncModelGenerator::set_model_name()
{
	model_name = "Basic_lognormal_async";
}

void BasicLognormalAsyncModelGenerator::initialise_params()
{
	this->param_names.resize(this->N_FRC_mod_params+4);
	this->param_min.resize(this->N_FRC_mod_params+4);
	this->param_max.resize(this->N_FRC_mod_params+4);

	this->param_names[this->N_FRC_mod_params+0] = FRC_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+1] = VD_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+2] = SIGMA_PARAM_NAME;
	this->param_names[this->N_FRC_mod_params+3] = DELAY_PARAM_NAME;

	this->param_min[this->N_FRC_mod_params+0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[this->N_FRC_mod_params+1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[this->N_FRC_mod_params+2] = SIGMA_MIN;
	this->param_min[this->N_FRC_mod_params+3] = 0;

	this->param_max[this->N_FRC_mod_params+0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[this->N_FRC_mod_params+1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[this->N_FRC_mod_params+2] = SIGMA_MAX;
	this->param_max[this->N_FRC_mod_params+3] = DELAY_MAX_SECS;
}

void BasicLognormalAsyncModelGenerator::generate_model(const std::vector<double> & params, 
		                std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = LOGNORMAL_CODE;
	opts.lung_unit_type = BASIC_UNIT_CODE;
	opts.breath_model_type = ASYNC_MODEL_CODE;
	m->build_model(opts, params, this->param_names);
}


void MBWModelBase::set_input_data(const MBWModelInputs * inputs)
{
	//pointers to constant quantities first (no need to copy locally)
	this->sim_step_durations = &(inputs->sim_step_durations);
	this->conc_measurement_steps = &(inputs->conc_measurement_steps);
	this->igvol_measurement_steps = &(inputs->igvol_measurement_steps);


	boost::random::normal_distribution<> volndist(0.0, VOL_NOISE*inputs->av_vol_step);  //fractional noise
	this->Cinit.resize(inputs->Ntests);
	this->sim_vol_steps.resize(inputs->Ntests);
	this->sim_volTO_cumul.resize(inputs->Ntests);
	for(int it = 0; it < inputs->Ntests; it++)
	{
		//perturb Cinit
		boost::random::normal_distribution<> Cinitndist(0, inputs->Cinit_std[it]);  //absolute noise
		this->Cinit[it] = inputs->Cinit[it] + Cinitndist(*(rng.get()));
		//perturb vol steps
		this->sim_vol_steps[it].resize(inputs->sim_vol_steps[it].size());
		this->sim_volTO_cumul[it].resize(inputs->sim_vol_steps[it].size());
		double last_volTO = 0;
		for(int is = 0; is < int(inputs->sim_vol_steps[it].size()); is++)
		{
			this->sim_vol_steps[it][is] = inputs->sim_vol_steps[it][is] + volndist(*(rng.get()));
			this->sim_volTO_cumul[it][is] = last_volTO + abs(this->sim_vol_steps[it][is]);
			last_volTO = this->sim_volTO_cumul[it][is];
		}
	}
	this->mouth_point = inputs->rebreathe_vol;
	this->conc_measurement_point = inputs->rebreathe_vol - inputs->machine_ds; 
}

void MBWModelBase::build_model(const MBWModelOptions & opt,
							   const std::vector<double> & params_h,
		                       const std::vector<std::string> param_names_h)
{
	using namespace std;
	map<std::string, double> param_dict;
	for(int ip = 0; ip < int(params_h.size()); ip++)
	{
		param_dict[param_names_h[ip]] = params_h[ip];
	}

	//process parameters
	double VDSfrac = 0;
	double VD = param_dict.at(VD_PARAM_NAME);
	if(param_dict.find(VDSFRAC_PARAM_NAME) != param_dict.end()) 
	{
		VDSfrac = param_dict.at(VDSFRAC_PARAM_NAME);
	}

	double Vbag = (param_dict.at(FRC_PARAM_NAME))/ ((double) opt.Nunits);
	//assign function pointers based on options
	//ventilation dist options
	if(opt.vent_dist_type == LOGNORMAL_CODE) this->generate_vent_dist = &generate_vent_dist_lognormal_rand;
	else if(opt.vent_dist_type == BIMODAL_CODE) this->generate_vent_dist = &generate_vent_dist_bimodal_rand;
	//lung unit options
	if(opt.lung_unit_type == BASIC_UNIT_CODE) this->initialise_lung_unit = &build_basic_lung_unit;
	else if(opt.lung_unit_type == ASYMM_UNIT_CODE) this->initialise_lung_unit = &build_asymm_lung_unit;

	//generate ventilation dist
	vector<double> Vratios;
	this->generate_vent_dist(param_dict, opt.Nunits, Vratios);
	//build lung units
	this->units.resize(opt.Nunits);
	for(int i = 0; i < opt.Nunits; i++)
	{
		//build based on options
		this->initialise_lung_unit(this->units[i], Vbag, Vbag, Vratios[i], param_dict);
	}
	this->build_airway_model(VD, VDSfrac, opt);

	//sync option
	if(opt.breath_model_type == SYNC_MODEL_CODE)  this->create_sync_vent_solver();
	else this->create_async_vent_solver(param_dict.at(DELAY_PARAM_NAME));
	this->washout_start_inflation.resize(this->Cinit.size());
	vector<double> Deltas;
	Deltas.resize(this->Cinit.size());
	double total_deltas = 0;
	for(size_t n = 1; n < this->Cinit.size(); n++)
	{
		std::stringstream ss;
		ss << FRC_TEST_MODIFIER_NAME << "_" << n;
		Deltas[n] = param_dict.at(ss.str().c_str());
		total_deltas += Deltas[n];
	}
	Deltas[0] = -total_deltas/(this->Cinit.size() + total_deltas);   
	//by definition, mean must be 0, so first test must cancel out the rest
	this->washout_start_inflation[0] = Deltas[0]*param_dict.at(FRC_PARAM_NAME);
	for(size_t n = 1; n < this->Cinit.size(); n++)
	{
		std::stringstream ss;
		ss << FRC_TEST_MODIFIER_NAME << "_" << n;
		this->washout_start_inflation[n] = (Deltas[0] + (1 + Deltas[0])*Deltas[n])*param_dict.at(FRC_PARAM_NAME);
	}

	//assume step is relative to DS volume
	this->params = params_h;
	this->simulate_washin = opt.simulate_washin;
	this->washout_start_timepoint = opt.washout_start_timepoint;
}

void MBWModelBase::run_washout_model(MBWModelOutputs* output)
{
	this->sim_conc.resize(this->sim_vol_steps.size());
	this->sim_igvol_cumul.resize(this->sim_vol_steps.size());
	for(size_t n = 0; n < this->sim_vol_steps.size(); n++) //loop over MBW tests
	{
		int npts = this->sim_vol_steps[n].size();
		this->sim_conc[n].resize(npts);
		this->sim_igvol_cumul[n].resize(npts);
		if(this->simulate_washin)
		{
			this->reset_model(0.0, this->washout_start_inflation[n]);
		}
		else
		{
			this->reset_model(this->Cinit[n], this->washout_start_inflation[n]);
		}
		for(int t = 0; t < npts; t++)    //loop over time points to simulate
		{
			/*if(t==0 || this->sim_vol_steps[n][t] *this->sim_vol_steps[n][t-1] < 0)
			{
				this->reassign_vent_ratios();
			}*/
			if(this->sim_vol_steps[n][t] > 0) //inhalation
			{
				if(t >= this->washout_start_timepoint)
				{
					this->set_inhaled_bc(0);   //concentration at mouth on inhalation
				}
				else
				{
					this->set_inhaled_bc(this->Cinit[n]);   //concentration at mouth on inhalation
				}
			}
			this->sim_igvol_cumul[n][t] = this->breath_step(this->sim_vol_steps[n][t], this->sim_step_durations->at(n)[t]);
			if(t > 0) this->sim_igvol_cumul[n][t] += this->sim_igvol_cumul[n][t-1];
			this->sim_conc[n][t] = this->get_mouth_conc();
		}
	}	
	this->measure_values(output);
}

void MBWModelBase::simulate(MBWModelOutputs* output)
{
	//run simulation
	//auto start = chrono::system_clock::now();

	//run MBW
	this->run_washout_model(output);
}

void MBWModelBase::create_sync_vent_solver()
{
	for(size_t i = 0; i < this->units.size(); i++)
	{
		this->units[i]->set_delay_ts(0.0);
	}
	this->vent_solver = std::make_shared<VentilationSolver>(&(this->units));
}

void MBWModelBase::create_async_vent_solver(const double & delay_ts)
{
	generate_async_delays(this->units, delay_ts);

	this->vent_solver = std::make_shared<AsyncVentilationSolver>(&(this->units), this->sim_vol_steps[0], 
		                                                         this->sim_step_durations->at(0));
}

void CompartmentalModelBase::build_airway_model(const double & VD, const double & VDSfrac,
												const MBWModelOptions & opt)
{
	using namespace std;
	vector<std::shared_ptr<FlexibleVolumeElement>> SDS(1);
	double VDS = this->mouth_point;
	double VDP = VD;
	if(VDSfrac > 0.0) //shared dead-space volume
	{
		VDS += VDSfrac*VD;
		VDP -= VDSfrac*VD;
		this->use_mouth_conc_generic = true;
	}
	else
	{
		this->use_mouth_conc_generic = false;
	}

	SDS[0] = std::make_shared<FlexibleVolumeElement>(VDS, 0, 0);
	this->shared_ds = std::make_shared<DSVolume>(SDS);

	//initialise model based on params
	this->private_ds.resize(opt.Nunits);
	double VDperunit = VDP / ((double) opt.Nunits);
	//build private dead-space objects
	for(int i = 0; i < opt.Nunits; i++)
	{
		vector<std::shared_ptr<FlexibleVolumeElement>> PDS(1);
		PDS[0] = std::make_shared<FlexibleVolumeElement>(VDperunit, 0, 0);
		this->private_ds[i] = std::make_shared<DSVolume>(PDS);
	}	

	//create mixing point
	double mp_vol_scale = opt.mix_vol_frac_step*(VDS + VDP);
	this->mixing_point = std::make_shared<MixingPoint>(mp_vol_scale);
}

double CompartmentalModelBase::get_mouth_conc()
{
	if(this->use_mouth_conc_generic)
	{
		return this->get_mouth_conc_generic();
	}
	else
	{
		return this->get_end_SDS_conc();
	}
}

void CompartmentalModelBase::measure_values(MBWModelOutputs* output)
{
	//count total number of data points
	int n_conc_measurement_points = 0, n_igvol_measurement_points = 0,
		tot_measurement_points = 0;
	for(size_t n = 0; n < this->conc_measurement_steps->size(); n++)
	{
		n_conc_measurement_points += int(this->conc_measurement_steps->at(n).size());
		n_igvol_measurement_points += int(this->igvol_measurement_steps->at(n).size());
	}
	tot_measurement_points = n_conc_measurement_points + n_igvol_measurement_points;

	//fill simulated vector
	output->simulated.resize(tot_measurement_points);
	int im_start = 0;
	for(size_t n = 0; n < this->conc_measurement_steps->size(); n++)
	{
		for(size_t im = 0; im < this->conc_measurement_steps->at(n).size(); im++)
		{
			int t = this->conc_measurement_steps->at(n)[im];
			output->simulated[im_start + im] = this->sim_conc[n][t];  //noise?
		}
		im_start += int(this->conc_measurement_steps->at(n).size());
		
		double last_cumul_volTO = 0;
		double last_cumul_igvol = 0;
		for(size_t im = 0; im < this->igvol_measurement_steps->at(n).size(); im++)
		{
			int t = this->igvol_measurement_steps->at(n)[im];

			output->simulated[im_start + im] = (this->sim_igvol_cumul[n][t] - last_cumul_igvol); //noise?
			last_cumul_igvol = this->sim_igvol_cumul[n][t];  
			last_cumul_volTO = this->sim_volTO_cumul[n][t];
		}
		im_start += int(this->igvol_measurement_steps->at(n).size());   
		
	}

	//fill Vrates vector
	output->Vrates.resize(NCOMPS);
	output->Vdelays.resize(NCOMPS);
	for(int iV = 0; iV < NCOMPS; iV++)
	{
		output->Vrates[iV] = this->units[iV]->get_vent_ratio();
		output->Vdelays[iV] = this->units[iV]->get_delay_ts();
	}
}

double CompartmentalModelBase::get_end_SDS_conc()
{
	return this->shared_ds->get_conc_bottom();
}

double CompartmentalModelBase::get_mouth_conc_generic()
{
	double vol = 0;
	size_t e = 0;
	while(vol < this->conc_measurement_point && e < this->shared_ds->VolumeElements.size())
	{
		vol += this->shared_ds->VolumeElements[e]->get_volume();
		e++;
	}
	if(vol < this->conc_measurement_point) 
	{
		std::cerr << "Warning: measurement not taken exactly at the mouth"
		          << std::endl;
	}

	double vover = vol - this->conc_measurement_point;
	if(vover > 0)
	{
		e--;   //go back to last edge
		double vnot_over = this->shared_ds->VolumeElements[e]->get_volume() - vover;
		double ctop_eff = (vover*this->shared_ds->VolumeElements[e]->get_ctop() 
						+ vnot_over*this->shared_ds->VolumeElements[e]->get_cbottom())
			            / (this->shared_ds->VolumeElements[e]->get_volume());

		return ctop_eff;
	}
	else
	{
		if(e > 0) return this->shared_ds->VolumeElements[e-1]->get_cbottom();
		else return this->shared_ds->VolumeElements[e]->get_ctop();
	}
}

void CompartmentalModelBase::reset_model(const double & C0, const double & inflation)
{
	//replace lung units
	for(size_t n = 0; n < this->units.size(); n++)
	{
			std::shared_ptr<LungUnit> old_unit = this->units[n];
			this->units[n]->reset(old_unit->get_FRC_volume(), old_unit->get_FRC_volume()*C0,
								  old_unit->get_vent_ratio());
	}

	//re-init ventilation model
	std::vector<double> initial_dvols(NCOMPS);   //initial vols not necessarily the same as FRC
	this->vent_solver->reset_solver(inflation, initial_dvols);
	for(size_t n = 0; n < this->units.size(); n++)
	{
		double vh = this->units[n]->get_FRC_volume() + initial_dvols[n];
		this->units[n]->set_volume(vh);
		this->units[n]->set_ig_volume(vh*C0);
	}

	//replace private ds
	for(size_t n = 0; n < this->private_ds.size(); n++)
	{
		double vol = this->private_ds[n]->get_volume();
		this->private_ds[n] = std::make_shared<DSVolume>(
			                       std::make_shared<FlexibleVolumeElement>(vol, C0, C0));
	}

	double vol = this->shared_ds->get_volume();
	//replace shared ds
	this->shared_ds = std::make_shared<DSVolume>(
		                   std::make_shared<FlexibleVolumeElement>(vol, C0, C0));
}

double CompartmentalModelBase::get_igvol_rebreathe()
{
	double vol = 0, igvol = 0;
	size_t e = 0;
	while(vol < this->conc_measurement_point && e < this->shared_ds->VolumeElements.size())
	{
		vol += this->shared_ds->VolumeElements[e]->get_volume();
		igvol += this->shared_ds->VolumeElements[e]->get_igvol();
		e++;
	}
	if(vol < this->conc_measurement_point) 
	{
		std::cerr << "Warning: measurement not taken exactly at the mouth"
		          << std::endl;
	}
	
	//mouth measurement taken somewhere in this volume
	double vover = vol - this->conc_measurement_point;
	if(vover > 0)
	{
		e--;   //go back to last edge
		double vnot_over = this->shared_ds->VolumeElements[e]->get_volume() - vover;
		double ctop_eff = (vover*this->shared_ds->VolumeElements[e]->get_ctop() 
						+ vnot_over*this->shared_ds->VolumeElements[e]->get_cbottom())
			            / (this->shared_ds->VolumeElements[e]->get_volume());

		FlexibleVolumeElement eff_vol(vover,ctop_eff,this->shared_ds->VolumeElements[e]->get_cbottom());
		igvol -= eff_vol.get_igvol();   //subtract ig vol from bit over
	}

	return igvol;
}

double CompartmentalModelBase::breath_step(const double & dvol, const double & dt)
{
	//vectors for mixing point computation
	std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> 
		          into_mixpoint_from_above, into_mixpoint_from_below;
	std::vector<double> dv_out_above, dv_out_below;
	double exhaled_igvol;
	into_mixpoint_from_above.reserve(1);
	into_mixpoint_from_below.reserve(this->units.size());
	dv_out_above.reserve(1);
	dv_out_below.reserve(this->units.size());
	exhaled_igvol  = -this->get_igvol_rebreathe();    //count ig volume in rebreathe vol initially 
	//(if exhaling, some or all of it will end up in exhaled vols so needs to be subtracted to avoid double counting
	//(if inhaling, some or all will end up in DS, so needs to be subtracted to avoid missing
	if(dvol > 0) //inhalation -- add vol to shared ds 
	{
		std::vector<std::shared_ptr<FlexibleVolumeElement>> injected(1);
		injected[0] = std::make_shared<FlexibleVolumeElement>(dvol, 
			                     this->inhaled_conc, this->inhaled_conc);
		exhaled_igvol -= injected[0]->get_igvol();  //rebreathed and inhaled
		std::vector<std::shared_ptr<FlexibleVolumeElement>> ejected;
		//inhale vol element into shared ds and store what comes out other end in above vector
		this->shared_ds->shunt_down(injected, ejected);
		into_mixpoint_from_above.push_back(ejected);
	}
	else //exhalation
	{
		dv_out_above.push_back(-dvol);
	}

	//get volume updates of units
	std::vector<double> dv;
	this->vent_solver->compute_volume_changes(dvol, dt, dv);

	//separate into inhaling and exhaling
	std::vector<size_t> dv_inh;
	dv_inh.reserve(this->units.size());
	for(size_t i = 0; i < dv.size(); i++)
	{
		if(dv[i] > 0) //inhaling
		{
			dv_inh.push_back(i);
			dv_out_below.push_back(dv[i]);
		}
		else  //exhaling
		{
			std::vector<std::shared_ptr<FlexibleVolumeElement>> injected,ejected;
			this->units[i]->exhale(-dv[i], injected, dt); 
			this->private_ds[i]->shunt_up(injected, ejected);
			into_mixpoint_from_below.push_back(ejected);
		}
	}

	//do mixing point calc
	std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> 
		                                out_of_mixpoint_above, out_of_mixpoint_below;
	this->mixing_point->mix_elements(into_mixpoint_from_above, 
		                             into_mixpoint_from_below, dv_out_above, 
		                             dv_out_below, out_of_mixpoint_above, 
									 out_of_mixpoint_below);

	//distribute vols out of mixing point
	if(dvol <= 0) //exhalation, first dvoutmp is shared ds
	{
		std::vector<std::shared_ptr<FlexibleVolumeElement>> exhaled;
		this->shared_ds->shunt_up(out_of_mixpoint_above[0], exhaled);
		for(size_t j = 0; j < exhaled.size(); j++) 
		{
			exhaled_igvol += exhaled[j]->get_igvol();
		}
	}
	//finish units inhaling
	for(size_t iu = 0; iu < dv_inh.size(); iu++)    
	{
		size_t i = dv_inh[iu];
		std::vector<std::shared_ptr<FlexibleVolumeElement>> ejected;
		this->private_ds[i]->shunt_down(out_of_mixpoint_below[iu], ejected);  
		this->units[i]->inhale(ejected,dt);
	}

	exhaled_igvol += this->get_igvol_rebreathe();   //add on igvol left in rebreathe  
	//(if exhaling, some or all of it will end up in exhaled vols so needs to be subtracted to avoid double counting
	//(if inhaling, some or all will end up in DS, so needs to be subtracted to avoid missing

	return exhaled_igvol;
}

void generate_vent_dist_lognormal_rand(const std::map<std::string,double> & params, 
								  const int & Ncomps, std::vector<double> & x_vals)
{
	x_vals.resize(Ncomps);
	double mu = -0.5*params.at(SIGMA_PARAM_NAME)*params.at(SIGMA_PARAM_NAME);
	boost::random::lognormal_distribution<> dist(mu, params.at(SIGMA_PARAM_NAME));
	for(int i = 0; i < Ncomps; i++)
	{
		x_vals[i] = dist(*(rng.get()));
		while(x_vals[i] > VENT_RATIO_LIMIT)
		{
			x_vals[i] = dist(*(rng.get()));
		}
	}

	rescale_to_unit_mean(x_vals);
}

void generate_vent_dist_lognormal_determin(const std::map<std::string,double> & params, 
								  const int & Ncomps, std::vector<double> & x_vals)
{
	double rt2sigma = sqrt(2)*params.at(SIGMA_PARAM_NAME);
	double sigmasq = params.at(SIGMA_PARAM_NAME)*params.at(SIGMA_PARAM_NAME);
	double cdf_up = 0, cdf_down = 0, cdf_mom_up = 0, cdf_mom_dwn = 0;
	for(int i = 0; i < Ncomps; i++)
	{
		if(i == Ncomps-1) //value of
		{
			cdf_up = 1.0;
			cdf_mom_up = 1.0;
		}
		else
		{
			//chop cdf into equal chunks
			cdf_up = double(i+1)/double(Ncomps);   //value of cdf at end of component
			//find corresponding x positions
			double xup = exp(-0.5*sigmasq - rt2sigma*boost::math::erf_inv(1.0 - 2*cdf_up));
			//calc cdf_moment at these points
			cdf_mom_up = 0.5*(1.0 - boost::math::erf((0.5*sigmasq - log(xup))/rt2sigma));	
		}
		//xval is cdf moment diff / cdf diff (= 1/Ncomps)
		double xval = (cdf_mom_up - cdf_mom_dwn)*Ncomps;
		cdf_down = cdf_up;
		cdf_mom_dwn = cdf_mom_up;
		x_vals[i] =  xval;
	}
	rescale_to_unit_mean(x_vals);
}

void generate_vent_dist_bimodal_rand(const std::map<std::string,double> & params, 
								  const int & Ncomps, std::vector<double> & x_vals)
{
	x_vals.resize(Ncomps);
	//slow compartment mean
	double s_mean = 1.0/(1.0 - params.at(VFASTFRAC_PARAM_NAME)*(1.0 - params.at(MURATIO_LS_PARAM_NAME)));
	//slow compartment sigma
	double s_sigma = params.at(SIGMA_PARAM_NAME) / params.at(SIGRATIO_LS_PARAM_NAME);
	//slow compartment mu
	double mu1 = log(s_mean) - 0.5*s_sigma*s_sigma;
	//fast compartment mu
	double mu2 = log(params.at(MURATIO_LS_PARAM_NAME)*s_mean) -
		         0.5*params.at(SIGMA_PARAM_NAME)*params.at(SIGMA_PARAM_NAME);

	boost::random::lognormal_distribution<> dist_slow(mu1, s_sigma);
	boost::random::lognormal_distribution<> dist_fast(mu2, params.at(SIGMA_PARAM_NAME));
	boost::random::uniform_01<> udist;

	for(int i = 0; i < Ncomps; i++)
	{
		double r = udist(*(rng.get()));
		if(r < params.at(VFASTFRAC_PARAM_NAME))
		{
			x_vals[i] = dist_fast(*(rng.get()));
			while(x_vals[i] > VENT_RATIO_LIMIT)
			{
				x_vals[i] = dist_fast(*(rng.get()));
			}
		}
		else
		{
			x_vals[i] = dist_slow(*(rng.get()));
			while(x_vals[i] > VENT_RATIO_LIMIT)
			{
				x_vals[i] = dist_slow(*(rng.get()));
			}
		}
	}

	rescale_to_unit_mean(x_vals);
}

void build_basic_lung_unit(std::shared_ptr<LungUnit> & unit, const double & V0, 
						   const double & IGV0, const double & DV, 
						   const std::map<std::string,double> & params)
{
	unit = std::make_shared<LungUnit>(V0, IGV0, DV);
}

void build_asymm_lung_unit(std::shared_ptr<LungUnit> & unit, const double & V0, 
						   const double & IGV0, const double & DV, 
						   const std::map<std::string,double> & params)
{
	unit = std::make_shared<AsymmLungUnit>(V0, IGV0, DV, params.at(ASYMM_PARAM_NAME), 
				                               params.at(DIFFSCALE_PARAM_NAME));
}

void rescale_to_unit_mean(std::vector<double> & x)
{
		//correct so that mean x = 1
	double mean_x = 0;
	int Ncomps = int(x.size());
	for(int i = 0; i < Ncomps; i++)
	{
		mean_x += x[i];
	}
	mean_x /= double(Ncomps);

	//double check_mean_x = 0;  //sanity check for debugging
	for(int i = 0; i < Ncomps; i++)
	{
		x[i] /= mean_x;
		//check_mean_x += x_vals[i];
	}
	//check_mean_x /= double(Ncomps);
}

void generate_async_delays(std::vector<std::shared_ptr<LungUnit>> & units, const double & delay_ts)
{
	double max_vr = 0;
	double min_vr = 1;
	for(size_t i = 0; i < NCOMPS; i++)
	{
		if(units[i]->get_vent_ratio() > max_vr) max_vr = units[i]->get_vent_ratio();
		if(units[i]->get_vent_ratio() < min_vr) min_vr = units[i]->get_vent_ratio();
	}

	for(size_t i = 0; i < NCOMPS; i++)
	{
		units[i]->set_delay_ts(delay_ts/(1 + units[i]->get_vent_ratio()));
		//delay ts is only achieved for vent ratio = 0
	}
}

void AsyncVentilationSolver::fit_warmup_data(const std::vector<double> & flow_data,
						      const std::vector<double> & flow_duration)
{
	//need to extract warm-up data from actual flow data
	//get mean exhalation size and duration
	double exh_vol_sum = 0, exh_dur_sum = 0;
	int exh_count = 0;
	//get mean inhalation size and duration
	double inh_vol_sum = 0, inh_dur_sum = 0;
	int inh_count = 0;
	
	for(size_t i = 0; i < flow_data.size(); i++)
	{
		if(flow_data[i] > 0)  //inhalation
		{
			if(i == 0 || flow_data[i-1] < 0)   //new inhalation
			{
				inh_count += 1;
			}
			inh_vol_sum += flow_data[i];
			inh_dur_sum += flow_duration[i];
		}
		else
		{
			if(i == 0 || flow_data[i-1] > 0)  //new exhalation   
			{
				exh_count += 1;
			}
			exh_vol_sum -= flow_data[i];
			exh_dur_sum += flow_duration[i];
		}
	}
	//assume it is just the mean volumes and durations to warm up
	double mean_exh_vol = exh_vol_sum/exh_count;
	double mean_inh_vol = inh_vol_sum/inh_count;
	double mean_exh_dur = exh_dur_sum/exh_count;
	double mean_inh_dur = inh_dur_sum/inh_count;
	double max_delay = this->delays.maxCoeff();

	int Nwarmup_breaths = int(ceil((mean_inh_dur + mean_exh_dur)/(10*max_delay)));
	this->warm_up_fluxes = std::vector<double>(2*Nwarmup_breaths);
	this->warm_up_flux_durations = std::vector<double>(2*Nwarmup_breaths);
	for(int i = 0; i < Nwarmup_breaths; i++)
	{
		this->warm_up_fluxes[2*i] = mean_inh_vol;
		this->warm_up_flux_durations[2*i] = mean_inh_dur;
		this->warm_up_fluxes[2*i+1] = -mean_exh_vol;
		this->warm_up_flux_durations[2*i+1] = mean_exh_dur;
	}
}

void AsyncVentilationSolver::run_warmup(const double & inflation_vol, std::vector<double> & initial_dvols)
{
	for(size_t i = 0; i < NCOMPS; i++)
	{
		this->Vols_old[i] = units->at(i)->get_volume();
		this->V0[i] = units->at(i)->get_FRC_volume();
	}
	this->Vp_old = this->Vols_old.sum()/NCOMPS;
	std::vector<double> wuf_copy = this->warm_up_fluxes;
	bool search_back = true;
	size_t nback = this->warm_up_fluxes.size()-1;
	if(inflation_vol != 0.0)
	{
		while(search_back)
		{
			if((this->warm_up_fluxes[nback] > 0 && inflation_vol > 0) ||
			   (this->warm_up_fluxes[nback] < 0 && inflation_vol < 0))
			{
				wuf_copy[nback] += inflation_vol;      //add extra inhalation/exhalation to end of warm up
				search_back = false;
			}
			nback--;
		}
	}

	for(size_t n = 0; n < wuf_copy.size(); n++)  //do warm up calc
	{
		this->compute_volume_changes(wuf_copy[n], this->warm_up_flux_durations[n], initial_dvols);
	}
}

AsyncVentilationSolver::AsyncVentilationSolver(std::vector<std::shared_ptr<LungUnit>> * units,
		                                      const std::vector<double> & flux_data,
					                          const std::vector<double> & flux_durations):VentilationSolver(units)
{
	this->dt_min = 0.01; 
	this->A.resize(NCOMPS+1,NCOMPS+1);
	this->Vols_old.resize(NCOMPS);
	this->V0.resize(NCOMPS);
	this->Vrates.resize(NCOMPS);
	this->delays.resize(NCOMPS);
	this->A.reserve(3*NCOMPS);
	for(size_t i = 0; i < NCOMPS; i++)
	{
		this->A.insert(i,i) = units->at(i)->get_delay_ts()/this->dt_min + 0.5;
		this->A.insert(NCOMPS,i) = 1.0/this->dt_min;
		this->A.insert(i,NCOMPS) = -0.5*units->at(i)->get_vent_ratio_original();
		this->Vrates[i] = units->at(i)->get_vent_ratio_original();
		this->delays[i] = units->at(i)->get_delay_ts();
	}
	// Compute the numerical factorization (A does not change)
	this->solver.analyzePattern(A); 
	this->solver.factorize(A);

	this->fit_warmup_data(flux_data,flux_durations);
}

void AsyncVentilationSolver::compute_volume_changes(const double & dvol, const double & dtime, 
						    std::vector<double> & dvols)
{
	if(dvols.size() != NCOMPS) dvols.resize(NCOMPS);
	int Nsteps = int(ceil(dtime/this->dt_min));
	double remainder = Nsteps*this->dt_min - dtime;
	Eigen::VectorXd b(NCOMPS+1);
	Eigen::VectorXd Vols_orig = this->Vols_old;
	Eigen::VectorXd Vols_new = this->Vols_old;
	double Vp_new = this->Vp_old;
	for(int n = 0; n < Nsteps; ++n)
	{
		this->Vols_old = Vols_new;
		this->Vp_old = Vp_new;
		b.head(NCOMPS) =  this->V0 + (1.0/this->dt_min)*this->delays.asDiagonal()*this->Vols_old 
						  - 0.5*this->Vols_old + 0.5*this->Vp_old*this->Vrates;
		b[NCOMPS] = this->Vols_old.sum()/this->dt_min + dvol/dtime;
		Eigen::VectorXd x = this->solver.solve(b);
		for(size_t i = 0; i < NCOMPS; i++)
		{
			Vols_new[i] = x[i];  //update volumes
		}
		Vp_new = x[NCOMPS];

	}
	//interpolate last step
	if(remainder > 0)
	{
		double dt_frac = remainder/this->dt_min;
		Vols_new = (1.0-dt_frac)*this->Vols_old + dt_frac*Vols_new;
		Vp_new = (1.0-dt_frac)*this->Vp_old + dt_frac*Vp_new;
	}
	//check no negative volumes and correct
	double vol_lost = 0.0;  
	for(size_t i = 0; i < NCOMPS; i++)
	{	
		if(this->Vols_old[i] < 0)
		{
			vol_lost -= Vols_new[i];
			Vols_new[i] = 0;
		}
	}
	if(vol_lost > 0)   //fudge correction if there are negatives
	{
		double vol_frac = vol_lost/Vols_new.sum();
		Vols_new = (1-vol_frac)*Vols_new;
		//std::cout << "Warning. Negative volumes encountered." << std::endl;
	}
	for(size_t i = 0; i < NCOMPS; i++) dvols[i] = Vols_new[i] - Vols_orig[i];
	this->Vols_old = Vols_new;
	this->Vp_old = Vp_new;
}

void TrumpetModelBase::build_airway_model(const double & VD, const double & VDSfrac,
										  const MBWModelOptions & opt)
{
	using namespace std;
	double VDPhys = VD - this->mouth_point;
	double VDMouth = TRUMPET_MOUTH_FRAC*VDPhys;
	double VDLung = (1-TRUMPET_MOUTH_FRAC)*VDPhys;

	//these can be defined in builder as they don't change
	int MaxCondGen = 15;
	double L2Dratio = 3.0;
	int NLmin = 20;
	//can we get max flow rate from the data? -- for now guess 0.5L/s
	double PeMax = 1.0;
	double LengthScaleFactor = pow(0.5,1.0/3.0);
	double MouthFrac = 0.2;
	//end
	
	//add in something to deal with mouth and DS volume here
	vector<shared_ptr<FlexibleVolumeElement>> ExtraDS(1);
	ExtraDS[0] = make_shared<FlexibleVolumeElement>(VD - VDLung, 0, 0);
	this->extra_ds = make_shared<DSVolume>(ExtraDS);

	double GenSep = VDSfrac*(MaxCondGen);
	double V0L = VDLung/(MaxCondGen+1);  //V per gen
	double R0dm = pow(3*V0L*(1-LengthScaleFactor)/(2*M_PI*L2Dratio*log(2.0)),1.0/3.0);
	//assume scherer diff in airways
	double dxapprox = 2*0.37*R0dm*PeMax;
	double Ltot = 2*R0dm*L2Dratio*(1.0-pow(LengthScaleFactor,MaxCondGen+1))/(1.0-LengthScaleFactor);
	double Lsep = 2*R0dm*L2Dratio*(1.0-pow(LengthScaleFactor,GenSep))/(1.0-LengthScaleFactor);
	int NL = max(NLmin,int(ceil(Ltot/dxapprox)));
	double dx = Ltot/double(NL);
	int NLsep = int(Lsep/dx);
	int NCondNodes = 1 + NLsep + (NL - NLsep)*NCOMPS;
	int NCondEdges = NLsep + (1 + NL - NLsep)*NCOMPS;
	std::vector<Eigen::Triplet<double>> IncTrips;
	this->EdgeCSA = Eigen::VectorXd::Zero(NCondEdges);
	this->NodeVol = Eigen::VectorXd::Zero(NCondNodes);
	IncTrips.reserve(2*NCondEdges);

	//
	double xh = 0;
	double Vint0 = 0;
	double invdx = 1.0/dx;
	double VolSF = log(2)/(3.0*(1 - LengthScaleFactor));
	double xSF = 0.5/(L2Dratio*R0dm);
	double V0 = 2*M_PI*R0dm*R0dm*R0dm*L2Dratio;
	for(int j = 0; j < NL; j++)
	{
		double z0 = 3*log(1.0-xh*xSF);  //replace with function
		double z1 = 3*log(1.0-(xh+0.5*dx)*xSF);
		double z2 = 3*log(1.0-(xh+dx)*xSF);
		double Vint1 = V0*VolSF*(z1-z0);  //from this node to +dx/2
		double Vint2 = V0*VolSF*(z2-z1);   //from +dx/2 to + dx
		double CSAh = invdx*(Vint1 + Vint2);
		double Nvolh = 0.5*(Vint0 + Vint1);
		if(j < NLsep)
		{
			IncTrips.push_back(Eigen::Triplet<double>(j,j, 1.0));
			IncTrips.push_back(Eigen::Triplet<double>(j,j+1,-1.0));
			this->EdgeCSA[j] = CSAh;
			this->NodeVol[j] = Nvolh;
		}
		else
		{
			if(j == NLsep)
			{
				int NodeIn = NLsep;
				this->NodeVol[NodeIn] = Nvolh;
				for(int n = 0; n < NCOMPS; n++)
				{
					int EdgeNo = NLsep + n;
					int NodeOut = NLsep + n + 1;
					IncTrips.push_back(Eigen::Triplet<double>(EdgeNo, NodeIn,  1.0));
					IncTrips.push_back(Eigen::Triplet<double>(EdgeNo, NodeOut,-1.0));
					this->EdgeCSA[EdgeNo] = CSAh/NCOMPS;
				}
			}
			else
			{
				for(int n = 0; n < NCOMPS; n++)
				{
					int EdgeNo = NLsep + (1 + j - NLsep)*n;
					int NodeIn = NLsep + (j - NLsep)*n + 1;
					int NodeOut = NLsep + (1 + j - NLsep)*n + 1;
					IncTrips.push_back(Eigen::Triplet<double>(EdgeNo, NodeIn,  1.0));
					IncTrips.push_back(Eigen::Triplet<double>(EdgeNo, NodeOut,-1.0));
					this->EdgeCSA[EdgeNo] = CSAh/NCOMPS;
					this->NodeVol[NodeIn] = Nvolh/NCOMPS;
					if(j==NL-1)
					{
						this->NodeVol[NodeOut] = Vint2;
					}
				}
			}
		}
		xh += dx;
		Vint0 = Vint2;
	}
	this->Incidence.setFromTriplets(IncTrips.begin(),IncTrips.end());
}


