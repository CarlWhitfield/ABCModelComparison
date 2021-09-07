#include"MBW_models.h"
#include"MBW_model_params.h"
#include"read_write_codes.h"
#include<boost\filesystem.hpp>
#include<boost\regex.hpp>
#include<file_manip.h>
#include<map>

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

	boost::filesystem::directory_iterator it_end; //yields beyond end iterator 
	for(boost::filesystem::directory_iterator it(filepath); it!=it_end; ++it)
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
	
	//read test data
	this->Ntests = int(MBWTestPaths.size());
	this->conc_measurements.resize(this->Ntests);
	this->measurement_steps.resize(this->Ntests);
	this->sim_vol_steps.resize(this->Ntests);
	this->sim_step_durations.resize(this->Ntests);
	this->dist_weights.resize(this->Ntests);
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
		for(int iline = 1; iline < int(infile.size()); iline++)
		{
			std::vector<std::string> line = string_split(infile[iline], ",");
			this->sim_vol_steps[n].push_back(StringToNumber<double>(line[col_numbers.at(VOLUME_NAME)]));
			this->av_vol_step += abs(this->sim_vol_steps[n][iline-1]);
			this->sim_step_durations[n].push_back(StringToNumber<double>(line[col_numbers.at(DURATION_NAME)]));
			if(StringToNumber<int>(line[col_numbers.at(MEASURED_NAME)]))
			{
				this->measurement_steps[n].push_back(iline-1);
				this->conc_measurements[n].push_back(StringToNumber<double>(line[col_numbers.at(CONC_NAME)]));
				this->dist_weights[n].push_back(1.0);  //automatically set to 1, can re-write to include weights
			}
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
	this->MRI_vbag = StringToNumber<double>(topline[col_numbers.at(VBAG_NAME)]);
	//our definition of FRC in this model excludes dead space
	this->FRC0 -= this->dead_space;

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
	std::shared_ptr<MBWOptionList> options = std::make_shared<MBWOptionList>();
	std::shared_ptr<MBWParameterList> params = std::make_shared<MBWParameterList>();

	std::vector<std::string> infile = get_all_lines(params_file);
	for(unsigned int i = 0; i < infile.size(); i++)
	{
		std::vector<std::string> line = string_split(infile[i]," ");
		
		if (!(options->parse(line[0],line[1]))) 
		{
			if(!(params->parse(line[0],line[1])))
			{
				std::cerr << "Error, did not recognise option " << line[0] << "\n";
			}
		}
	}
	params->check_validity(options.get());
	//now use these to create input file

	//generate flux -- can be moved to another function
	//first generate breath sizes
	double bp = params->get_param<double>(BPERIOD_PARAM_NAME)->get_value();
	double vt = params->get_param<double>(VT_PARAM_NAME)->get_value();
	double dvt = params->get_param<double>(VTFRANGE_PARAM_NAME)->get_value();

	int Nwashin = MBW_BMAX;
	int Nwashin_steps = int(Nwashin*bp*100);
	int Nwashout = MBW_BMAX;
	int Nwashout_steps = int(Nwashout*bp*100);
	this->Ntests=3;
	boost::random::uniform_01<double> u01;
	for(int n = 0; n < this->Ntests; n++)
	{
		std::vector<double> breath_vols(Nwashin+Nwashout, 0.0);
		for(int j = 0; j < Nwashin; j++)
		{
			breath_vols[j] = 2.0*(vt + 2*dvt*(u01(*(rng.get())) - 0.5)); 
			//washin breath size is double
		}
		for(int j = Nwashin; j < Nwashin + Nwashout; j++)
		{
			breath_vols[j] = vt + 2*dvt*(u01(*(rng.get())) - 0.5);
		}
		this->sim_step_durations.push_back(std::vector<double>(Nwashin_steps+Nwashout_steps, 0.01));
		this->sim_vol_steps.push_back(std::vector<double>(Nwashin_steps+Nwashout_steps, 0.0));
		this->measurement_steps.push_back(std::vector<int>(Nwashin_steps+Nwashout_steps, 0));
		for(int i = 0; i < Nwashin_steps+Nwashout_steps; i++)
		{
			int nb = int(i*0.01/bp);
			if(i*0.01/bp - nb < 0.5) //inhale
			{
				this->sim_vol_steps[n][i] = 0.01*breath_vols[nb]/(0.5*bp);
			}
			else //exhale
			{
				this->sim_vol_steps[n][i] = -0.01*breath_vols[nb]/(0.5*bp);
			}
			this->measurement_steps[n][i] = i;
		}
		this->Cinit.push_back(0.2);  //in %
		this->Cinit_std.push_back(SF6_NOISE);
	}
	this->machine_ds = params->get_param<double>(MACHINE_DS_PARAM_NAME)->get_value();
	this->av_vol_step = 2*0.01*vt/bp;
	//create and return options file
	opts.lung_unit_type = options->get_option<char>(LUNG_UNIT_KEY)->get_value();
	opts.vent_dist_type = options->get_option<char>(VDIST_KEY)->get_value();
	opts.NMR_samples = NMRI_SAMPLES;            //could be user inputted
	opts.Nunits = NCOMPS;                       //could be user inputted
	opts.mix_vol_frac_step = NMV_FRAC_STEP;     //could be user inputted
	opts.simulate_washin = true;
	opts.washout_start_timepoint = Nwashin_steps;


	//create params for MBW model
	MBWModelParams.push_back(params->get_param<double>(FRC_PARAM_NAME)->get_value());
	MBWModelParamNames.push_back(FRC_PARAM_NAME);
	MBWModelParams.push_back(params->get_param<double>(VD_PARAM_NAME)->get_value());
	MBWModelParamNames.push_back(VD_PARAM_NAME);
	MBWModelParams.push_back(params->get_param<double>(SIGMA_PARAM_NAME)->get_value());
	MBWModelParamNames.push_back(SIGMA_PARAM_NAME);
	bool SDS = (params->get_param<double>(VDSFRAC_PARAM_NAME)->get_value() > 0);
	if(SDS)
	{
		MBWModelParams.push_back(params->get_param<double>(VDSFRAC_PARAM_NAME)->get_value());
		MBWModelParamNames.push_back(VDSFRAC_PARAM_NAME);
	}
	if(opts.vent_dist_type == BIMODAL_CODE)
	{
		MBWModelParams.push_back(params->get_param<double>(MURATIO_LS_PARAM_NAME)->get_value());
		MBWModelParamNames.push_back(MURATIO_LS_PARAM_NAME);
		MBWModelParams.push_back(params->get_param<double>(SIGRATIO_LS_PARAM_NAME)->get_value());
		MBWModelParamNames.push_back(SIGRATIO_LS_PARAM_NAME);
		MBWModelParams.push_back(params->get_param<double>(VFASTFRAC_PARAM_NAME)->get_value());
		MBWModelParamNames.push_back(VFASTFRAC_PARAM_NAME);
	}
	if(opts.lung_unit_type == ASYMM_UNIT_CODE)
	{
		MBWModelParams.push_back(params->get_param<double>(ASYMM_PARAM_NAME)->get_value());
		MBWModelParamNames.push_back(ASYMM_PARAM_NAME);
		MBWModelParams.push_back(params->get_param<double>(DIFFSCALE_PARAM_NAME)->get_value());
		MBWModelParamNames.push_back(DIFFSCALE_PARAM_NAME);
	}
	//add in asymm unit params here
}

bool MBWModelOutputs::extra_outputs() const
{
	return true;
}

void MBWModelOutputs::print_extra_outputs(std::string & line) const
{
	std::stringstream ss;
	for(int imri = 0; imri < int(this->MRI_sample.size()); imri++)
	{
		if(imri > 0) ss << ",";
		ss << this->MRI_sample[imri];
	}
	std::string line1 = ss.str().c_str();
	ss.clear();
	ss.str("");
	for(int isim = 0; isim < int(this->simulated.size()); isim++)
	{
		ss << ",";
		ss << this->simulated[isim];
	}
	std::string line2 = ss.str().c_str();
	line = line1 + line2;
}

void MBWModelOutputs::get_headers(std::string & line) const
{
	std::stringstream ss;
	for(int imri = 0; imri < int(this->MRI_sample.size()); imri++)
	{
		if(imri > 0) ss << ",";
		ss << "MRI_sample_" << imri;
	}
	std::string line1 = ss.str().c_str();
	ss.clear();
	ss.str("");
	for(int isim = 0; isim < int(this->simulated.size()); isim++)
	{
		ss << ",";
		ss << "MBW_sample_" << isim;
	}
	std::string line2 = ss.str().c_str();
	line = line1 + line2;
}

MBWModelOptions::MBWModelOptions()   //default model options
{
	this->NMR_samples = NMRI_SAMPLES;
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

void BasicLognormalModelGenerator::initialise_params()
{
	this->param_names.resize(3);
	this->param_min.resize(3);
	this->param_max.resize(3);

	this->param_names[0] = FRC_PARAM_NAME;
	this->param_names[1] = VD_PARAM_NAME;
	this->param_names[2] = SIGMA_PARAM_NAME;

	this->param_min[0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[2] = SIGMA_MIN;

	this->param_max[0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[2] = SIGMA_MAX;
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
	m->build_model(opts, params, this->param_names);
}

void LognormalModelSDSGenerator::set_model_name()
{
	model_name = "Lognormal_SDS";
}

void LognormalModelSDSGenerator::initialise_params()
{
	this->param_names.resize(4);
	this->param_min.resize(4);
	this->param_max.resize(4);

	this->param_names[0] = FRC_PARAM_NAME;
	this->param_names[1] = VD_PARAM_NAME;
	this->param_names[2] = SIGMA_PARAM_NAME;
	this->param_names[3] = VDSFRAC_PARAM_NAME;

	this->param_min[0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[2] = SIGMA_MIN;
	this->param_min[3] = 0;

	this->param_max[0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[2] = SIGMA_MAX;
	this->param_max[3] = 1;
}

void BimodalModelSDSGenerator::set_model_name()
{
	model_name = "Bimodal_SDS";
}

void BimodalModelSDSGenerator::initialise_params()
{
	this->param_names.resize(7);
	this->param_min.resize(7);
	this->param_max.resize(7);

	this->param_names[0] = FRC_PARAM_NAME;
	this->param_names[1] = VD_PARAM_NAME;
	this->param_names[2] = SIGMA_PARAM_NAME;
	this->param_names[3] = VDSFRAC_PARAM_NAME;
	this->param_names[4] = MURATIO_LS_PARAM_NAME;
	this->param_names[5] = SIGRATIO_LS_PARAM_NAME;
	this->param_names[6] = VFASTFRAC_PARAM_NAME;


	this->param_min[0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[2] = SIGMA_MIN;
	this->param_min[3] = 0;
	this->param_min[4] = 1.0;   //ratio of large to small, cannot be < 1
	this->param_min[5] = 1.0/SIGRATIO_MAX;
	this->param_min[6] = 0.0;

	this->param_max[0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[2] = SIGMA_MAX;
	this->param_max[3] = 1;
	this->param_max[4] = MURATIO_MAX;   //ratio of large to small, cannot be < 1
	this->param_max[5] = SIGRATIO_MAX;
	this->param_max[6] = 1.0;
}

void BimodalModelSDSGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = BIMODAL_CODE;
	opts.lung_unit_type = BASIC_UNIT_CODE;
	m->build_model(opts, params, this->param_names);
}

void LognormalAsymmSDSModelGenerator::set_model_name()
{
	model_name = "Lognormal_SDS_Asymm";
}

void LognormalAsymmSDSModelGenerator::initialise_params()
{
	this->param_names.resize(6);
	this->param_min.resize(6);
	this->param_max.resize(6);

	this->param_names[0] = FRC_PARAM_NAME;
	this->param_names[1] = VD_PARAM_NAME;
	this->param_names[2] = SIGMA_PARAM_NAME;
	this->param_names[3] = VDSFRAC_PARAM_NAME;
	this->param_names[4] = ASYMM_PARAM_NAME;
	this->param_names[5] = DIFFSCALE_PARAM_NAME;

	this->param_min[0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[2] = SIGMA_MIN;
	this->param_min[3] = 0;
	this->param_min[4] = 0.0;
	this->param_min[5] = DIFFSCALE_MIN;

	this->param_max[0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[2] = SIGMA_MAX;
	this->param_max[3] = 1;
	this->param_max[4] = 1.0;
	this->param_max[5] = DIFFSCALE_MAX;
}

void LognormalAsymmSDSModelGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = LOGNORMAL_CODE;
	opts.lung_unit_type = ASYMM_UNIT_CODE;
	m->build_model(opts, params, this->param_names);
}

void BasicLognormalAsymmModelGenerator::set_model_name()
{
	model_name = "Lognormal_NoSDS_Asymm";
}

void BasicLognormalAsymmModelGenerator::initialise_params()
{
	this->param_names.resize(5);
	this->param_min.resize(5);
	this->param_max.resize(5);

	this->param_names[0] = FRC_PARAM_NAME;
	this->param_names[1] = VD_PARAM_NAME;
	this->param_names[2] = SIGMA_PARAM_NAME;
	this->param_names[3] = ASYMM_PARAM_NAME;
	this->param_names[4] = DIFFSCALE_PARAM_NAME;

	this->param_min[0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[2] = SIGMA_MIN;
	this->param_min[3] = 0.0;
	this->param_min[4] = DIFFSCALE_MIN;

	this->param_max[0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[2] = SIGMA_MAX;
	this->param_max[3] = 1.0;
	this->param_max[4] = DIFFSCALE_MAX;
}

void BasicLognormalAsymmModelGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = LOGNORMAL_CODE;
	opts.lung_unit_type = ASYMM_UNIT_CODE;
	m->build_model(opts, params, this->param_names);
}

void BimodalAsymmSDSModelGenerator::set_model_name()
{
	model_name = "Bimodal_SDS_Asymm";
}

void BimodalAsymmSDSModelGenerator::initialise_params()
{
	this->param_names.resize(9);
	this->param_min.resize(9);
	this->param_max.resize(9);

	this->param_names[0] = FRC_PARAM_NAME;
	this->param_names[1] = VD_PARAM_NAME;
	this->param_names[2] = SIGMA_PARAM_NAME;
	this->param_names[3] = VDSFRAC_PARAM_NAME;
	this->param_names[4] = MURATIO_LS_PARAM_NAME;
	this->param_names[5] = SIGRATIO_LS_PARAM_NAME;
	this->param_names[6] = VFASTFRAC_PARAM_NAME;
	this->param_names[7] = ASYMM_PARAM_NAME;
	this->param_names[8] = DIFFSCALE_PARAM_NAME;

	this->param_min[0] = FRC_MIN_FRAC*inputs->FRC0;
	this->param_min[1] = VD_MIN_FRAC*inputs->dead_space;
	this->param_min[2] = SIGMA_MIN;
	this->param_min[3] = 0;
	this->param_min[4] = 1.0;   //ratio of large to small, cannot be < 1
	this->param_min[5] = 1.0/SIGRATIO_MAX;
	this->param_min[6] = 0.0;
	this->param_min[7] = 0.0;
	this->param_min[8] = DIFFSCALE_MIN;

	this->param_max[0] = FRC_MAX_FRAC*inputs->FRC0;
	this->param_max[1] = VD_MAX_FRAC*inputs->dead_space;
	this->param_max[2] = SIGMA_MAX;
	this->param_max[3] = 1;
	this->param_max[4] = MURATIO_MAX;   //ratio of large to small, cannot be < 1
	this->param_max[5] = SIGRATIO_MAX;
	this->param_max[6] = 1.0;
	this->param_max[7] = 1.0;
	this->param_max[8] = DIFFSCALE_MAX;
}

void BimodalAsymmSDSModelGenerator::generate_model(const std::vector<double> & params,
		                                      std::shared_ptr<CompartmentalModelBase> & m) const
{
	m = std::make_shared<CompartmentalModelBase>();
	m->set_input_data(inputs);
	MBWModelOptions opts;
	opts.vent_dist_type = BIMODAL_CODE;
	opts.lung_unit_type = ASYMM_UNIT_CODE;
	m->build_model(opts, params, this->param_names);
}

void CompartmentalModelBase::set_input_data(const MBWModelInputs * inputs)
{
	//pointers to constant quantities first (no need to copy locally)
	this->sim_step_durations = &(inputs->sim_step_durations);
	this->measurement_steps = &(inputs->measurement_steps);


	boost::random::normal_distribution<> volndist(0.0, VOL_NOISE*inputs->av_vol_step);  //fractional noise
	this->Cinit.resize(inputs->Ntests);
	this->sim_vol_steps.resize(inputs->Ntests);
	for(int it = 0; it < inputs->Ntests; it++)
	{
		//perturb Cinit
		boost::random::normal_distribution<> Cinitndist(0, inputs->Cinit_std[it]);  //absolute noise
		this->Cinit[it] = inputs->Cinit[it] + Cinitndist(*(rng.get()));
		//perturb vol steps
		this->sim_vol_steps[it].resize(inputs->sim_vol_steps[it].size());
		for(int is = 0; is < int(inputs->sim_vol_steps[it].size()); is++)
		{
			this->sim_vol_steps[it][is] = inputs->sim_vol_steps[it][is] + volndist(*(rng.get()));
		}
	}
	this->mouth_point = inputs->machine_ds;
	this->MRI_bag_vol = inputs->MRI_vbag;
}

void CompartmentalModelBase::build_model(const MBWModelOptions & opt,
										 const std::vector<double> & params_h,
		                                 const std::vector<std::string> param_names_h)
{
	using namespace std;
	map<std::string, double> param_dict;
	for(int ip = 0; ip < int(params_h.size()); ip++)
	{
		param_dict[param_names_h[ip]] = params_h[ip];
	}

	vector<std::shared_ptr<FlexibleVolumeElement>> SDS(1);
	double VDS = this->mouth_point;
	double VDP = param_dict.at(VD_PARAM_NAME);
	if(param_dict.find(VDSFRAC_PARAM_NAME) != param_dict.end()) //shared dead-space volume
	{
		VDS += param_dict.at(VDSFRAC_PARAM_NAME)*param_dict.at(VD_PARAM_NAME);
		VDP -= param_dict.at(VDSFRAC_PARAM_NAME)*param_dict.at(VD_PARAM_NAME);
		this->get_mouth_conc = &CompartmentalModelBase::get_mouth_conc_generic;
	}
	else
	{
		this->get_mouth_conc = &CompartmentalModelBase::get_end_SDS_conc;
	}

	SDS[0] = std::make_shared<FlexibleVolumeElement>(VDS, 0, 0);
	this->shared_ds = std::make_shared<DSVolume>(SDS);

	//initialise model based on params
	this->private_ds.resize(opt.Nunits);
	this->units.resize(opt.Nunits);
	double VDperunit = VDP / ((double) opt.Nunits);
	//build private dead-space objects
	for(int i = 0; i < opt.Nunits; i++)
	{
		vector<std::shared_ptr<FlexibleVolumeElement>> PDS(1);
		PDS[0] = std::make_shared<FlexibleVolumeElement>(VDperunit, 0, 0);
		this->private_ds[i] = std::make_shared<DSVolume>(PDS);
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
	for(int i = 0; i < opt.Nunits; i++)
	{
		//build based on options
		this->initialise_lung_unit(this->units[i], Vbag, Vbag, Vratios[i], param_dict);
	}
	//create mixing point
	double mp_vol_scale = opt.mix_vol_frac_step*(VDS + VDP);
	//assume step is relative to DS volume
	this->mixing_point = std::make_shared<MixingPoint>(mp_vol_scale);
	this->params = params_h;
	this->NMRIPoints = opt.NMR_samples;
	this->washout_start_inflation = opt.washout_start_inflation;
	this->washout_start_inflation_duration = opt.washout_start_inflation_duration;
	this->simulate_washin = opt.simulate_washin;
	this->washout_start_timepoint = opt.washout_start_timepoint;
}

void CompartmentalModelBase::simulate(MBWModelOutputs* output)
{
	//run simulation
	//auto start = chrono::system_clock::now();

	//run MBW
	this->run_washout_model(output);
	//run MRI measurement
	this->run_MRI_model(output);
}

void CompartmentalModelBase::run_washout_model(MBWModelOutputs* output)
{
	this->sim_conc.resize(this->sim_vol_steps.size());
	for(size_t n = 0; n < this->sim_vol_steps.size(); n++) //loop over MBW tests
	{
		int npts = this->sim_vol_steps[n].size();
		this->sim_conc[n].resize(npts);
		if(this->simulate_washin)
		{
			this->reset_model(0.0, this->washout_start_inflation,
			              this->washout_start_inflation_duration);
		}
		else
		{
			this->reset_model(this->Cinit[n], this->washout_start_inflation,
			              this->washout_start_inflation_duration);
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
			this->breath_step(this->sim_vol_steps[n][t], this->sim_step_durations->at(n)[t]);
			this->sim_conc[n][t] = (this->*get_mouth_conc)();
		}
	}	
	this->measure_values(output);
}

void CompartmentalModelBase::measure_values(MBWModelOutputs* output)
{
	//count total number of data points
	int tot_measurement_points = 0;
	for(size_t n = 0; n < this->measurement_steps->size(); n++)
	{
		tot_measurement_points += int(this->measurement_steps->at(n).size());
	}
	//fill simulated vector
	output->simulated.resize(tot_measurement_points);
	int im_start = 0;
	for(size_t n = 0; n < this->measurement_steps->size(); n++)
	{
		for(size_t im = 0; im < this->measurement_steps->at(n).size(); im++)
		{
			int t = this->measurement_steps->at(n)[im];
			output->simulated[im_start + im] = this->sim_conc[n][t];
		}
		im_start += int(this->measurement_steps->at(n).size());
	}
}

void CompartmentalModelBase::run_MRI_model(MBWModelOutputs* output)
{
	this->reset_model(0);   //sets all concs to zero
	this->set_inhaled_bc(1.0);
	this->breath_step(this->MRI_bag_vol, 5.0);
	this->measure_MRI_dist(output);
}

void CompartmentalModelBase::measure_MRI_dist(MBWModelOutputs* output)
{
	//assume private DS is not measured
	//sample randomly from model
	double tot_gas_vol = 0, tot_ig_vol = 0;
	double non_gas_frac = 0.1;
	std::vector<double> cumulative_vol(this->units.size() + this->private_ds.size());
	for(size_t i = 0; i < this->units.size(); i++)
	{
		tot_gas_vol += this->units[i]->get_volume();
		tot_ig_vol += this->units[i]->get_conc()*this->units[i]->get_volume();
		cumulative_vol[i] = tot_gas_vol;
	}
	//for(size_t i = 0; i < this->private_ds.size(); i++)
	//{
	//	tot_gas_vol += this->private_ds[i]->get_volume();
	//	tot_ig_vol += this->private_ds[i]->sum_ig_volume();
	//	cumulative_vol[this->units.size() + i] = tot_gas_vol;
	//}

	boost::random::uniform_01<> udist;
	boost::random::normal_distribution<> ndist(0, MRI_NOISE_FRAC);
	output->MRI_sample.resize(this->NMRIPoints);
	for(int i = 0; i < this->NMRIPoints; i++)
	{
		double vol_point;
		double measure_error;
		if(MRI_NOISE_FRAC > 0)  //if random
		{
			vol_point = tot_gas_vol*udist(*(rng.get()));
			measure_error = (tot_ig_vol/tot_gas_vol)*ndist(*(rng.get()));
		}
		else   //if deterministic -- evenly spaced samples
		{
			vol_point = tot_gas_vol*(i+0.5)/this->NMRIPoints;
			measure_error = 0;
		}
		size_t n = 0;
		while(cumulative_vol[n] < vol_point)
		{
			n++;
		}
		if(n < this->units.size())
		{
			output->MRI_sample[i] = this->units[n]->get_ig_volume() 
				             / this->units[n]->get_volume() + measure_error;
		}

		if(output->MRI_sample[i] < 0) output->MRI_sample[i] = 0;  //no negative values allowed
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
	while(vol < this->mouth_point && e < this->shared_ds->VolumeElements.size())
	{
		vol += this->shared_ds->VolumeElements[e]->get_volume();
		e++;
	}
	if(vol < this->mouth_point) 
	{
		std::cerr << "Warning: measurement not taken exactly at the mouth"
		          << std::endl;
	}
	e--;
	//mouth measurement taken somewhere in this volume
	
	return ((vol - this->mouth_point) * this->shared_ds->VolumeElements[e]->get_ctop() 
		    + (this->mouth_point - vol + this->shared_ds->VolumeElements[e]->get_volume())
			* this->shared_ds->VolumeElements[e]->get_cbottom())
			/(this->shared_ds->VolumeElements[e]->get_volume());
}

void CompartmentalModelBase::compute_volume_changes(const double & dvol, const double & dtime, 
												    std::vector<double> & vols)
{
	//general function to calc volume updates assuming no interdependence of units
	vols.resize(this->units.size());
	double mean_dvol = dvol/this->units.size();
	for(size_t i = 0; i < this->units.size(); i++)
	{
		vols[i] = this->units[i]->return_volume_change(mean_dvol, dtime);
	}
}

void CompartmentalModelBase::reset_model(const double & C0, const double & inflation, const double & dt)
{
	//replace lung units
	std::vector<double> dv;
	this->compute_volume_changes(inflation, dt, dv);
	for(size_t n = 0; n < this->units.size(); n++)
	{
		std::shared_ptr<LungUnit> old_unit = this->units[n];
		this->units[n]->reset(old_unit->get_FRC_volume(), old_unit->get_FRC_volume()*C0,
			                  old_unit->get_vent_ratio());
		if(inflation > 0)
		{
			std::vector<std::shared_ptr<FlexibleVolumeElement>> inhaled;
			inhaled.push_back(std::make_shared<FlexibleVolumeElement>(dv[n],C0,C0));
			this->units[n]->inhale(inhaled, dt);
		}
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

void CompartmentalModelBase::breath_step(const double & dvol, const double & dt)
{
	//vectors for mixing point computation
	std::vector<std::vector<std::shared_ptr<FlexibleVolumeElement>>> 
		          into_mixpoint_from_above, into_mixpoint_from_below;
	std::vector<double> dv_out_above, dv_out_below;
	into_mixpoint_from_above.reserve(1);
	into_mixpoint_from_below.reserve(this->units.size());
	dv_out_above.reserve(1);
	dv_out_below.reserve(this->units.size());
	if(dvol > 0) //inhalation -- add vol to shared ds 
	{
		std::vector<std::shared_ptr<FlexibleVolumeElement>> injected(1);
		injected[0] = std::make_shared<FlexibleVolumeElement>(dvol, 
			                     this->inhaled_conc, this->inhaled_conc);
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
	this->compute_volume_changes(dvol, dt, dv);

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
	}
	//finish units inhaling
	for(size_t iu = 0; iu < dv_inh.size(); iu++)    
	{
		size_t i = dv_inh[iu];
		std::vector<std::shared_ptr<FlexibleVolumeElement>> ejected;
		this->private_ds[i]->shunt_down(out_of_mixpoint_below[iu], ejected);  
		this->units[i]->inhale(ejected,dt);
	}
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