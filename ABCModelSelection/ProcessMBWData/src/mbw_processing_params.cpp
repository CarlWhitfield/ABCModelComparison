#include"read_mbw_data.h"
#include"mbw_processing_params.h"

using namespace::inlist;

LCIOptionList::LCIOptionList(){  //constructor for options list
	//mutliple choice options
	this->add(CUTOFF_KEY, std::make_shared<Option<char>>(DEFAULT_CUTOFF, std::string(CUTOFF_KEY), 
		                       Cutoff_option_list, Cutoff_option_name_list, CUTOFF_OPTION_COUNT));

	//bool options
	this->add(CORRECT_FOR_BIAS_KEY, std::make_shared<Option<bool>>(DEFAULT_BIAS_CORRECTION, 
		std::string(CORRECT_FOR_BIAS_KEY)));

	this->add(MEASURE_SUBSET_KEY, std::make_shared<Option<bool>>(DEFAULT_MEASURE_SUBSET, 
		std::string(MEASURE_SUBSET_KEY)));

	this->add(FIXED_VOLSTEP_KEY, std::make_shared<Option<bool>>(DEFAULT_FIXED_STEP_SIZE, 
		std::string(FIXED_VOLSTEP_KEY)));

	this->add(TIDAL_WASHIN_KEY, std::make_shared<Option<bool>>(DEFAULT_TIDAL_WASHIN, 
		std::string(TIDAL_WASHIN_KEY)));

	this->add(MEASURE_INSPIRED_KEY, std::make_shared<Option<bool>>(DEFAULT_MEASURE_INSPIRED,
		std::string(MEASURE_INSPIRED_KEY)));
}

LCIParameterList::LCIParameterList(const int & Ntests){  //constructor for options list
	this->set_Ntests(Ntests);

	std::stringstream ss;
	for(int i = 0; i < this->Ntests; i++){     //add params that can change between tests
		ss << TEMPERATURE_PARAM_NAME << "_" << i+1;
		this->add(ss.str().c_str(), std::make_shared<DoubleParam>(DEFAULT_TEMP, 
		std::string(ss.str().c_str())));
		ss.str("");
		ss.clear();
		
		ss << HUMIDITY_PARAM_NAME << "_" << i+1;
		this->add(ss.str().c_str(), std::make_shared<PositiveDoubleParam>(DEFAULT_HUMIDITY,
		std::string(ss.str().c_str())));
		ss.str("");
		ss.clear();

		ss << PRESSURE_PARAM_NAME << "_" << i+1;
		this->add(ss.str().c_str(), std::make_shared<PositiveDoubleParam>(DEFAULT_PRESSURE,
		std::string(ss.str().c_str())));
		ss.str("");
		ss.clear();

		ss << CO2_DELAY_PARAM_NAME << "_" << i+1;
		this->add(ss.str().c_str(), std::make_shared<PositiveDoubleParam>(DEFAULT_CO2_DELAY,
		std::string(ss.str().c_str())));
		ss.str("");
		ss.clear();

		ss << O2_DELAY_PARAM_NAME << "_" << i+1;
		this->add(ss.str().c_str(), std::make_shared<PositiveDoubleParam>(DEFAULT_O2_DELAY,
		std::string(ss.str().c_str())));
		ss.str("");
		ss.clear();		

		ss << CINIT_PARAM_NAME << "_" << i+1;
		this->add(ss.str().c_str(), std::make_shared<PositiveDoubleParam>(DEFAULT_CINIT,
		std::string(ss.str().c_str())));
		ss.str("");
		ss.clear();		
	}

	//add params that cannot change between tests
	this->add(MACHINE_DS_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_MACHINE_DS,
		std::string(MACHINE_DS_PARAM_NAME)));

	this->add(REBREATHE_VOL_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_REBREATHE_VOL,
		std::string(REBREATHE_VOL_PARAM_NAME)));

	this->add(MRI_VBAG_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_MRI_VBAG,
		std::string(MRI_VBAG_PARAM_NAME)));

	this->add(LCI_CUTOFF_PARAM_NAME, std::make_shared<FractionParam>(DEFAULT_LCI_CUTOFF,
		std::string(LCI_CUTOFF_PARAM_NAME)));

	this->add(WASHIN_MIN_CONC_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_WASHIN_MIN_CONC,
		std::string(WASHIN_MIN_CONC_NAME)));

	this->add(WASHOUT_MAX_CONC_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_WASHIN_MAX_CONC,
		std::string(WASHOUT_MAX_CONC_NAME)));

	this->add(NMEASURE_PARAM_NAME,  std::make_shared<IntParam>(DEFAULT_NMEASURE,
		std::string(NMEASURE_PARAM_NAME)));

	this->add(NBREATH_CUTOFF_NAME, std::make_shared<IntParam>(DEFAULT_NBREATH_CUTOFF,
		std::string(NBREATH_CUTOFF_NAME)));

	this->add(NPHASEII_PARAM_NAME, std::make_shared<IntParam>(DEFAULT_NPHASEII,
		std::string(NPHASEII_PARAM_NAME)));

	this->add(NPHASEIII_PARAM_NAME, std::make_shared<IntParam>(DEFAULT_NPHASEIII,
		std::string(NPHASEIII_PARAM_NAME)));

	this->add(SIM_STEP_FRAC_NAME, std::make_shared<FractionParam>(DEFAULT_SIM_STEP_FRAC,
		std::string(SIM_STEP_FRAC_NAME)));
}

void LCIParameterList::check_validity(LCIOptionList *o)
{
	for(auto & it: this->dict1)
	{
		if(!(it.second->isOK()))
		{
			std::cout << "Setting " << it.first << " to default value." << std::endl;
			it.second->set_to_default();
		}
	}

	for(auto & it: this->dict2)
	{
		if(!(it.second->isOK()))
		{
			std::cout << "Setting " << it.first << " to default value." << std::endl;
			it.second->set_to_default();
		}
	}

	//check for ones that have not been changed from default...
	std::stringstream ss;
	for(int n = 0; n < this->Ntests; n++){
		ss << CO2_DELAY_PARAM_NAME << "_" << n;  //CO2 delay is needed, and therefore important
		std::string CO2_Delay = ss.str().c_str();
		ss.str("");
		ss.clear();
		if(!(this->get_param<double>(CO2_Delay)->has_changed())){   //if not changed from default
			std::cerr << "Warning: CO2 Delay in test " << n << " not changed from default." << std::endl;
		}
		
		ss << O2_DELAY_PARAM_NAME << "_" << n;  //O2 delay is not needed generally, but should be set to CO2 delay if not given
		std::string O2_Delay = ss.str().c_str();
		ss.str("");
		ss.clear();
		if(!(this->get_param<double>(O2_Delay)->has_changed())){   //if not changed from default
			std::cerr << "Setting O2 Delay in test " << n << " to same as CO2 delay." << std::endl;
			this->get_param<double>(O2_Delay)->update_value(this->get_param<double>(CO2_Delay)->get_value());
		}
	}

	if(!(this->get_param<double>(MACHINE_DS_PARAM_NAME)->has_changed())){
		std::cerr << "Warning: Machine DS not changed from default." << std::endl;
	}
	
	if(!(this->get_param<double>(REBREATHE_VOL_PARAM_NAME)->has_changed())){
		std::cerr << "Warning: Rebreathe volume not changed from default." << std::endl;
	}
	

}

bool DoubleParam::isOK()
{
	return true;
}

bool PositiveDoubleParam::isOK()
{
	if(this->value >= 0.0) return true;
	else return false;
}

bool FractionParam::isOK()
{
	if(this->value < 0.0 || this->value > 1.0) return false;
	else return true;
}