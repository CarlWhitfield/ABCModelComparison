#include"MBW_model_params.h"

using namespace inlist;

MBWOptionList::MBWOptionList():OptionList()
{
	//mutliple choice options
	this->add(VDIST_KEY, std::make_shared<Option<char>>(DEFAULT_VDIST_TYPE, std::string(VDIST_KEY), 
		                       Vdist_option_list, Vdist_option_name_list, VDIST_OPTION_COUNT));
	this->add(LUNG_UNIT_KEY, std::make_shared<Option<char>>(DEFAULT_LUNG_UNIT, 
		std::string(LUNG_UNIT_KEY), Lung_unit_option_list, Lung_unit_option_name_list, 
		LUNG_UNIT_OPTION_COUNT));
}

MBWParameterList::MBWParameterList():ParameterList()
{
	//float or int options
	this->add(FRC_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_FRC, 
		                                            std::string(FRC_PARAM_NAME)));
	this->add(VD_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_VD, 
		                                            std::string(VD_PARAM_NAME)));
	this->add(VDSFRAC_PARAM_NAME, std::make_shared<FractionParam>(DEFAULT_VDSFRAC, 
		                                            std::string(VDSFRAC_PARAM_NAME)));
	this->add(SIGMA_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_SIGMA, 
		                                            std::string(SIGMA_PARAM_NAME)));
	this->add(MURATIO_LS_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_MURATIO, 
		                                            std::string(MURATIO_LS_PARAM_NAME)));
	this->add(SIGRATIO_LS_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_SIGRATIO, 
		                                            std::string(SIGRATIO_LS_PARAM_NAME)));
	this->add(VFASTFRAC_PARAM_NAME, std::make_shared<FractionParam>(DEFAULT_VFASTFRAC, 
		                                            std::string(VFASTFRAC_PARAM_NAME))); 
	this->add(ASYMM_PARAM_NAME, std::make_shared<FractionParam>(DEFAULT_ASYMM,   
		                                            std::string(ASYMM_PARAM_NAME)));   
	this->add(DIFFSCALE_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_DIFFSCALE,   
		                                            std::string(DIFFSCALE_PARAM_NAME)));
	this->add(MACHINE_DS_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_MACHINE_DS,
													std::string(MACHINE_DS_PARAM_NAME)));
	this->add(VT_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_VT,   
		                                            std::string(VT_PARAM_NAME)));
	this->add(VTFRANGE_PARAM_NAME, std::make_shared<FractionParam>(DEFAULT_VTFRANGE,
													std::string(VTFRANGE_PARAM_NAME)));
	this->add(BPERIOD_PARAM_NAME, std::make_shared<PositiveDoubleParam>(DEFAULT_BPERIOD,   
		                                            std::string(BPERIOD_PARAM_NAME)));

}

void MBWParameterList::check_validity(MBWOptionList *o)
{
	std::string generic_params[6] = {FRC_PARAM_NAME, VD_PARAM_NAME, SIGMA_PARAM_NAME, VDSFRAC_PARAM_NAME, 
		                             VT_PARAM_NAME, VTFRANGE_PARAM_NAME};
	double default_params[6] = {DEFAULT_FRC, DEFAULT_VD, DEFAULT_SIGMA, DEFAULT_VDSFRAC, DEFAULT_VT, 
		                             DEFAULT_VTFRANGE};
	for(int n = 0; n < 6; n++)
	{
		if(this->get_param<double>(generic_params[n])->isOK() == false)
		{
			std::cerr << generic_params[n] << " value invalid: setting to default\n";
			this->get_param<double>(generic_params[n])->set_phys_value(default_params[n]);
		}
	}
	//extra check for VD
	if(this->get_param<double>(VD_PARAM_NAME)->get_value() >= 
	   this->get_param<double>(FRC_PARAM_NAME)->get_value())
	{
		std::cerr << VD_PARAM_NAME << " value larger than FRC: setting to 5% of FRC\n";
		this->get_param<double>(VD_PARAM_NAME)->set_phys_value(0.05*
			             this->get_param<double>(FRC_PARAM_NAME)->get_value());
	}
	//check bimodal params if relevant
	if(o->get_option<char>(VDIST_KEY)->get_value() == BIMODAL_CODE)
	{
		std::string bimodal_params[3] = {MURATIO_LS_PARAM_NAME, SIGRATIO_LS_PARAM_NAME, VFASTFRAC_PARAM_NAME};
		double default_bm_params[3] = {DEFAULT_MURATIO, DEFAULT_SIGRATIO, DEFAULT_VFASTFRAC};
		for(int n = 0; n < 3; n++)
		{
			if(this->get_param<double>(bimodal_params[n])->isOK() == false)
			{
				std::cerr << bimodal_params[n] << " value invalid: setting to default\n";
				this->get_param<double>(bimodal_params[n])->set_phys_value(default_bm_params[n]);
			}
		}
	}
	//check asymm params if relevant
	if(o->get_option<char>(LUNG_UNIT_KEY)->get_value() == ASYMM_UNIT_CODE)
	{
		std::string asymm_params[2] = {ASYMM_PARAM_NAME, DIFFSCALE_PARAM_NAME};
		double default_as_params[2] = {DEFAULT_ASYMM, DEFAULT_DIFFSCALE};
		for(int n = 0; n < 2; n++)
		{
			if(this->get_param<double>(asymm_params[n])->isOK() == false)
			{
				std::cerr << asymm_params[n] << " value invalid: setting to default\n";
				this->get_param<double>(asymm_params[n])->set_phys_value(default_as_params[n]);
			}
		}
	}
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