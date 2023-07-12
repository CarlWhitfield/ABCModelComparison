#ifndef	MBW_MODEL_PARAMS_H
#define MBW_MODEL_PARAMS_H

#include"list_template.h"
#include "MBW_models.h"

#define DEFAULT_VDIST_TYPE LOGNORMAL_CODE
#define DEFAULT_LUNG_UNIT BASIC_UNIT_CODE
#define DEFAULT_BREATH_MODEL SYNC_MODEL_CODE
#define DEFAULT_FRC 3.0 //L
#define DEFAULT_VD 0.1  //L
#define DEFAULT_VDSFRAC 0.0   //0--1
#define DEFAULT_SIGMA 0.5     //Unitless
#define DEFAULT_MURATIO 1.0    //Unitless >= 1
#define DEFAULT_SIGRATIO 1.0    //Unitless > 0
#define DEFAULT_VFASTFRAC 1.0   //Unitless 0--1
#define DEFAULT_ASYMM 0.0     //0--1
#define DEFAULT_DIFFSCALE 0.5   //s
#define DEFAULT_DELAY_SECS 1.0  //s

#define MACHINE_DS_PARAM_NAME "Machine_deadspace"
#define DEFAULT_MACHINE_DS 0.036   //L
#define VT_PARAM_NAME "Tidal_Volume"
#define DEFAULT_VT 1.0   //L
#define VTFRANGE_PARAM_NAME "Tidal_Volume_Frac_Range"
#define DEFAULT_VTFRANGE 0.2         //fractional range of variation (i.e. U(0.8,1.2)*VT)
#define BPERIOD_PARAM_NAME "Breath_period"
#define DEFAULT_BPERIOD 5.0     //s
#define MBW_BMAX 30

class MBWOptionList: public inlist::OptionList<char, bool>
{
public:
	MBWOptionList();
};

class MBWParameterList: public inlist::ParameterList<int, double>
{
public:
	MBWParameterList();   //constructor for simulation
	void add_start_inflation_params(const int & Ntests);
	void check_validity(MBWOptionList *o);
};

class DoubleParam: public inlist::Parameter<double>
{
public:
	DoubleParam(const double & val, const std::string & nam):Parameter<double>(val, nam){};
	bool isOK();
};

class PositiveDoubleParam: public inlist::Parameter<double>
{
public:
	PositiveDoubleParam(const double & val, const std::string & nam):Parameter<double>(val, nam){};
	bool isOK();
};

class FractionParam: public inlist::Parameter<double>
{
public:
	FractionParam(const double & val, const std::string & nam):Parameter<double>(val, nam){};
	bool isOK();
};

#endif