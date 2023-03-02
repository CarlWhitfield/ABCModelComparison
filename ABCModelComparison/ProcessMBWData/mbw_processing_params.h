#ifndef	MBW_PROCESS_PARAMS_H
#define MBW_PROCESS_PARAMS_H

#include"list_template.h"

//Option definitions
#define CUTOFF_KEY "Cutoff_mode"
#define CUTOFF_OPTION_COUNT 2
#define LCI_CUTOFF_CODE 'l'    //LCI threshold
#define NBREATHS_CUTOFF_CODE 'b'   // no. of breaths
const char Cutoff_option_list[] = {LCI_CUTOFF_CODE,NBREATHS_CUTOFF_CODE};
const std::string Cutoff_option_name_list[] = {"LCI", "NBreaths"};

//bool options
#define CORRECT_FOR_BIAS_KEY "Bias_correction"
#define MEASURE_INSPIRED_KEY "Measure_inspired"
#define MEASURE_SUBSET_KEY "Measure_subset"
#define FIXED_VOLSTEP_KEY "Fixed_volstep"
#define TIDAL_WASHIN_KEY "Tidal_washin"

//Param definitions
#define TEMPERATURE_PARAM_NAME "Temperature"
#define HUMIDITY_PARAM_NAME "Humidity"
#define PRESSURE_PARAM_NAME "Pressure"
#define CO2_DELAY_PARAM_NAME "CO2Delay"
#define O2_DELAY_PARAM_NAME "CO2Delay"
#define MACHINE_DS_PARAM_NAME "MachineDS"
#define REBREATHE_VOL_PARAM_NAME "RBVol"
#define MRI_VBAG_PARAM_NAME "Vbag"
#define LCI_CUTOFF_PARAM_NAME "LCI_cutoff"
#define WASHIN_MIN_CONC_NAME "Washin_min_conc"   
#define WASHOUT_MAX_CONC_NAME "Washout_max_conc"  
#define NMEASURE_PARAM_NAME "NMeasure"
#define NBREATH_CUTOFF_NAME "NBreaths"
#define NPHASEII_PARAM_NAME "NPhaseII"   // 4 //in measured breath, number of phase II points to use
#define NPHASEIII_PARAM_NAME "NPhaseIII"  // 4 //in measured breath, number of phase III points to use
#define SIM_STEP_FRAC_NAME "Vol_step"   //if fixed step size true, simulation step size in frac of VT mean//change depending on whether rebreathe or tidal washin 

//defaults
#define DEFAULT_CUTOFF LCI_CUTOFF_CODE
#define DEFAULT_BIAS_CORRECTION true
#define DEFAULT_MEASURE_INSPIRED true
#define DEFAULT_MEASURE_SUBSET false   //if false all breaths are used
#define DEFAULT_FIXED_STEP_SIZE false   //use fixed vol step size in sims?
#define DEFAULT_TIDAL_WASHIN false   //use washin too for breath volume correction
#define DEFAULT_TEMP 20.0   //Celcius
#define DEFAULT_HUMIDITY 50.0   //Relative %
#define DEFAULT_PRESSURE 760.0   //Atmospheric mmHg
#define DEFAULT_CO2_DELAY 1000.0  //milliseconds
#define DEFAULT_O2_DELAY 1000.0  //milliseconds
#define DEFAULT_MACHINE_DS 48.0  //millilitres
#define DEFAULT_REBREATHE_VOL 65.0  //millilitres 
#define DEFAULT_MRI_VBAG 1.0 //litres
#define DEFAULT_LCI_CUTOFF 0.025   //fraction of Cinit where LCI point is defined
#define DEFAULT_WASHIN_MIN_CONC 0.15    //threshold conc for observations dureing inhalation above which determines washout has definitely begun
#define DEFAULT_WASHIN_MAX_CONC 0.01    //equivalent threshold below which determines washout has begun 
#define DEFAULT_NMEASURE 20          //if measure subset is true, (max) number of breaths to use in measurements
#define DEFAULT_NBREATH_CUTOFF 6   //if cutoff option is 'b', number of breaths to stop at
#define DEFAULT_NPHASEII 4   // 4 //in measured breath, number of phase II points to use
#define DEFAULT_NPHASEIII 4  // 4 //in measured breath, number of phase III points to use
#define DEFAULT_SIM_STEP_FRAC 0.04   //if fixed step size true, simulation step size in frac of VT mean//change depending on whether rebreathe or tidal washin 

class LCIOptionList: public inlist::OptionList<char, bool>
{
public:
	LCIOptionList();
};

class LCIParameterList: public inlist::ParameterList<int, double>
{
protected:
	int Ntests;
public:
	LCIParameterList(const int & Ntests);   //constructor for simulation
	void check_validity(LCIOptionList *o);
	inline void set_Ntests(const int & nt){ this->Ntests = nt; }
	inline int get_Ntests() const { return this->Ntests; }
};

class IntParam: public inlist::Parameter<int>
{
public:
	IntParam(const int & val, const std::string & nam):Parameter<int>(val, nam){};
	bool isOK();
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




