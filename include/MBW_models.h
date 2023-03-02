#ifndef MBW_MODELS_H
#define MBW_MODELS_H

#define _USE_MATH_DEFINES
#include"ABC_model_selection.h"
#include"compartmental.h"
#include<Eigen/Sparse>
#include<limits>
#include<map>
#include<math.h>

//Define which data to fit to
#define USE_SF6_CONC true
#define USE_SF6_VOL false

//quantities that are neither model nor patient specific
const double SF6_NOISE = 0.0002;   //largest estimate from SNR in Horsley et al. Thorax 2007
const double VOL_NOISE = 0.01;    //SNR for **volume** from same paper, note this probably depends on flow rate
const double MRI_NOISE_FRAC = 0.02;    //1/SNR expected for MRI for mean -- ballpark figure
const double VENT_RATIO_LIMIT = 10;   //maximum ventilation ratio (i.e. truncation of ventilation distribution)

//default model options
#define NCOMPS 50  //
#define NMRI_SAMPLES 1000
#define NMV_FRAC_STEP 0.2

//prior limits
const double FRC_MIN_FRAC = 0.5;   //as a fraction of MBW measured FRC
const double FRC_MAX_FRAC = 2.0;
const double VD_MIN_FRAC = 0.5;   //as a fraction of MBW measured FDS
const double VD_MAX_FRAC = 3.0;
const double SIGMA_MIN = 0;   //as a fraction of MBW measured FDS
const double SIGMA_MAX = 2.0;
const double MURATIO_MAX = 20.0;  //increase if posterior suggests higher
const double SIGRATIO_MAX = 20.0;
const double DIFFSCALE_MIN = 1E-02;   //min time resolution
const double DIFFSCALE_MAX = 20;    //largest plausible diffusion timescale (scale of 2 breaths)
const double DELAY_MAX_SECS = 5.0;
const double FRC_TEST_MODIFIER_MIN = -0.25;
const double FRC_TEST_MODIFIER_MAX = 0.25;

#define FRC_PARAM_NAME "VAcin"
#define VD_PARAM_NAME "VD"
#define VDSFRAC_PARAM_NAME "VD_shared_fraction"
#define SIGMA_PARAM_NAME "Sigma"
#define MURATIO_LS_PARAM_NAME "Mu_ratio"
#define SIGRATIO_LS_PARAM_NAME "Sigma_ratio"
#define VFASTFRAC_PARAM_NAME "V_fast_frac"
#define ASYMM_PARAM_NAME "Acin_Asymm"
#define DIFFSCALE_PARAM_NAME "Acin_diff_tscale"
#define DELAY_PARAM_NAME "Delay_secs"
#define FRC_TEST_MODIFIER_NAME "FRC_diff_test"

//option codes for tree structure
#define VDIST_KEY "VentDist"
#define VDIST_OPTION_COUNT 2
#define LOGNORMAL_CODE 'l'
#define BIMODAL_CODE 'b'
const char Vdist_option_list[] = {LOGNORMAL_CODE,BIMODAL_CODE};
const std::string Vdist_option_name_list[] = {"Lognormal", "Bimodal"};

#define LUNG_UNIT_KEY "Unit"
#define LUNG_UNIT_OPTION_COUNT 2
#define BASIC_UNIT_CODE 'u'
#define ASYMM_UNIT_CODE 'a'
const char Lung_unit_option_list[] = {BASIC_UNIT_CODE,ASYMM_UNIT_CODE};
const std::string Lung_unit_option_name_list[] = {"Symm_unit", "Asymm_unit"};

#define BREATHING_MODEL_KEY "BreathModel"
#define BREATHING_MODEL_OPTION_COUNT 2
#define SYNC_MODEL_CODE 'e'
#define ASYNC_MODEL_CODE 'r'
const char Breathing_model_option_list[] = {SYNC_MODEL_CODE,ASYNC_MODEL_CODE};
const std::string Breathing_model_option_name_list[] = {"Sync_breathing_model", "Async_breathing_model"};

class CompartmentalModelBase;

class RMSEDistance: public DistanceFunctionBase{};

struct MBWModelOptions
{
public:
	//model specific options
	char vent_dist_type, lung_unit_type, breath_model_type;
	int Nunits, NMR_samples;
	double mix_vol_frac_step;
	double washout_start_inflation, washout_start_inflation_duration;
	bool simulate_washin;
	int washout_start_timepoint;

	MBWModelOptions();
};

class MBWModelInputs: public ModelInputBase
{
public:
	std::vector<std::vector<double>> sim_vol_steps, sim_step_durations,
		conc_measurements, igvol_diff_measurements, conc_weights, igvol_weights;
	std::vector<std::vector<int>> conc_measurement_steps, igvol_measurement_steps;
	std::vector<double> Cinit, Cinit_std;
	double dead_space, FRC0, machine_ds, MRI_vbag, rebreathe_vol; 
	double av_vol_step;
	//timestep and median fowler dead-space
	int Ntests;
	//Eigen::MatrixXd weights;
	std::string subject_name;

	MBWModelInputs(){};
	MBWModelInputs(const std::string & filepath):ModelInputBase(filepath)
	{
		this->read_inputs(filepath);
	}
	void read_inputs(const std::string & filepath);
	void generate_inputs(const std::string & params_file, MBWModelOptions & opts, 
						 std::vector<double> & MBWModelParams,  std::vector<std::string> & MBWModelParamNames);
};

class MBWModelOutputs: public ModelOutputBase
{
public:
	std::vector<double> MRI_sample;
	bool extra_outputs() const;
	void print_extra_outputs(std::string & line) const;
	void get_headers(std::string & line) const;
};

class AsymmLungUnit: public LungUnit
{
protected:
	double A, DT;
	double vol1, vol2, conc1, conc2, IGvol1, IGvol2;
public:
	AsymmLungUnit():LungUnit(){};
	AsymmLungUnit(const double & frc, const double & igvol, const double & x, 
		          const double & Asymm, const double & Dscale):LungUnit(frc, igvol, x)
	{
		this->A = Asymm;
		this->DT = Dscale;
		vol1 = 0.5*(1-A)*frc;
		vol2 = 0.5*(1+A)*frc;
		IGvol1 = 0.5*(1-A)*igvol;
		IGvol2 = 0.5*(1+A)*igvol;
		conc1 = IGvol1/vol1;
		conc2 = IGvol2/vol2;
	}   //lung unit constructor

	void inhale(const std::vector<std::shared_ptr<FlexibleVolumeElement>> & fv, 
		        const double & dt);   //absorb and destroy volume element
	void exhale(const double & dv, 
		     std::vector<std::shared_ptr<FlexibleVolumeElement>> & exhaled, 
			 const double & dt);

	inline void reset(const double & frc, const double & igvol, const double & x)
	{
		*(this) = AsymmLungUnit(frc,igvol,x,A,DT);  //call copy constructor
	}
};

class VentilationSolver
{
protected:
	std::vector<std::shared_ptr<LungUnit>> *units;
public:
	VentilationSolver(std::vector<std::shared_ptr<LungUnit>> * units)
	{
		this->units = units;
	}

	virtual void reset_solver(const double & inflation_vol, std::vector<double> & initial_dvols)
	{
		this->compute_volume_changes(inflation_vol, 0.0, initial_dvols);
	}

	virtual void compute_volume_changes(const double & dvol, const double & dtime, 
						            std::vector<double> & vols)
	{
		//general function to calc volume updates assuming no interdependence of units
		//pass units to function
		if(vols.size() != NCOMPS) vols.resize(NCOMPS);
		double mean_dvol = dvol/NCOMPS;
		for(size_t n = 0; n < NCOMPS; n++)
		{
			vols[n] = units->at(n)->return_volume_change(mean_dvol, dtime);
		}
	}
};

class AsyncVentilationSolver: public VentilationSolver
{
protected:
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd Vols_old, Vrates, V0, delays;
	double Vp_old, dt_min;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	std::vector<double> warm_up_fluxes, warm_up_flux_durations;
	void fit_warmup_data(const std::vector<double> & flow_data,
						      const std::vector<double> & flow_duration);
	void run_warmup(const double & inflation_vol, std::vector<double> & initial_dvols);
public:
	AsyncVentilationSolver(std::vector<std::shared_ptr<LungUnit>> * units,
		                   const std::vector<double> & flow_data,
						   const std::vector<double> & flow_duration);
	void reset_solver(const double & inflation_vol, std::vector<double> & initial_dvols)
	{
		this->run_warmup(inflation_vol, initial_dvols);
	}

	void compute_volume_changes(const double & dvol, const double & dtime, 
						    std::vector<double> & vols);
};

class MBWModelBase: public ModelBase<MBWModelOutputs>
{
protected:
	std::vector<double> washout_start_inflation;
	bool simulate_washin;
	int washout_start_timepoint;
	double mouth_point, conc_measurement_point, inhaled_conc;

	//pointers to data from MBW input and options (stored in generator)
	//non-pointers to data that is perturbed before use
	std::vector<std::vector<double>> sim_conc, sim_igvol_cumul, sim_vol_steps, sim_volTO_cumul;
	std::vector<double> Cinit;
	const std::vector<std::vector<double>> *sim_step_durations;
	const std::vector<std::vector<int>> *conc_measurement_steps, *igvol_measurement_steps;
	double MRI_bag_vol;
	int NMRIPoints;

	void run_washout_model(MBWModelOutputs* output);
	void run_MRI_model(MBWModelOutputs* output);

	//this is a dummy model, so some of these virtual functions are not defined
	//functions to be defined in derived classes
	virtual void reset_model(const double & C0, const double & start_inflation=0.0){};
	virtual double breath_step(const double & dvol, const double & dt){ return 0.0; }
	virtual void measure_values(MBWModelOutputs* output){};
	virtual void measure_MRI_dist(MBWModelOutputs* output){};
	virtual double get_mouth_conc(){ return 0.0; }

public:
	inline void set_inhaled_bc(const double & c){ this->inhaled_conc = c; }
	virtual void set_input_data(const MBWModelInputs * inputs);
	void simulate(MBWModelOutputs* output);
	//to be defined in derived classes
	virtual void build_model(const MBWModelOptions & opt, const std::vector<double> & paramsh,
		const std::vector<std::string> param_names_h){};
};

class CompartmentalModelBase: public MBWModelBase
{
protected:
	//base class for compartmental models of MBW
	std::shared_ptr<VentilationSolver> vent_solver;
	std::vector<std::shared_ptr<LungUnit>> units;
	std::vector<std::shared_ptr<DSVolume>> private_ds;
	std::shared_ptr<DSVolume> shared_ds;
	std::shared_ptr<MixingPoint> mixing_point;
	bool is_random, use_mouth_conc_generic;

	void reset_model(const double & C0, const double & start_inflation=0.0);
	double breath_step(const double & dvol, const double & dt);
	void measure_values(MBWModelOutputs* output);
	void measure_MRI_dist(MBWModelOutputs* output);
	double get_mouth_conc();
	void (*generate_vent_dist)(const std::map<std::string,double> &, 
							   const int &, std::vector<double> & );
	void (*initialise_lung_unit)(std::shared_ptr<LungUnit> & unit, const double & V0, 
						         const double & IGV0, const double & DV, 
						         const std::map<std::string,double> & params);

public:
	CompartmentalModelBase(){};
	void build_model(const MBWModelOptions & opt, const std::vector<double> & paramsh,
		                     const std::vector<std::string> param_names_h);
	double get_end_SDS_conc();
	double get_mouth_conc_generic();
	double get_igvol_rebreathe();
	void create_sync_vent_solver();
	void create_async_vent_solver(const double & delay_ts);
};

class CompartmentalModelGeneratorBase: 
	        public ModelGeneratorBase<CompartmentalModelBase, MBWModelInputs>
{
protected:
	std::vector<double> param_min, param_max;
	size_t N_FRC_mod_params;
	void initialise_FRC_modifier_params(MBWModelInputs* inputs);
public:
	CompartmentalModelGeneratorBase(MBWModelInputs* inputs):
		ModelGeneratorBase<CompartmentalModelBase, MBWModelInputs>(inputs)
	{
		this->initialise_FRC_modifier_params(inputs);
	}

	virtual void generate_from_prior(std::vector<double> & params_generated) const;
	virtual double prior_density(const std::vector<double> & params) const;
	virtual void generate_model(const std::vector<double> & params, 
		                        std::shared_ptr<CompartmentalModelBase> & m) const;
};

class BasicLognormalModelGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		BasicLognormalModel:
		3 parameters: FRC, Deadspace, sigma
		No shared DS, no asymm, no delays
	-------------------------------------------------------*/
protected:
	void initialise_params();
	void set_model_name();
public:
	BasicLognormalModelGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
};

class LognormalModelSDSGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		LognormalModelSDSGenerator:
		4 parameters: FRC, Deadspace, sigma, Shared deadspace fraction
		No asymm, no delays
	-------------------------------------------------------*/
protected:
	void set_model_name();
	void initialise_params();
public:
	LognormalModelSDSGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
};

class BimodalModelSDSGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		BimodalModelSDSGenerator:
		Ventilation drawn from a bimodal lognormal distribution
		7 parameters: FRC, Deadspace, sigma, Shared deadspace fraction, 
		ratio of mean vent in mode 1 to mode 2, ratio of sigma parameters,
	    fraction of lung in fast ventilating mode
		No asymm, no delays
	-------------------------------------------------------*/
protected:
	void set_model_name();
	void initialise_params();
public:
	BimodalModelSDSGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
	void generate_model(const std::vector<double> & params, 
		                std::shared_ptr<CompartmentalModelBase> & m) const;
};

class BasicLognormalAsymmModelGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		BasicLognormalAsymmModelGenerator:
		Ventilation drawn from a lognormal distribution
		5 parameters: FRC, Deadspace, sigma,
		asymmetry parameter, scale of diffusion parameter
		no delays
	-------------------------------------------------------*/
protected:
	void set_model_name();
	void initialise_params();
public:
	BasicLognormalAsymmModelGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
	void generate_model(const std::vector<double> & params, 
		                std::shared_ptr<CompartmentalModelBase> & m) const;
};

class LognormalAsymmSDSModelGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		LognormalAsymmSDSModelGenerator:
		Ventilation drawn from a lognormal distribution
		6 parameters: FRC, Deadspace, sigma, shared ds frac,
		asymmetry parameter, scale of diffusion parameter
		no delays
	-------------------------------------------------------*/
protected:
	void set_model_name();
	void initialise_params();
public:
	LognormalAsymmSDSModelGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
	void generate_model(const std::vector<double> & params, 
		                std::shared_ptr<CompartmentalModelBase> & m) const;
};

class BimodalAsymmSDSModelGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		BimodalAsymmSDSModelGenerator:
		Ventilation drawn from a lognormal distribution
		9 parameters: FRC, Deadspace, sigma, shared ds frac,
		ratio of mean vent in mode 1 to mode 2, ratio of sigma parameters,
	    fraction of lung in fast ventilating mode,
		asymmetry parameter, scale of diffusion parameter
		no delays
	-------------------------------------------------------*/
protected:
	void set_model_name();
	void initialise_params();
public:
	BimodalAsymmSDSModelGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
	void generate_model(const std::vector<double> & params, 
		                std::shared_ptr<CompartmentalModelBase> & m) const;
};

class BasicLognormalAsyncModelGenerator: public CompartmentalModelGeneratorBase
{
	/*------------------------------------------------------
		BasicLognormalAsyncModel
		Ventilation/conductance drawn from a lognormal distribution
		4 parameters: FRC, Deadspace, sigma, delay
		no asymm
	-------------------------------------------------------*/
protected:
	void set_model_name();
	void initialise_params();
public:
	BasicLognormalAsyncModelGenerator(MBWModelInputs* inputs):
		CompartmentalModelGeneratorBase(inputs)
	{
		this->set_model_name();
		this->initialise_params();
	}
	void generate_model(const std::vector<double> & params, 
		                std::shared_ptr<CompartmentalModelBase> & m) const;
};

void generate_vent_dist_lognormal_rand(const std::map<std::string,double> & params, 
								  const int & Ncomps, std::vector<double> & x_vals);

void generate_vent_dist_lognormal_determin(const std::map<std::string,double> & params, 
								  const int & Ncomps, std::vector<double> & x_vals);

void generate_vent_dist_bimodal_rand(const std::map<std::string,double> & params, 
								  const int & Ncomps, std::vector<double> & x_vals);

void build_basic_lung_unit(std::shared_ptr<LungUnit> & unit, const double & V0, 
						   const double & IGV0, const double & DV, 
						   const std::map<std::string,double> & params);

void build_asymm_lung_unit(std::shared_ptr<LungUnit> & unit, const double & V0, 
						   const double & IGV0, const double & DV, 
						   const std::map<std::string,double> & params);

void rescale_to_unit_mean(std::vector<double> & x);

void generate_async_delays(std::vector<std::shared_ptr<LungUnit>> & units, const double & delay_ts);

double get_mouth_conc_compartmental_generic(CompartmentalModelBase *m);
double get_end_SDS_conc_compartmental(CompartmentalModelBase *m);


class TrumpetModelBase: public MBWModelBase  
{
//network model based on trumpets to deal with diffusion
protected:
	Eigen::SparseMatrix<double> Incidence, ResLap, DifLap, AdvLap;

public:
	TrumpetModelBase(){};
	void build_model(const MBWModelOptions & opt, const std::vector<double> & paramsh,
		             const std::vector<std::string> param_names_h);

//	std::shared_ptr<VentilationSolver> vent_solver;
//	std::vector<std::shared_ptr<LungUnit>> units;
//	std::vector<std::shared_ptr<DSVolume>> private_ds;
//	std::shared_ptr<DSVolume> shared_ds;
//	std::shared_ptr<MixingPoint> mixing_point;
//	bool is_random, use_mouth_conc_generic;
//
//	void reset_model(const double & C0, const double & start_inflation=0.0);
//	double breath_step(const double & dvol, const double & dt);
//	void measure_values(MBWModelOutputs* output);
//	void measure_MRI_dist(MBWModelOutputs* output);
//	double get_mouth_conc();
//	void (*generate_vent_dist)(const std::map<std::string,double> &, 
//							   const int &, std::vector<double> & );
//	void (*initialise_lung_unit)(std::shared_ptr<LungUnit> & unit, const double & V0, 
//						         const double & IGV0, const double & DV, 
//						         const std::map<std::string,double> & params);
//
//public:
//	CompartmentalModelBase(){};
//	void build_model(const MBWModelOptions & opt, const std::vector<double> & paramsh,
//		                     const std::vector<std::string> param_names_h);
//	double get_end_SDS_conc();
//	double get_mouth_conc_generic();
//	double get_igvol_rebreathe();
//	void create_sync_vent_solver();
//	void create_async_vent_solver(const double & delay_ts);





};

#endif