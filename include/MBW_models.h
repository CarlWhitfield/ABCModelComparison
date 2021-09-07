#ifndef MBW_MODELS_H
#define MBW_MODELS_H

#include"ABC_model_selection.h"
#include"compartmental.h"
#include<limits>
#include<map>

//quantities that are neither model nor patient specific
const double SF6_NOISE = 0.0002;   //largest estimate from SNR in Horsley et al. Thorax 2007
const double VOL_NOISE = 0.01;    //SNR for **volume** from same paper, note this probably depends on flow rate
const double MRI_NOISE_FRAC = 0.02;    //1/SNR expected for MRI for mean -- ballpark figure
const double VENT_RATIO_LIMIT = 100;   //maximum ventilation ratio (i.e. truncation of ventilation distribution)

//default model options
const int NCOMPS = 50;  //
const int NMRI_SAMPLES = 1000;
const double NMV_FRAC_STEP = 0.2;

//prior limits
const double FRC_MIN_FRAC = 0.5;   //as a fraction of MBW measured FRC
const double FRC_MAX_FRAC = 2.0;
const double VD_MIN_FRAC = 0.5;   //as a fraction of MBW measured FDS
const double VD_MAX_FRAC = 3.0;
const double SIGMA_MIN = 0;   //as a fraction of MBW measured FDS
const double SIGMA_MAX = 4.0;
const double MURATIO_MAX = 20.0;  //increase if posterior suggests higher
const double SIGRATIO_MAX = 20.0;
const double DIFFSCALE_MIN = 1E-02;   //min time resolution
const double DIFFSCALE_MAX = 20;    //largest plausible diffusion timescale (scale of 2 breaths)

#define FRC_PARAM_NAME "FRC"
#define VD_PARAM_NAME "VD"
#define VDSFRAC_PARAM_NAME "VD_shared_fraction"
#define SIGMA_PARAM_NAME "Sigma"
#define MURATIO_LS_PARAM_NAME "Mu_ratio"
#define SIGRATIO_LS_PARAM_NAME "Sigma_ratio"
#define VFASTFRAC_PARAM_NAME "V_fast_frac"
#define ASYMM_PARAM_NAME "Acin_Asymm"
#define DIFFSCALE_PARAM_NAME "Acin_diff_tscale"

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

class RMSEDistance: public DistanceFunctionBase{};

struct MBWModelOptions
{
public:
	//model specific options
	char vent_dist_type, lung_unit_type;
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
		conc_measurements, dist_weights;
	std::vector<std::vector<int>> measurement_steps;
	std::vector<double> Cinit, Cinit_std;
	double dead_space, FRC0, machine_ds, MRI_vbag; 
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

class CompartmentalModelBase: public ModelBase<MBWModelOutputs>
{
protected:
	//base class for compartmental models of MBW
	std::vector<std::shared_ptr<LungUnit>> units;
	std::vector<std::shared_ptr<DSVolume>> private_ds;
	std::shared_ptr<DSVolume> shared_ds;
	std::shared_ptr<MixingPoint> mixing_point;
	bool is_random;
	double mouth_point, inhaled_conc;
	double washout_start_inflation, washout_start_inflation_duration;
	bool simulate_washin;
	int washout_start_timepoint;

	//pointers to data from MBW input and options (stored in generator)
	//non-pointers to data that is perturbed before use
	std::vector<std::vector<double>> sim_conc, sim_vol_steps;
	std::vector<double> Cinit;
	const std::vector<std::vector<double>> *sim_step_durations;
	const std::vector<std::vector<int>> *measurement_steps;
	double MRI_bag_vol;
	int NMRIPoints;

	virtual void measure_values(MBWModelOutputs* output);
//	//functions to be defined in derived classes
	virtual void measure_MRI_dist(MBWModelOutputs* output);
	void run_washout_model(MBWModelOutputs* output);
	void run_MRI_model(MBWModelOutputs* output);
	//virtual void reassign_vent_ratios();

	void compute_volume_changes(const double & dvol, const double & dtime, 
								std::vector<double> & vols);
	void reset_model(const double & C0, const double & start_inflation=0.0, 
		             const double & start_inflation_duration=0.0);
	void breath_step(const double & dvol, const double & dt);
	void (*generate_vent_dist)(const std::map<std::string,double> &, 
							   const int &, std::vector<double> & );
	void (*initialise_lung_unit)(std::shared_ptr<LungUnit> & unit, const double & V0, 
						         const double & IGV0, const double & DV, 
						         const std::map<std::string,double> & params);
	double (CompartmentalModelBase::*get_mouth_conc)();

public:
	CompartmentalModelBase(){};
	inline void set_inhaled_bc(const double & c){ this->inhaled_conc = c; }
	virtual void set_input_data(const MBWModelInputs * inputs);
	virtual void build_model(const MBWModelOptions & opt, const std::vector<double> & paramsh,
		                     const std::vector<std::string> param_names_h);
	double get_end_SDS_conc();
	double get_mouth_conc_generic();
	virtual void simulate(MBWModelOutputs* output);
};

class CompartmentalModelGeneratorBase: 
	        public ModelGeneratorBase<CompartmentalModelBase, MBWModelInputs>
{
protected:
	std::vector<double> param_min, param_max;
public:
	CompartmentalModelGeneratorBase(MBWModelInputs* inputs):
		ModelGeneratorBase<CompartmentalModelBase, MBWModelInputs>(inputs){};

	virtual void generate_from_prior(std::vector<double> & params_generated) const;
	virtual double prior_density(const std::vector<double> & params) const;
	virtual void generate_model(const std::vector<double> & params, 
		                        std::shared_ptr<CompartmentalModelBase> & m) const;
};

class BasicLognormalModelGenerator: public CompartmentalModelGeneratorBase
{
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

#endif