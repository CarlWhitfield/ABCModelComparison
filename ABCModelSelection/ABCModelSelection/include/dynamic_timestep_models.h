#ifndef MODELS2_H
#define MODELS2_H

#include<vector>
#include<string>
#include<sstream>
#include<map>
#include<chrono>
#include<Windows.h>
#include"compartmental.h"
#include"read_mbw_data.h"
#include<boost/multiprecision/random.hpp>
#include<boost/random/lognormal_distribution.hpp>
#include<boost/random/normal_distribution.hpp>
#include<boost/random/uniform_01.hpp>
#include<boost/shared_array.hpp>
#include<boost/shared_ptr.hpp>
#include<boost/make_shared.hpp>








//classes for running improved models
//base class for models
class ResCompartmentalModel: public CompartmentalModel
{
protected:
	void compute_volume_changes(const double & dvol, const double & dtime, 
		                        std::vector<double> & vols);
	void setup_warmup_breaths(boost::shared_ptr<MBWData> data);
	std::vector<double> tau, tau_orig, warm_up_durations, warm_up_steps;
public:	
	ResCompartmentalModel():CompartmentalModel(){};
	ResCompartmentalModel(boost::shared_ptr<MBWData> data):CompartmentalModel(data)
	{
		this->setup_warmup_breaths(data);
	 }
	ResCompartmentalModel(const std::map<std::string,double> & p, const double & VT, 
		               const double & Btime, boost::shared_ptr<MBWData> data)
	{
		this->build_test_case(p,VT,Btime,data);
		this->setup_warmup_breaths(data);
	}
	~ResCompartmentalModel(){};

	void build_model(const std::map<std::string,double> & params,
		             const double & machine_ds);
	void reset_model(const double & C0);
	void reassign_vent_ratios();
	void run_warm_up();
};

class BimodalCompartmentalModel: public CompartmentalModel
{
protected:
	void generate_vent_dist(const std::map<std::string,double> & params, 
							const int & Ncomps, const double & Vbag);
	double trunc_normal_mu(const double & mean, const double & sig);
public:
	BimodalCompartmentalModel():CompartmentalModel(){};
	BimodalCompartmentalModel(boost::shared_ptr<MBWData> data):CompartmentalModel(data){};
	BimodalCompartmentalModel(const std::map<std::string,double> & p, const double & VT, 
		                      const double & Btime, boost::shared_ptr<MBWData> mbw_to_fill)
	{
		this->is_random = true;
		this->build_test_case(p,VT,Btime,mbw_to_fill);
	}
	~BimodalCompartmentalModel(){};
};

class SDSCompartmentalModel: public CompartmentalModel
{
public:
	SDSCompartmentalModel():CompartmentalModel(){};
	SDSCompartmentalModel(boost::shared_ptr<MBWData> data):CompartmentalModel(data){};
	SDSCompartmentalModel(const std::map<std::string,double> & p, const double & VT, 
		               const double & Btime, boost::shared_ptr<MBWData> mbw_to_fill)
	{
		this->is_random = true;
		this->build_test_case(p,VT,Btime,mbw_to_fill);
	};
	~SDSCompartmentalModel(){};

	double get_mouth_conc();
	void build_model(const std::map<std::string,double> & params,
		                     const double & machine_ds);
};

class AcinCompartmentalModel: public CompartmentalModel
{
public:
	AcinCompartmentalModel():CompartmentalModel(){};
	AcinCompartmentalModel(boost::shared_ptr<MBWData> data):CompartmentalModel(data){};
	AcinCompartmentalModel(const std::map<std::string,double> & p, const double & VT, 
		               const double & Btime, boost::shared_ptr<MBWData> mbw_to_fill)
	{
		this->is_random = true;
		this->build_test_case(p,VT,Btime,mbw_to_fill);
	}
	~AcinCompartmentalModel(){};

	//void build_model(const std::map<std::string,double> & params,
	//	                     const double & machine_ds);
};


#endif