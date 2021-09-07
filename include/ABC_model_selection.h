#ifndef ABC_MODEL_SELECTION_H
#define ABC_MODEL_SELECTION_H

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<memory>
#include<unordered_map>
#include<mpi.h>
#include<boost/multiprecision/random.hpp>
#include<boost/math/distributions/normal.hpp>
#include<boost/random/uniform_01.hpp>
#include<boost/random/normal_distribution.hpp>
#include<chrono>
#include<signal.h>
#include<algorithm>

double ks_statistic_normalised(const std::vector<double> & v0, const std::vector<double> & v1);

//main program object
template<class ModelGenerator, class Model, class DistanceFunction, 
         class ModelInputs, class ModelOutputs>
class ABCModelSelection
{
//class for running generic ABC model selection algorithm
protected:
	std::vector<std::shared_ptr<ModelGenerator>> ModelGen;
	std::shared_ptr<DistanceFunction> DistFunc;
	std::shared_ptr<ModelInputs> ModelInput;

	std::vector<int> model_counts, model_counts_previous; 
	std::vector<double> distances, weights, weights_previous;
	std::vector<std::vector<int>> ia_map, ia_map_previous;
	std::vector<std::vector<double>> params, params_previous, param_vars,
		                             prior_vars, param_tau;
	//[i][j] (i = model, j = sim_no or j = sim_no*param_count + param_no)
	
	int t_gen, n_cores, i_core, n_accept, n_per_proc, n_gens, n_cut, n_times_cutoff;
	double cutoff_distance, min_distance, ksstat_lim, ks_stat, ks_stat_old;

	//params is vector of parameters to be filled, return value is model no.
	int generate_model_id() const;
	void generate_from_previous(const int & model_id, std::vector<double> & params) const;
	double perturb_density(const int & im, const std::vector<double> & params_1, 
		                   const std::vector<double> & params_2) const;
	void print_gen_summary(const std::vector<std::vector<int>> ia_map_h,
						   const std::vector<std::vector<double>> params_h) const;
	void overwrite_finalgen_summary(const std::vector<std::vector<int>> ia_map_h,
						   const std::vector<std::vector<double>> params_h,
	                       const std::string & extra_headers,
	                       const std::vector<std::vector<std::string>> & extra_outputs) const;

	void update_weights();
	void update_param_summary();
	void update_distance_threshold();
	void kill_dead_models();
	void update_previous();

public:
	ABCModelSelection(){};
	ABCModelSelection(const std::vector<std::shared_ptr<ModelGenerator>> & MG,
		              std::shared_ptr<DistanceFunction> DF, 
					  std::shared_ptr<ModelInputs> MI)
	{
		this->setup(MG, DF, MI);
	}

	void initialise(const int & na, const int & ng, const int & ncut,
		            const double & min_d, const double & kslim, const int & ntlim);

	inline void setup(const std::vector<std::shared_ptr<ModelGenerator>> & MG,
		              std::shared_ptr<DistanceFunction> DF, 
					  std::shared_ptr<ModelInputs> MI)
	{
		this->ModelGen = MG;
		this->DistFunc = DF;
		this->ModelInput = MI;
	}

	void run();
};

//base classes

class ModelInputBase
{
	//base class for inputs to each model that remain unchanged
public:
	std::vector<double> measured;  //stores measured data
	std::string filepath;  //filepath to data directory
	ModelInputBase()
	{
		this->filepath = ".";
	}
	ModelInputBase(const std::string & fpath)
	{
		this->filepath = fpath;
	}
};

class ModelOutputBase
{
	//base class for outputs other than simulated data to be compared
public:
	std::vector<double> simulated;  //stores measured data
	ModelOutputBase(){};
	virtual bool extra_outputs() const;
	virtual void print_exta_outputs(std::string & line) const {};
	void get_headers(std::string & line) const {};
};  

class DistanceFunctionBase
{

public:
	virtual double distance(const std::vector<double> & measured, 
		                    const std::vector<double> & simulated);
};

template<class ModelOutputs> 
class ModelBase
{
//base class for model realisations
protected:
	std::vector<double> params;
public:
	inline double get_param(const int & i) const
	{ 
		if(i < this->params.size()) return (this->params[i]); 
		else
		{
			std::cerr << "Error: parameter does not exist, returning 0." << std::endl;
			return 0.;
		}
	}
	virtual void simulate(ModelOutputs* output){};
}; 

template<class Model, class ModelInputs>
class ModelGeneratorBase
{
//Base class for generation of model realisations. The model is built and returned by a function
//generate model, this takes 
protected:
	std::vector<std::string> param_names;
	std::string model_name;
	const ModelInputs *inputs;
public:
	ModelGeneratorBase(ModelInputs * in)
	{
		inputs = in;
	}
	virtual void generate_model(const std::vector<double> & params, 
		                        std::shared_ptr<Model> & m) const {};
	inline int get_param_count() const { return int(this->param_names.size()); }
	inline std::string get_param_name(const int & i) const
	{ 
		if(i < int(this->param_names.size())) return (this->param_names[i]); 
		else
		{
			std::cerr << "Error: parameter name does not exist, returning Null string." 
				      << std::endl;
			return std::string("");
		}
	}
	inline std::string get_model_name() const { return this->model_name; }
	//should be defined in derived class, but default form is given (U(0,1))
	virtual void generate_from_prior(std::vector<double> & params_generated ) const;
	virtual double prior_density(const std::vector<double> & params) const;
	
};

//function declarations

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     initialise(const int & na, const int & ng,
				const int & ncut, const double & min_d, const double & kslim, const int & ntlim)
{
	//get core info
	MPI_Comm_rank(MPI_COMM_WORLD, &i_core);
	MPI_Comm_size(MPI_COMM_WORLD, &n_cores);
	//setup random number generator
	//long unsigned int seed = static_cast<long unsigned int>
	//(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	//seed += i_core;   //to ensure seeds are different on different cores
	//rng = std::make_shared<boost::random::mt19937>();
	//rng->seed(seed);

	//initialise ABC settings
	this->n_accept = na;
	this->n_per_proc = int(na/n_cores);
	if(na%n_cores > 0) 
	{
		this->n_per_proc += 1;
		this->n_accept =  this->n_per_proc*n_cores;
		if(i_core == 0) std::cout << "N_ACCEPT updated to: " << this->n_accept << std::endl;
	}
	this->n_gens = ng;
	this->n_cut = ncut;
	this->min_distance = min_d;
	this->ksstat_lim = kslim;
	this->n_times_cutoff = ntlim;
	this->ks_stat = 2*kslim;
	this->ks_stat_old = 2*kslim;

	//initialise variables and storage
	this->distances.resize(n_accept);
	this->weights.resize(n_accept);
	this->weights_previous.resize(n_accept);
	this->ia_map.resize(this->ModelGen.size());
	this->params.resize(this->ModelGen.size());
	this->model_counts.resize(this->ModelGen.size());
	this->param_tau.resize(this->ModelGen.size());
	this->param_vars.resize(this->ModelGen.size());
	this->prior_vars.resize(this->ModelGen.size());
	MPI_Barrier(MPI_COMM_WORLD);
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     run()
{
	bool go = true;
	this->t_gen = 0;
	//std::vector<std::string> other_outputs;
	//std::vector<int> other_output_lengths;
	std::string other_headers;
	while(go)
	{
		std::vector<std::vector<double>> local_params(this->ModelGen.size());
		std::vector<std::vector<int>> local_ia_map(this->ModelGen.size());
		std::vector<double> local_distances(this->n_per_proc);
		std::vector<int> local_model_counts(this->ModelGen.size());
		std::vector<std::vector<std::string>> local_other_outputs(this->ModelGen.size());
		for(size_t im = 0; im < this->ModelGen.size(); im++)
		{
			local_params[im].reserve(this->n_accept);
			local_ia_map[im].reserve(this->n_accept);
			local_other_outputs[im].reserve(this->n_accept);
		}


		int local_tot_sims = 0;
		/*std::vector<int> local_other_output_lengths(n_per_proc);*/
		if(i_core == 0) std::cout << "---------------Initiating generation " << t_gen << "---------------" << std::endl;
		for(int ia = 0; ia < this->n_per_proc; ia++)
		{
			bool repeat = true;
			ModelOutputs output;
			double dist;
			int id;
			std::vector<double> pgen;
			while(repeat)
			{
				id = this->generate_model_id();
				if(this->t_gen == 0) 
				{
					 this->ModelGen[id]->generate_from_prior(pgen);
				}
				else
				{
					this->generate_from_previous(id, pgen);  //draw based on previous
				}
				if(this->ModelGen[id]->prior_density(pgen) > 0)   //redraw if not in prior
				{
					std::shared_ptr<Model> model;
					this->ModelGen[id]->generate_model(pgen, model);
					model->simulate(&output);
					local_tot_sims++;
					dist = this->DistFunc->distance(this->ModelInput->measured, output.simulated);
					if(t_gen == 0 || dist < this->cutoff_distance)
					{
						repeat = false;
					}
				}
			}
			if(output.extra_outputs()) 
			{
				std::string outstring;
				output.print_extra_outputs(outstring);
				local_other_outputs[id].push_back(outstring);
				if(i_core == 0 && ia == 0)
				{
					output.get_headers(other_headers);
				}
				//local_other_output_lengths[ia] = local_other_outputs[ia].size();
			}
			local_distances[ia] = dist;
			local_params[id].insert(local_params[id].end(),pgen.begin(),pgen.end());
			local_ia_map[id].push_back(i_core*n_per_proc + ia);
		}
		std::cout << "Generation " << t_gen << " complete on core " << i_core << "." << std::endl;
		for(int im = 0; im < int(local_params.size()); im++)
		{
			local_model_counts[im] = local_params[im].size() / this->ModelGen[im]->get_param_count();  
			//std::cout << "Core " << i_core << " model " << im << " count = " <<  local_model_counts[im] << std::endl;
		}

		//sum total model counts
		int total_sims;
		MPI_Allreduce(&local_tot_sims, &total_sims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(local_model_counts.data(), model_counts.data(), local_model_counts.size(), MPI_INT,
			          MPI_SUM, MPI_COMM_WORLD);
		MPI_Allgather(local_distances.data(), n_per_proc, MPI_DOUBLE, distances.data(),
			          n_per_proc, MPI_DOUBLE, MPI_COMM_WORLD);

		//-----print gen summary to screen-----//
		if(i_core == 0)
		{
			std::vector<std::string> outputs;
			int max_len = 0;
			for(int im = 0; im < int(this->ModelGen.size()); im++)
			{
				std::stringstream ss;
				ss << this->ModelGen[im]->get_model_name() << ": " 
				      << model_counts[im]; 
				outputs.push_back(ss.str());
				if(int(outputs[im].length()) > max_len) max_len = int(outputs[im].length());
			}
			std::stringstream ss;
			ss << "Generation " << t_gen << " models accepted";
			std::string hstr(ss.str());
			if(int(hstr.length()) > max_len) max_len = int(hstr.length());
			std::string hout;
			hout.resize(max_len + 4 - hstr.length(), '-');
			hout.insert((max_len + 4 - int(hstr.length()))/2,hstr);
			std::cout << hout << std::endl;
			for(int im = 0; im < int(this->ModelGen.size()); im++)
			{
				std::string ws;
				ws.resize(max_len - int(outputs[im].length()) + 1, ' ');
				std::cout << "| " << outputs[im] << ws << '|' << std::endl;
			}
			std::string ds;
			ds.resize(max_len + 4, '-');
			std::cout << ds << std::endl;
		}
		//----------------------------------//

		//resize arrays
		for(int im = 0; im < int(this->ModelGen.size()); im++)
		{
			this->ia_map[im].resize(this->model_counts[im]);
			this->params[im].resize(this->model_counts[im]*this->ModelGen[im]->get_param_count());
		}
		
		//gather ia_map and params to root node
		//std::cout << "Core" << i_core << ": gather params" << std::endl;
		//std::cout << this->n_cores << std::endl;
		for(int im = 0; im < int(this->ModelGen.size()); im++)
		{
			std::vector<int> counts(this->n_cores);
			std::vector<int> displs(this->n_cores); //resize
			std::vector<int> pcounts(this->n_cores);
			std::vector<int> pdispls(this->n_cores);
			int ppm = this->ModelGen[im]->get_param_count();  //params per model
			MPI_Gather(&local_model_counts[im], 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			//gather total counts for this model on core 0
			if(i_core == 0)   //only on core 0
			{
				displs[0] = 0;
				pdispls[0] = 0;
				pcounts[0] = ppm*counts[0];
				for(int ic = 1; ic < n_cores; ic++) 
				{
					displs[ic] = displs[ic-1] + counts[ic-1];   //cumulative model counts
					pcounts[ic] = ppm*counts[ic];   //parameter count
					pdispls[ic] = pdispls[ic-1] + pcounts[ic-1];  //cumulative param count
				}
			}
			//gather ia_map vectors into one
			//std::cout << i_core << ' ' << local_ia_map[im].size() << std::endl;
			//MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gatherv(local_ia_map[im].data(), local_ia_map[im].size(), MPI_INT,
						ia_map[im].data(), counts.data(), displs.data(), 
						MPI_INT, 0, MPI_COMM_WORLD);
			//std::cout << i_core << ' ' << local_params[im].size() << std::endl;
			//MPI_Barrier(MPI_COMM_WORLD);
			//gather param vectors into one
			MPI_Gatherv(local_params[im].data(), local_params[im].size(), MPI_DOUBLE, 
						params[im].data(), pcounts.data(), pdispls.data(), 
						MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		//if(i_core == 0) std::cout << "Params gathered." << std::endl;

		//broadcast params and ia maps
		for(int im = 0; im < int(this->ModelGen.size()); im++)
		{
			MPI_Bcast(this->ia_map[im].data(), this->ia_map[im].size(), MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(this->params[im].data(), this->params[im].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		this->update_weights();
		this->print_gen_summary(local_ia_map, local_params);
		this->overwrite_finalgen_summary(local_ia_map, local_params, other_headers, local_other_outputs);
		this->update_param_summary();
		if(t_gen == 0)
		{
			this->prior_vars = this->param_vars;
		}

		if(i_core == 0) 
		{
			std::cout << "Total sims generated in generation " << t_gen << ": " << total_sims << std::endl;
			std::cout << "Normalised K-S statistic: " << this->ks_stat << std::endl;
		}
		//conditions for loop termination
		bool term_cond1 = (this->t_gen > 0 && this->cutoff_distance <= this->min_distance);
		bool term_cond2 = (this->t_gen > 0 && this->t_gen == this->n_gens);
		bool term_cond3 = (this->t_gen > 1 && this->ks_stat <= this->ksstat_lim &&
			               this->ks_stat_old <= this->ksstat_lim);
		bool term_cond4 = (total_sims > this->n_times_cutoff*this->n_accept);
		if(term_cond1 || term_cond2 || term_cond3 || term_cond4)
		{
			go = false;
		}
		else
		{
			this->update_distance_threshold();
			if(i_core == 0) 
			{
				std::cout << "New distance threshold: " << this->cutoff_distance << std::endl;
			}
			this->kill_dead_models();
			this->update_previous();
		}

		this->t_gen++;
	}
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
int ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     generate_model_id() const
{
	//pick model at random
	boost::random::uniform_01<> udist;
	double rand = udist(*(rng.get()));
	for(int im = 0; im < int(this->ModelGen.size()); im++)
	{
		if(rand < double(im+1) / double(this->ModelGen.size()))
		{
			return im;
		}
	}
	return (int(this->ModelGen.size())-1);
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     generate_from_previous(const int & model_id, std::vector<double> & params) const
{
	//pick param set from model at random
	boost::random::uniform_01<> udist;
	double pset = udist(*(rng.get()));
	int nparams = this->ModelGen[model_id]->get_param_count();
	double wsum = 0;
	int ipset = -1;
	while(wsum < pset)
	{
		ipset++;
		wsum += this->weights_previous[ia_map_previous[model_id][ipset]];
	}
	params.assign(this->params_previous[model_id].begin()+nparams*ipset,
		          this->params_previous[model_id].begin()+nparams*(ipset+1));
	//perturb
	for(int ip = 0; ip < nparams; ip++)
	{
		boost::random::normal_distribution<> ndist(0,this->param_tau[model_id][ip]);
		params[ip] += ndist(*(rng.get()));
	}
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
double ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     perturb_density(const int & im, const std::vector<double> & params_1, const std::vector<double> & params_2) const
{
	double density = 1;
	for(int ip = 0; ip < int(params_1.size()); ip++)
	{
		boost::math::normal_distribution<> ndist(params_1[ip],this->param_tau[im][ip]);
		density *= boost::math::pdf(ndist, params_2[ip]);
	}
	return density;
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     print_gen_summary(const std::vector<std::vector<int>> ia_map_h,
					   const std::vector<std::vector<double>> params_h) const
{
	int number;
	//std::cout << "Core " << i_core << " ready to print" << std::endl;
	if(i_core > 0)
	{
		MPI_Recv(&number, 1, MPI_INT, i_core-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//std::cout << "Core " << i_core << "received ping from core " << i_core - 1 << std::endl;
	}

	for(int im = 0; im < int(this->ModelGen.size()); im++)
	{
		std::ofstream fo;
		std::stringstream ss;
		int nparam = this->ModelGen[im]->get_param_count();

		ss << this->ModelInput->filepath << "\\ABC_model_" << this->ModelGen[im]->get_model_name() 
			<< "_generation_" << this->t_gen << ".csv";
		//std::cout << "Core " << i_core << " Filepath :" << ss.str().c_str() << std::endl;
		if(i_core == 0) fo.open(ss.str().c_str(),std::fstream::out);
		else fo.open(ss.str().c_str(),std::fstream::app);

		//std::cout <<  "Core " << i_core << " File opened" << std::endl;
		//print headers
		if(i_core==0)
		{
			fo << "Distance,Weight";
			for(int ip = 0; ip < nparam; ip++)
			{
				fo << "," << this->ModelGen[im]->get_param_name(ip);
			}
			//fo << "," << extra_headers; 
			fo << std::endl;
		}
		for(int io = 0; io < int(ia_map_h[im].size()); io++)
		{
			fo << this->distances[ia_map_h[im][io]] << "," << this->weights[ia_map_h[im][io]];
			for(int ip = 0; ip < nparam; ip++)
			{
				fo << "," << params_h[im][io*nparam + ip];
			}
			//fo << "," << extra_outputs[im][io];
			fo << std::endl;
		}
		fo.close();
	}

	if(i_core < n_cores - 1) 
	{
		number = -1;
		MPI_Send(&number, 1, MPI_INT, i_core+1, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
 
template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     overwrite_finalgen_summary(const std::vector<std::vector<int>> ia_map_h,
					   const std::vector<std::vector<double>> params_h, const std::string & extra_headers,
	                   const std::vector<std::vector<std::string>> & extra_outputs) const
{
	int number;
	//std::cout << "Core " << i_core << " ready to print" << std::endl;
	if(i_core > 0)
	{
		MPI_Recv(&number, 1, MPI_INT, i_core-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//std::cout << "Core " << i_core << "received ping from core " << i_core - 1 << std::endl;
	}

	for(int im = 0; im < int(this->ModelGen.size()); im++)
	{
		std::ofstream fo;
		std::stringstream ss;
		int nparam = this->ModelGen[im]->get_param_count();

		ss << this->ModelInput->filepath << "\\ABC_model_" << this->ModelGen[im]->get_model_name() 
			<< "_summary.csv";
		//std::cout << "Core " << i_core << " Filepath :" << ss.str().c_str() << std::endl;
		if(i_core == 0) fo.open(ss.str().c_str(),std::fstream::out);
		else fo.open(ss.str().c_str(),std::fstream::app);

		//std::cout <<  "Core " << i_core << " File opened" << std::endl;
		//print headers
		if(i_core==0)
		{
			fo << "Distance,Weight";
			for(int ip = 0; ip < nparam; ip++)
			{
				fo << "," << this->ModelGen[im]->get_param_name(ip);
			}
			fo << "," << extra_headers << std::endl;
		}
		for(int io = 0; io < int(ia_map_h[im].size()); io++)
		{
			fo << this->distances[ia_map_h[im][io]] << "," << this->weights[ia_map_h[im][io]];
			for(int ip = 0; ip < nparam; ip++)
			{
				fo << "," << params_h[im][io*nparam + ip];
			}
			fo << "," << extra_outputs[im][io] << std::endl;
		}
		fo.close();
	}

	if(i_core < n_cores - 1) 
	{
		number = -1;
		MPI_Send(&number, 1, MPI_INT, i_core+1, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
	 update_param_summary()
{
	for(int im = 0; im < int(this->ModelGen.size()); im++)  //loop over models (im)
	{
		int model_pc = this->ModelGen[im]->get_param_count();
		this->param_tau[im].resize(model_pc); 
		this->param_vars[im].resize(model_pc);
		if(i_core==0)
		{
			if(this->model_counts[im] > 1)      //if model is not extinct
			{
				for(int ip = 0; ip < model_pc; ip++)   //loop over params (ip)
				{
					double m1 = 0, m2 = 0;
					for(int ia = 0; ia < this->model_counts[im]; ia++)
					{
						int iap = ia*model_pc + ip;
						m1 += this->params[im][iap];
						m2 += this->params[im][iap]*this->params[im][iap];
					}
					m1 /= double(this->model_counts[im]);
					m2 /= double(this->model_counts[im]);
					this->param_vars[im][ip] = m2 - m1*m1;
					this->param_tau[im][ip] = sqrt(2*this->param_vars[im][ip]);    
				}
			}
			else     //if model does not have enough points to compute variance
			{
				for(int ip = 0; ip < model_pc; ip++)   //loop over params (ip)
				{
					this->param_vars[im][ip] = 0.0;
					this->param_tau[im][ip] = 1.0;//sqrt(2*this->prior_vars[im][ip]);
				}
			}
		}
		MPI_Bcast(param_vars[im].data(), model_pc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(param_tau[im].data(), model_pc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	if(t_gen > 0 && i_core==0)
	{
		this->ks_stat = 0;
		//for each param, integrate over all others
		for(int im = 0; im < int(this->ModelGen.size()); im++)  //loop over models (im)
		{
			int model_pc = this->ModelGen[im]->get_param_count();
			if(this->model_counts[im] > 0)      //if model is not extinct
			{
				for(int ip = 0; ip < model_pc; ip++)   //loop over params (ip)
				{
					std::vector<double> p_new(this->model_counts[im]), p_old(this->model_counts_previous[im]);
					for (int ia = 0; ia < this->model_counts[im]; ia++)
					{
						int iap = ia*model_pc + ip;
						p_new[ia] = this->params[im][iap];
					}
					for (int ia = 0; ia < this->model_counts_previous[im]; ia++)
					{
						int iap = ia*model_pc + ip;
						p_old[ia] = this->params_previous[im][iap];
					}
					double ks_stath = ks_statistic_normalised(p_new,p_old);
					if (ks_stath > this->ks_stat) this->ks_stat = ks_stath;
				}
			}
			double dmc = abs(double(this->model_counts[im] - this->model_counts_previous[im]));
			double ks_stat_ms = dmc/sqrt(2*n_accept);
			if (ks_stat_ms > this->ks_stat) this->ks_stat = ks_stat_ms;
		}
	}
	MPI_Bcast(&ks_stat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     update_weights()
{
	for(int im = 0; im < int(this->ModelGen.size()); im++)
	{
		if(i_core==0)
		{
			int nparams = this->ModelGen[im]->get_param_count();
			if(t_gen > 0)
			{
				double tot_weight = 0;
				for(int ia = 0; ia < this->model_counts[im]; ia++)
				{
					std::vector<double> p2;
					p2.assign(params[im].begin() + ia*nparams, 
							  params[im].begin() + (ia+1)*nparams);
					int iw = ia_map[im][ia];
					this->weights[iw] = 0;
					for(int iao = 0; iao < this->model_counts_previous[im]; iao++)
					{
						std::vector<double> p1;
						p1.assign(params_previous[im].begin() + iao*nparams, 
							      params_previous[im].begin() + (iao+1)*nparams);
						
						this->weights[iw] += (weights_previous[ia_map_previous[im][iao]]
						                      *perturb_density(im,p1,p2));
					}
					this->weights[iw] = this->ModelGen[im]->prior_density(p2) / this->weights[iw];
					tot_weight += this->weights[iw];
				}
				//normalise
				for(int ia = 0; ia < this->model_counts[im]; ia++)
				{
					int iw = ia_map[im][ia];
					this->weights[iw] /= tot_weight;
				}
			}
			else
			{
				for(int ia = 0; ia < this->model_counts[im]; ia++)
				{
					int iw = ia_map[im][ia];
					this->weights[iw] = 1.0 / double(this->model_counts[im]);
				}
			}
		}
	}
	MPI_Bcast(weights.data(), weights.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     update_distance_threshold()
{
	if(i_core == 0)
	{
		std::vector<double> dist_copy = distances;
		std::sort(dist_copy.begin(),dist_copy.end(),std::less<double>());
		this->cutoff_distance = dist_copy[n_cut];
	}
	MPI_Bcast(&cutoff_distance,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     kill_dead_models()
{
	int nmodels = this->ModelGen.size();
	std::vector<int> to_kill;
	for(int im = 0; im < nmodels; im++)
	{
		if(this->model_counts[im] == 0) to_kill.push_back(im);
	}
	for(int ikill = int(to_kill.size())-1; ikill >= 0; ikill--)
	{
		this->ModelGen.erase(this->ModelGen.begin()+to_kill[ikill]);
		this->model_counts.erase(this->model_counts.begin()+to_kill[ikill]);
		this->params.erase(this->params.begin()+to_kill[ikill]);
		this->ia_map.erase(this->ia_map.begin()+to_kill[ikill]);
		this->param_tau.erase(this->param_tau.begin()+to_kill[ikill]);
		this->param_vars.erase(this->param_vars.begin()+to_kill[ikill]);
		this->prior_vars.erase(this->prior_vars.begin()+to_kill[ikill]);
	}
}

template<class ModelGenerator, class Model, class DistanceFunction,
         class ModelInputs, class ModelOutputs>
void ABCModelSelection<ModelGenerator,Model,DistanceFunction,ModelInputs,ModelOutputs>::
     update_previous()
{
	params_previous = params;
	weights_previous = weights;
	ia_map_previous = ia_map;
	model_counts_previous = model_counts;
	ks_stat_old = ks_stat;
	//std::cout << "weights count " << weights.size() << std::endl;
	/*for(int im = 0; im < int(this->ModelGen.size()); im++)
	{
		std::cout << "Core: " << i_core << " Model " << im 
			      << " param count " << params[im].size()
				  << " ia_map_count " << ia_map[im].size()
				  << " model_count " << model_counts[im] << std::endl;

	}*/
}

template<class Model, class ModelInputs>
void ModelGeneratorBase<Model,ModelInputs>::
	 generate_from_prior(std::vector<double> & params_generated) const
{
	extern std::shared_ptr<boost::random::mt19937> rng;
	//default prior is U(0,1) for all parameters
	boost::random::uniform_01<> udist;
	int np = int(this->param_names.size());
	params_generated.resize(np);
	for(int ip = 0; ip < np; ip++)
	{
		params_generated[ip] = udist(*(rng.get()));
	}
}

template<class Model, class ModelInputs>
double ModelGeneratorBase<Model,ModelInputs>::
	 prior_density(const std::vector<double> & params) const
{
	//returns (non-normalsied) probability density of priors
	for(int ip = 0; ip < int(params.size()); ip++)
	{
		if(params[ip] < 0 || params[ip] > 1) return 0.;
	}
	return 1.;
}


//template<class Model, class ModelInputs>
//void ModelGeneratorBase<Model,ModelInputs>::
//     update_param_summary(const std::vector<std::vector<double>> & param_list)
//{
//	//default param summary is vairance of parameters
//	int nsim = int(param_list.size());
//	int np = int(this->param_names.size());
//	this->param_vars.resize(np);
//	for(int ip = 0; ip < np; ip++)
//	{
//		double m1 = 0, m2 = 0;
//		for(int isim = 0; isim < nsim; isim++)
//		{
//			m1 += param_list[isim][ip];
//			m2 += param_list[isim][ip]*param_list[isim][ip];
//		}
//		m1 /= nsim;
//		m2 /= nsim;
//
//		this->param_vars[ip] = m2 - m1*m1;
//	}
//}
//	
//template<class Model, class ModelInputs>
//void ModelGeneratorBase<Model,ModelInputs>::
//	 perturb_param_set(std::vector<double> & params)
//{	
//	extern std::shared_ptr<boost::random::mt19937> rng;  //pseudorandom number generator
//	//perturb parameters with uniform distribution
//	if(params.size() != this->param_vars.size())
//	{
//		std::cerr << "Error, incorrect parameter count, perturbation failed." << std::endl;
//	}
//	else
//	{
//		boost::random::uniform_01<> udist;
//		for(int ip = 0; ip < int(this->param_vars.size()); ip++)
//		{
//			params[ip] += (udist(*(rng.get())) - 0.5)*sqrt(24*this->param_vars[ip]);
//		}
//	}
//}
//	
//template<class Model, class ModelInputs>
//double ModelGeneratorBase<Model,ModelInputs>::
//       perturbation_density(const std::vector<double> & new_params, 
//		                    const std::vector<double> & old_params)
//{
//	//returns (non-normalsied) probability density of perturbation
//	for(int ip = 0; ip < int(this->param_vars.size()); ip++)
//	{
//		if(abs(new_params[ip] - old_params[ip]) > 0.5*sqrt(24*this->param_vars[ip]))
//		{
//			return 0.;
//		}
//	}
//	return 1.;
//}


#endif