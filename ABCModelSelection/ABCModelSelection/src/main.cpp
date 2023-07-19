//Need to consider best way to organise the models:
//I think each model should have options specified by an input file, which is
//then sorted out by the program.
//Options: Ventilation distribution type, shared-DS parameter, phase-III mechanisms, 
//gas exchange, pulmonary blood compartment
//Priors always assigned in the same way for each model
//Where applicable these options should be contained within the model class object,
//which remains unchanged and is used to generate model relaisations.
//model realisation is a separate object consisting of dead-space vol objects, lung unit objects,
//(essentially the object we already have)

//Create this first based on AL3C stuff
//ABCModelComparison object -- manages parameter generation, priors etc. Generic
//contains model[1] model[2] etc. each of which contain paramA, paramB, priorA, priorB, pertA, 
//pertB, lastA, lastB etc. and model generator. Manages all parallel stuff 

//Data reader and MBW specfic set up -- objects to read &  store data
//functions to use data and input files and generate priors and perturbation kernels for all
//parameters of all models

//ModelGenerator -- global -- parameter set sent to it and model object produced. Stores data
//input required for model too, and tranfers to generated model + noise

//Model objects -- local -- run washout, measure distance etc.

#include"ABC_model_selection.h"
#include"ABCsettings.h"

std::shared_ptr<boost::random::mt19937> rng;  //pseudorandom number generator

void assign_model_gens(MBWModelInputs *inputs,
	                   std::vector<std::shared_ptr<ModelGenBaseClassName>> & ModelGens)
{
ModelGens.push_back(std::make_shared<ModelGen1Class>(inputs));
#if N_MODELS > 1
ModelGens.push_back(std::make_shared<ModelGen2Class>(inputs));
#endif
#if N_MODELS > 2
ModelGens.push_back(std::make_shared<ModelGen3Class>(inputs));
#endif
#if N_MODELS > 3
ModelGens.push_back(std::make_shared<ModelGen4Class>(inputs));
#endif
#if N_MODELS > 4
ModelGens.push_back(std::make_shared<ModelGen5Class>(inputs));
#endif
#if N_MODELS > 5
std::cerr << "Too many models, igoring all models above model 5" << std::endl;
#endif
}

int main(int argc, char *argv[])
{
	//init MPI
	//signal(SIGINT, signal_callback_handler);
	MPI_Init(NULL,NULL);

	int i_core, n_cores;
	MPI_Comm_rank(MPI_COMM_WORLD, &i_core);
	MPI_Comm_size(MPI_COMM_WORLD, &n_cores);
	
	//setup rng (different on each core)
	unsigned int seed = static_cast<unsigned int>
			   (std::chrono::high_resolution_clock::now().time_since_epoch().count());
	seed += i_core;   //to ensure seeds are different on different cores
	rng = std::make_shared<boost::random::mt19937>();
	rng->seed(seed);

	//read inputs
	std::string input("");
	if( argc > 1 ) input = std::string(argv[1]);
	else
	{
		std::cout << "Please enter path for input files:";
		std::cin >> input;
	}
		
	std::shared_ptr<ModelInputName> inputs = std::make_shared<ModelInputName>(input);
	std::shared_ptr<DistanceName> dist = std::make_shared<DistanceName>();

	std::vector<std::shared_ptr<ModelGenBaseClassName>> model_gens;
	assign_model_gens(inputs.get(),model_gens);

	ABCModelSelection<ModelGenBaseClassName,ModelBaseClassName,DistanceName,
						ModelInputName, ModelOutputName> algorithm(model_gens,dist,inputs);
	//intialise algortihm for each processor
	algorithm.initialise(N_ACCEPT,N_GENS,N_CUT,MIN_DIST,KS_STAT_CUTOFF,N_TIMES_CUTOFF);
	//start
	algorithm.run();

	//print summary?

	MPI_Finalize();
	return 0;
}
