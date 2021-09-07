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
#include<Windows.h>
#include<signal.h>
#include<process.h>

std::shared_ptr<boost::random::mt19937> rng;  //pseudorandom number generator

inline static HANDLE start_mpicode(const std::vector<std::string> & args)
{
	std::stringstream ss;
	int argc = int(args.size());
	ss << args[0];
	for(int i = 1; i < argc; i++)
	{
		ss << " " << args[i];
	}
	std::string cmd = ss.str().c_str();
    PROCESS_INFORMATION pi = {0};
    STARTUPINFO si = {0};
    si.cb = sizeof(STARTUPINFO);
	size_t size = cmd.size();
	char* orig = new char [size+1];
	copy(cmd.begin(), cmd.end(), orig);
	orig[size] = '\0';
	std::cout << cmd << std::endl;
	wchar_t* input = new wchar_t [size+1];
	size_t convertedChars = 0;
    mbstowcs_s(&convertedChars, input, size+1, orig, _TRUNCATE);
	delete [] orig;
    if (CreateProcess(NULL, input, NULL, NULL, TRUE, 0, NULL, NULL, &si, &pi))
    {
		delete [] input;
        CloseHandle(pi.hThread);
        return pi.hProcess;
    }
	else
	{
		delete [] input;
		return NULL;
	}
}

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
	
	//check MPI has been used, if not, use it
	if(n_cores != N_PROCS)
	{
		//restart
		std::vector<std::string> args(3 + argc);
		args[0]="mpiexec";// -debug";
		args[1]="-n";
		std::stringstream ss;
		ss << N_PROCS;
		args[2]=ss.str().c_str();

		for (int i=0;i<argc;i++) 
		{
			args[3+i]=argv[i];
		}

		//function to set mpi executable running in separate process (windows only)
		HANDLE ff = start_mpicode(args);
		if(ff != NULL)    //do not exit this code until child process is complete
		{
			// wait with ten-second checks
			while (WAIT_TIMEOUT == WaitForSingleObject(ff, 10000))
			{
				// your wait code goes here.
			}

			// close the handle no matter what else.
			CloseHandle(ff);
			return 0;
		}
		else
		{
			std::cerr<<"warning: mpiexec failed"<<std::endl;
			return 1;
		}
	}

	//setup rng (different on each core)
	long unsigned int seed = static_cast<long unsigned int>
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

	//test code
	/*vector<ModelGenBaseClassName*> model_gens;
	model_gens.resize(2);
	std::shared_ptr<ModelGen1Class> modelgen1 = std::make_shared<ModelGen1Class>(inputs.get());
	std::shared_ptr<ModelGen2Class> modelgen2 = std::make_shared<ModelGen2Class>(inputs.get());
	model_gens[0] = modelgen1.get();
	model_gens[1] = modelgen2.get();

	std::shared_ptr<CompartmentalModelBase> model1 = std::make_shared<CompartmentalModelBase>();
	std::vector<double> params;
	params.push_back(2.0);
	params.push_back(0.15);
	params.push_back(0.8);
	model_gens[0]->generate_model(params, model1.get());

	std::shared_ptr<CompartmentalModelBase> model2 = std::make_shared<CompartmentalModelBase>();
	params.push_back(0.5);
	model_gens[1]->generate_model(params, model2.get());*/



	

	////get input path
	//string input_path = argv[argc-1]; 
	//if(i_core == 0) cout << "Input file path is: " << input_path << endl;

	////initialise ABC object
	//vector<shared_ptr<ModelGenBaseClassName>> ModelGens;
	//assign_model_gens(ModelGens);
	//shared_ptr<DistanceName> DistFunc = make_shared<DistanceName>();
	//shared_ptr<ModelInputName> ModelIn = make_shared<ModelInputName>(input_path);

	//ABCModelSelection<ModelGenBaseClassName, ModelBaseClassName, DistanceName,
	//	    ModelInputName, ModelOutputName> ABCalgorithm(ModelGens, DistFunc, ModelIn);
	//
	//ABCalgorithm.initialise(N_ACCEPT, N_GENS, N_CUT, MIN_DIST, N_TIMES_CUTOFF);
	//
	////run ABC algorithm
	//MPI_Barrier(MPI_COMM_WORLD);

	//ABCalgorithm.run();

	MPI_Finalize();
	return 0;
}