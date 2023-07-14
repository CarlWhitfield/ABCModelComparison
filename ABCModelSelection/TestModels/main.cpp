#include"ABC_model_selection.h"
#include"ABCsettings.h"
#include<Windows.h>
#include<signal.h>
#include<process.h>
#include <chrono>

std::shared_ptr<boost::random::mt19937> rng;  //pseudorandom number generator

int main(int argc, char *argv[])
{
	long unsigned int seed = static_cast<long unsigned int>
			   (std::chrono::high_resolution_clock::now().time_since_epoch().count());
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

	std::vector<ModelGenBaseClassName*> model_gens;
	model_gens.resize(2);
	std::shared_ptr<ModelGen1Class> modelgen1 = std::make_shared<ModelGen1Class>(inputs.get());
	std::shared_ptr<ModelGen2Class> modelgen2 = std::make_shared<ModelGen2Class>(inputs.get());
	model_gens[0] = modelgen1.get();
	model_gens[1] = modelgen2.get();

	std::shared_ptr<CompartmentalModelBase> model1 = std::make_shared<CompartmentalModelBase>();
	std::vector<double> params;
	params.push_back(2.0);
	params.push_back(0.15);
	params.push_back(1.5);
	model_gens[0]->generate_model(params, model1);

	std::shared_ptr<CompartmentalModelBase> model2 = std::make_shared<CompartmentalModelBase>();
	params.push_back(5.0);
	model_gens[1]->generate_model(params, model2);

	std::shared_ptr<MBWModelOutputs> outputs1, outputs2;
	outputs1 = std::make_shared<MBWModelOutputs>();
	outputs2 = std::make_shared<MBWModelOutputs>();

	auto start = std::chrono::high_resolution_clock::now();
	model1->simulate(outputs1.get());
	auto endt = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endt - start);
	std::cout << "Model 1 took: " << duration.count()/1E6 << "s." << std::endl;

	start = std::chrono::high_resolution_clock::now();
	model2->simulate(outputs2.get());
	endt = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(endt - start);
	std::cout << "Model 2 took: " << duration.count()/1E6 << "s." << std::endl;

	//write to csv
	size_t imtot = 0;
	for(size_t i = 0; i < inputs->sim_vol_steps.size(); i++)
	{
		std::ofstream outfile;
		std::stringstream fname;
		fname << "test_model_comp_T" << i+1 << ".csv";
		outfile.open(fname.str().c_str());
		outfile << "Time (s), Model 1, Model 2" << std::endl;
		double t = 0;
		size_t im = 0;
		for(size_t n = 0; n < inputs->sim_vol_steps[i].size(); n++)
		{
			t += inputs->sim_step_durations[i][n];
			if(im < inputs->measurement_steps[i].size() && n == inputs->measurement_steps[i][im])
			{
				outfile << t << ", " << outputs1->simulated[imtot + im] << ", " << outputs2->simulated[imtot + im] << std::endl;
				im++;
			}
		}
		imtot += im;
		outfile.close();
	}
	
	return 0;
}