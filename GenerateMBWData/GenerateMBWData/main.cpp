#include<Windows.h>
#include<signal.h>
#include<read_write_codes.h>
#include"MBW_models.h"

std::shared_ptr<boost::random::mt19937> rng;  //pseudorandom number generator

int main(int argc, char *argv[])
{
	//argv[1] should be a .params file

	//setup rng
	long unsigned int seed = static_cast<long unsigned int>
			   (std::chrono::high_resolution_clock::now().time_since_epoch().count());
	//seed += i_core;   //to ensure seeds are different on different cores
	rng = std::make_shared<boost::random::mt19937>();
	rng->seed(seed);

	//read params file
	std::string param_file("");
	if( argc > 1 ) param_file = std::string(argv[1]);
	else
	{
		std::cout << "Please enter path to .params file:";
		std::cin >> param_file;
	}
	std::cout << "Param file name: " << param_file << std::endl;

	std::shared_ptr<MBWModelInputs> inputs = std::make_shared<MBWModelInputs>();
	MBWModelOptions model_options;
	std::vector<double> model_params;
	std::vector<std::string> model_param_names;
	inputs->generate_inputs(param_file, model_options, model_params, model_param_names);

	std::cout << "Inputs generated" << std::endl;
	//generate model
	CompartmentalModelBase model_h;
	model_h.set_input_data(inputs.get());
	model_h.build_model(model_options, model_params, model_param_names);

	std::cout << "Model built" << std::endl;
	std::shared_ptr<MBWModelOutputs> model_outputs = std::make_shared<MBWModelOutputs>();
	//run model
	model_h.simulate(model_outputs.get());

	std::cout << "Model built" << std::endl;
	//print outputs
	std::string fhead = param_file;
	fhead.erase(param_file.length() - 7, 7);
	int i_offset = 0;
	for(int n = 0; n < inputs->Ntests; n++)
	{
		std::stringstream ss;
		ss << fhead << "-" << n << ".txt";
		std::ofstream fout(ss.str().c_str());
		fout << "Raw data	******************************************************\n";
		fout << "Time[ms]\t" << "Flow\t" << "SF6\t" << "O2\t" << "CO2\t\n";
		
		for(unsigned int i = 0; i < inputs->sim_vol_steps[n].size(); i++)
		{
			fout << i*10 << '\t' << inputs->sim_vol_steps[n][i]/(BTPSout*0.01) << '\t'
				 << model_outputs->simulated[i_offset + i] << "\t0\t0\n";
		}
		i_offset += inputs->sim_vol_steps[n].size();

		fout.close();
	}



	return 0;
}