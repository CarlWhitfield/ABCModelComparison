# ABCModelComparison

The code in this repository uses ABC-SMC (Approximate Bayesian Computation - Sequential Monte Carlo) to fit models of ventilation in the human lung to MBW (Multiple Breath Washout) data. 

## External Dependencies

The following libraries are required to build the executables:

- Eigen C++ Libraries -- https://eigen.tuxfamily.org/
- Boost C++ Libraries: (filesystem.hpp and random.hpp) -- https://www.boost.org/
- Microsoft MPI
- Header libraries for lung modelling -- https://github.com/CarlWhitfield/Tube-network-headers/

The repository is separated into 3 different executables.

## GenerateMBWData.exe

Takes a '.params' file (see example in 'example_input_files') as command-line input in order to produce a single realisation of the ventilation model and generate simulated MBW data (3 separate '.txt' files).

### Internal Include Files

- ABC_model_selection.h
- MBW_models.h
- MBW_model_params.h
- read_mbw_data.h
- read_write_codes.h
- read_write_settings.h
- compartmental.h

### Source Files

- ABC_model_selection.cpp
- MBW_models.cpp
- MBW_model_params.cpp
- read_mbw_data.cpp
- compartmental.cpp
- GenerateMBWData/main.cpp

### Command Line Execution

The only argument required is the path to the .params file

.\GenerateMBWData.exe Working_dir\test.params

Output will be written to the same directory as the .params file.

## ProcessMBWData.exe

Takes raw MBW data ('.txt' files) and extracts the inputs for the ABC-SMC algorithm. 

### Internal Include Files

- read_mbw_data.h
- read_write_codes.h
- read_write_settings.h

### Source Files

- read_mbw_data.cpp
- compartmental.cpp
- read_write_MBW_files_main.cpp

### Command Line Execution

Several command line arguments are required:

1) (Integer) the number of MBW test files.
2) The paths to the input files
3) CO2 delay in seconds
4) Bag volume in L used for MRI image
5) MBW Machine deadspace in L

E.g. .\ProcessMBWData.exe 3 MBW_test1.txt MBW_test2.txt MBW_test3.txt 0.738 1.0 0.036

Produces a processed MBWTest file for each input file and one summary file 

## ABCModelComparison.exe

Runs the main ABC-SMC algorithm, with settings determined in 'ABCsettings.h'.

### Internal Include Files

- ABC_model_selection.h
- ABCsettings.h
- MBW_models.h
- MBW_model_params.h
- read_write_codes.h
- compartmental.h

### Source Files

- ABC_model_selection.cpp
- MBW_models.cpp
- MBW_model_params.cpp
- compartmental.cpp
- ABCModelComparison/main.cpp

### Command Line Execution

The only command line argument required is the path of the directory containing the files produced by ProcessMBWData.exe. There should only be one set of these files in the directory. 
