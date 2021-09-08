# ABCModelComparison

The code in this repository uses ABC-SMC (Approximate Bayesian Computation - Sequential Monte Carlo) to fit models of ventilation in the human lung to MBW (Multiple Breath Washout) data. 

The code required Boost C++ Libraries and OpenMPI

Instructions for building, compiling, and executing to follow. Briefly, the code generates 3 executables:

## GenerateMBWData.exe

Takes a '.params' file (see example in 'example_input_files') as command-line input in order to produce a single realisation of the ventilation model and generate simulated MBW data (3 separate '.txt' files)

## ProcessMBWData.exe

Takes raw MBW data ('.txt' files) and extracts the inputs for the ABC-SMC algorithm. Several command line arguments are required:

1) (Integer) the number of MBW test files.
2) The paths to the input files
3) CO2 delay in seconds
4) Bag volume in L used for MRI image
5) MBW Machine deadspace in L

E.g. .\ProcessMBWData.exe 3 MBW_test1.txt MBW_test2.txt MBW_test3.txt 0.738 1.0 0.036

Produces a processed MBWTest file for each input file and one summary file 

## ABCModelComparison.exe

Runs the main ABC-SMC algorithm, with settings determined in 'ABCsettings.h'. The only command line argument required is the path of the directory containing the files produced by ProcessMBWData.exe. There should only be one set of these files in the directory. 
