#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include<string>
#include<limits>

#define MACHINE_DS 0.036
#define SF6_NOISE 0.0002   //largest estimate from SNR in Horsley et al. Thorax 2007
#define VOL_NOISE 0.01    //SNR for **volume** from same paper, note this probably depends on flow rate
#define MRI_NOISE_FRAC 0.02    //1/SNR expected for MRI for mean -- ballpark figure
#define NMR_SAMPLES 1000

const int INT_MASKED_VALUE = std::numeric_limits<int>::infinity();  //-1000
const double DOUBLE_MASKED_VALUE = std::numeric_limits<double>::infinity();  //-1000 //used to identify missing data



#endif