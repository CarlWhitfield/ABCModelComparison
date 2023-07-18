#ifndef SETTINGS_H
#define SETTINGS_H

//ABC settings
#define N_PROCS	1 //number of processors to use
#define N_GENS 15   //0 to just keep going
#define N_ACCEPT 10 //1120 //1400// 1000  //min number of accepted
#define N_CUT 4   //670 //800 //500 //fraction to use as cutoff
#define MIN_DIST 0.0 //cutoff distance, define distance so that this is reasonable
#define USE_DIST_WEIGHTS false    //whether or not to use distance weights
#define USE_SCALE_WEIGHTS false   //whether or not to scale breath weight by distance
#define KS_STAT_CUTOFF 1.00  //normalised K-S statistic (K-S * sqrt(n*m/(n+m)))
#define N_TIMES_CUTOFF 100

#include"MBW_models.h"   //include headers for model libraries

//define base class for model generator
typedef CompartmentalModelGeneratorBase ModelGenBaseClassName;

//define base class for models
typedef CompartmentalModelBase ModelBaseClassName;

//define distance function to use
typedef RMSEDistance DistanceName;

//define model inputs class to use 
typedef MBWModelInputs ModelInputName;

//define model inputs class to use 
typedef MBWModelOutputs ModelOutputName;

//define models to compare (up to 5)
#define N_MODELS 1
typedef BasicLognormalModelGenerator ModelGen1Class;
//typedef BasicLognormalAsyncModelGenerator ModelGen2Class;
//typedef LognormalModelSDSGenerator ModelGen2Class;
//typedef BasicLognormalAsymmModelGenerator ModelGen3Class;
//typedef LognormalAsymmSDSModelGenerator ModelGen4Class;
//typedef BimodalModelSDSGenerator ModelGen2Class;
//typedef BimodalAsymmSDSModelGenerator ModelGen4Class;

#endif