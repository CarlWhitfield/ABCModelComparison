#ifndef READ_WRITE_SETTINGS_H
#define READ_WRITE_SETTINGS_H

const char breath_cutoff_option = 'l'; //'l' = LCI, 'b' = no. of breaths
const double LCI_frac = 0.025;   //fraction of Cinit where LCI point is defined
const int Nbreaths_keep = 6;  //if using option 'b', this is no. of breaths to keep

//threshold conc for observations dureing inhalation
// above which determines washout has definitely begun
const double washin_min_conc = 0.15;
//equivalent threshold below which determines washout has begun
const double washout_max_conc = 0.01;

//options
const bool MEASURE_SUBSET = false;   //if false all breaths are used
const size_t Nmeasure = 20;   //if true, number of breaths to use in measurements
const size_t NphaseII = 4;   // 4 //in measured breath, number of phase II points to use
const size_t NphaseIII = 4;  // 4 //in measured breath, number of phase III points to use

//FIXED STEP SIZE METHOD NEEDS CHECKING
const bool FIXED_STEP_SIZE = false;   //use fixed vol step size in sims?
const double vol_sim_step_frac = 0.04;   //if true, simulation step size in frac of VT mean

#endif
