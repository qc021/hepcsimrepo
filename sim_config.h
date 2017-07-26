#ifndef __SIM_CONFIG_H__
#define __SIM_CONFIG_H__
#include <iostream>
#include <string>
using namespace std;


const double CYCLE_LENGTH = 1.0; // cycle length in weeks
const double CYCLES_PER_YEAR = 52.17857143; // number of weeks in 1 year = 52.17857143
extern double DISCOUNT_FACT_C ; //annual discount factor for costs
extern double DISCOUNT_FACT_Q ; //anual discount factor for QALYs
const double DISCOUNT_FACT_LY = 0.0; // annual discount factor for LYs
extern double DISCOUNT_FACT_L;
extern double DISCOUNT_FACT_DALY; //anual discount factor for QALYs


extern bool EXPECTED_LIFE_ONLY ; // making it equal to 1 will calculate expected life years instead of QALYs
const bool ADVERSE_EVENTS = 1; // to exclude adverse event use ADVERSE_EVENTS = 0
 
const int TREATMENT_ARM = 0; //select which arm to run: 0=ALL,
const int USE_RANDOM_SEED = 0; //if 1, random seed is used based on clock time, otherwise, same random sequence will be generated each time.
const int TREATMENT_NAIVE = 0; // use 1 for tx-naive and 0 for tx-experienced
const int PROFILE_BOOTSTRAP = 0; // use 1 to bootsrap patient profiles from trial data

const int PRINT_PATLEVEL_OUT = 0; //use 1 to print individual patient level output
const int PRINT_RAWOUTPUT = 0 ; //use 1 to print raw output files for WebModel / ExcelModel
extern int PRINT_SCREEN; // use 1 to print output on screen
extern int PROGRESSION_AGE_GENDER ; //use 1 for fibrosis progression rates to vary with age and gender
const double HEALTH_INFLATION =  423.815/400.258 ;//using medical care component of CPI from 2011 to 2013 April
extern double TIME_HORIZON; //150.0;//30.0;//use 150.0 for lifetime or as defined [default=150]

// [QC] Added by Qiushi, 2/24/2015
const bool USE_COMMON_RANDOM_NUMBER = true;

// [QC] following variables are not "const" any more
extern bool TREATMENT; // to swicth-off treatment, change value to 0
extern long SIM_SEED;
extern long NUM_PATIENTS;

// [QC] Added by Qiushi, 9/5/2014
const int MAX_TXWEEK = 73;
const int OPTION_CREATE_PATIENT_PROFILE = 3;	//	1. [may be problematic] Use a distribution observed in trial to generate age,gender, fibrosis state
												//	2. Bootstrap from the cohort used in the trial 
												//	3. [base case] Use a single base profile 
// [QC] moved from trial-evaluation-function:
const double calibrationFactor = 1.0;
const double maxLivTrAge = 120.0;



extern int PSA_OPTION;			//use 1 to run PSA option and 0 for not running PSA
extern int ONEWAYSA_OPTION_PROJ_ACUTE;		//use 1 to run Oneway SA option and 0 for not running oneway SA

// [QC] added by Qiushi, 9/20/2014
extern long NUM_RUNS_PSA; // = 10000;//5000; // sample size for PSA
extern unsigned long PSA_SEED; // = 1;
const bool PSA_SAME_SEED_FOR_ALL_PARAM = false;
const bool PSA_SAME_QOL_F0_F3 = true;
const long PSA_MAX_SAMPLE_IN_SINGLE_OUTPUTFILE = 100000;
extern long NUM_RUNS_PSA_ONE_BATCH;


const string FILE_LHS = "LHSinput.in";	//  "LHSinput1000.in"
										//	"LHSinputOld.in"

const string FILE_base_cohort_profile = "PatientBaseProfilesP05101.in";
										// "PatientBaseProfiles.in"
										// "PatientBaseProfilesP0501F0.in"
										// "PatientBaseProfilesP05216.in"
										// "PatientBaseProfilesP05216AA.in"
										// "PatientBaseProfilesP05216nonAA.in"

extern string FILE_background_mortality;
//const string FILE_background_mortality = "Input_mortality_male_female_India.in";
//const string FILE_background_mortality = "Input_mortality_male_female.in";

extern bool ALLOW_RETREATMENT;
extern bool OVERRIDE_NOTREATMENT;
extern bool OVERRIDE_SYNC_RND_STRM;
const int MAX_LINES_TREATMENT = 2;

#endif