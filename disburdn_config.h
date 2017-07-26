#ifndef __DISBURDN_CONFIG_H__
#define __DISBURDN_CONFIG_H__


#define __DEFAULT_SETTING__	// for US, EU5 analysis
//#define __PAKISTAN_SETTING__	// for Pakistan analysis

#include <iostream>
#include "sim_config.h"
#include "SysUtil.h"

using namespace std;



const bool RECOVER_SALVAGE_POSTER_VERSION = false;
const bool FIXED_CAP_BIRTH_COHORT_SCREENING = false;

const int SCREEN_BIRTH1_YR = 1945;
const int SCREEN_BIRTH2_YR = 1965;

const int SCREEN_BIRTHCOHORT_START_YR = 2013;
const int SCREEN_BIRTHCOHORT_END_YR = 2019;




#if defined(__DEFAULT_SETTING__)
// default value:
const int START_YR = 1995;
const int END_YR = 2040; 
const long DISBDEN_INITIAL_NUM_PATIENTS = 4000000;
//const long DISBDEN_INITIAL_NUM_PATIENTS = 40000;
const int START_YR_SCREEN = 2013;
const int TX_START_YR = 2001;

const int START_YR_PEGRBV = 2001;
const int START_YR_PI1 = 2012;
const int START_YR_DAA1 = 2014;
const int START_YR_DAA2 = 2015;
const int START_YR_DAA3 = 2018;

//const int START_YR_DAA1 = 2040;
//const int START_YR_DAA2 = 2040;
//const int START_YR_DAA3 = 2040;

const int START_YR_TX_PRIORITIZATION = START_YR_DAA1; // START_YR_PI1;
const int END_YR_TX_PRIORITIZATION = 2030;
const int YEAR_RECORD_AGE_DISTR = 2001;


#elif defined(__PAKISTAN_SETTING__)

// values for pakistan disease burden analysis
const int START_YR = 1995;
const int END_YR = 2040;
const long DISBDEN_INITIAL_NUM_PATIENTS = 10000000;

const int START_YR_SCREEN = 2018;
const int TX_START_YR = 2004;

const int START_YR_PEGRBV = 1998;
const int START_YR_PI1 = 2016;
const int START_YR_DAA1 = 2016;
const int START_YR_DAA2 = 2016;
const int START_YR_DAA3 = 2016;
const int START_YR_TX_PRIORITIZATION = START_YR_DAA3; // START_YR_PI1;
const int END_YR_TX_PRIORITIZATION = END_YR;
const int YEAR_RECORD_AGE_DISTR = 2008;
#endif



const int AGE_CUTOFF_MEDICARE = 65;
const int ACA_START_YR = 2014;

const int MAX_NUM_TX_PEGRBV = 2;
const int MAX_NUM_TX_PI = 1; 
const int MAX_NUM_TX_DAA = 3;
const int MAX_NUM_TX_PRETRIPLE_contraind_mod = 1;
const int MAX_NUM_TX_CONTRA_FOR_PI1 = 2;
// optional; need to check with Jag
const int MAX_NUM_TX_ELIGIBLE_FOR_WAITING_PREWAVE1 = 2;

//const int MAX_NUM_TX_PRETRIPLE = 2;
//const int MAX_NUM_TX_PRETRIPLE_contraind_mod = 1;
//const int MAX_NUM_TX_TRIPLE = 3;
//const int MAX_NUM_TX_TRIPLE_contraind_mod = 2;
//const int MAX_NUM_TX_WAVE1_WAVE2 = 4;
//const int MAX_NUM_DAA_WAVE1_WAVE2 = 3;

const int NUM_SCREEN_BIRTH_COHORT_EACH_YEAR = 66572;


//const int WAVE3_START_YR = 2020;

const double LIV_TRANSP_MAX_AGE = 75;
const int LIV_TRANSP_CAP_DC = 1364;
const int LIV_TRANSP_CAP_HCC = 1216;


const bool ENABLE_HCC_FROM_F3 = false;
const bool ENABLE_HCC_FROM_F3SVR = false;
const bool ENABLE_REGRESSION = false;
const int YEAR_REGRESSION_AFTER_SVR = 2;

extern bool FLAG_SMOOTH_RESULTS;

extern bool FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION;
#endif

