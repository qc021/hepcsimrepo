#ifndef __SIM_H__
#define __SIM_H__

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<cassert>

#include"sim_config.h"
#include"data_type.h"
#include"rand_Gen.h"
#include"SysUtil.h"
#include"StatisticUtil.h"
#include"trials.h"
#include"rndVarGen.h"	// Qiushi's random number generator

using namespace std;




class HepCSim{
public:
	HepCSim();
	~HepCSim(){_tempOutf.close();};
	inline void SetRandomSeed(long argVal){_randomSeed=argVal;}

	int Run(modelParamType  argModelParam);
	int Run_HealthyPopulation(modelParamType  argModelParam);


	int OutputBaseCase(ofstream & outf);
	int ReadSmpParamForPSA(ifstream & inf, modelParamType & argModelParam);
	int PSA_initializeSampler(const vector<string> & vecVarName, const vector<string> & vecDistrName, const vector<double> & vecP1, const vector<double> & vecP2);
	int PSA_initializeSampler(const psaDistrType & argPSADistr);
	int PSA_sampleModelParamValue(modelParamType & argModelParam);
	int PSA_sampleModelParamValue_Univariate(modelParamType & argModelParam, string argVarName);
	int PSA_sampleCohort(baseCohortType & argCohort);	// only change age/state/gender value

	inline double GetAvgQALY(){return _sim_QALY;}
	inline double GetAvgLY() { return _sim_LY; }
	inline double GetAvgCost(){return _sim_cost;}
	inline counterType GetCounter(){return _simCounter;}
	inline double GetTxCost()	{return _sim_cost_tx;}
	inline double GetTestCost() { return _sim_cost_test; }
	inline double GetDALY_YLL() { return _sim_DALY_YLL; }
	inline double GetDALY_YLD() { return _sim_DALY_YLD; }
	inline double GetDALY_Total() { return _sim_DALY_YLL+_sim_DALY_YLD; }

	inline void SetTxDelayedWeeks(int argV){_tx_delayed_weeks=argV;}
	inline void SetPrintResultForEveryPatient(bool v){_print_results_every_patient=v;}

private:
	int SampleRandomNumbers_Tx(double * rand_pr0_CRN, double * rand_prContinueTx_CRN, double & rand_prSVR, double & rand_prETR, double & rand_prAE_Anemia,
		double & rand_prEpoUse, double & rand_prTxFail, double & rand_ETRgivenTxDiscontinue, double & rand_prTxResponse);

	int Run_NaturalHistory(patientType & pat, int & cycleNum, double & rand_pr, double & pr_Death, double * pr_mDeathAllCause, double * pr_fDeathAllCause,
		double & q_AgeNormal, double * q_mNormal, double * q_fNormal, modelParamType & argModelParam, double & eQALY, double & eLY, double & eCost, double & discountFactAdjQ, double & discountFactAdjC,
		int delayedWeek=-1, int horz = -1);

	bool Run_AcutePhaseHCV(patientType & argPat, int & argCycleNum, modelParamType & argModelParam, double & eQALY, double & eLY, double & eCost, 
		double * pr_mDeathAllCause, double * pr_fDeathAllCause, double * q_mNormal, double * q_fNormal, double & discountFactAdjQ, double & discountFactAdjC, double & c_drug);

	
	int ConvertParameters(modelParamType & argModelParam);
	int ResetCounters();
	
	long _randomSeed;	// =0 by default

	counterType _simCounter;
	double _sim_QALY, _sim_cost, _sim_cost_tx, _sim_cost_test, _sim_LY;

	double _sim_DALY_YLD, _sim_DALY_YLL;	// added for DALY results @ 8/26/2016

	ofstream _tempOutf;
	map<string, CRndVarGen> _psa_sampler;
	CRndVarGen _psa_sampler_age, _psa_sampler_gender, _psa_sampler_state;
	CRndVarGen _rnd_extra;

	mt19937 _psa_seed;
	
	int _tx_delayed_weeks;
	bool _print_results_every_patient;



};


#endif