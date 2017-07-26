#ifndef __SIM_SCREENING_H__
#define __SIM_SCREENING_H__

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
#include"SysUtil.h"
#include"StatisticUtil.h"
#include"trials.h"
#include"rndVarGen.h"	// Qiushi's random number generator

using namespace std;



class HCCScreenSim {
public:
	HCCScreenSim();
	~HCCScreenSim() { _tempOutf.close(); };
	inline void SetRandomSeed(long argVal) { _randomSeed = argVal; }

	void SetScreeningInterval(int argInt) { _scr_interval_months = argInt; }
	int Run(modelParamType  argModelParam);





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
	inline double GetTestCost() { return _sim_cost_test; }
	//inline double GetDALY_YLL() { return _sim_DALY_YLL; }
	//inline double GetDALY_YLD() { return _sim_DALY_YLD; }
	//inline double GetDALY_Total() { return _sim_DALY_YLL+_sim_DALY_YLD; }

	inline void SetTxDelayedWeeks(int argV){_tx_delayed_weeks=argV;}
	inline void SetPrintResultForEveryPatient(bool v){_print_results_every_patient=v;}

private:

	int _scr_interval_months; // default value = 3 months
	

	double _cycles_of_one_year; // model cycle = 1 month

	int InitializePatient(HCCScrPatientType & pat, const modelParamType& argModelParam);

	int StateTransDeath(double pProgHCC, double pProgFib, double bgMort, double livMort, 
		type_scr_fib_state &  state_next_fib, type_scr_hcc_state & state_next_hcc, type_scr_absb_state & state_next_abs);

	type_scr_hcc_state StateTransHCC(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram);

	double GetMortLiv(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram);
	double GetProbHCCProg(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram);
	double GetTxBaselineQoL(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram);
	string GetFibHCCStateStr(const HCCScrPatientType & argPat);

	int Run_NaturalHistory(HCCScrPatientType & argPat, const modelParamType & argModelParam);
	int Run_HCCTreatment(HCCScrPatientType & argPat, const modelParamType & argModelParam);
	
	int ConvertParameters(modelParamType & argModelParam);
	
	int ResetCounters();
	
	long _randomSeed;	// =0 by default

	counterType _simCounter;
	int _sim_cycle;
	
	// cumulative outcomes
	double _cumQALY;
	double _cumLY;
	double _cumCost, _cumCost_screening;

	// average outcomes
	double _sim_QALY, _sim_cost, _sim_LY , _sim_cost_tx, _sim_cost_test;

	//double _sim_DALY_YLD, _sim_DALY_YLL;	// added for DALY results @ 8/26/2016

	ofstream _tempOutf;
	map<string, CRndVarGen> _psa_sampler;
	CRndVarGen _psa_sampler_age, _psa_sampler_gender, _psa_sampler_state;
	CRndVarGen _rnd_pt;
	CRndVarGen _rnd_extra;

	mt19937 _psa_seed;
	
	int _tx_delayed_weeks;
	bool _print_results_every_patient;



};


#endif