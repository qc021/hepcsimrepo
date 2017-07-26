#ifndef __PROJECT_AGENOTX_H__
#define __PROJECT_AGENOTX_H__

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<cassert> //header file for function "assert()"
#include"sim.h"
using namespace std;


class Project_AgeNoTx{
public:
	Project_AgeNoTx();
	~Project_AgeNoTx(){};

	int BaseCase_detailed();	// need to calculate age threshold outside this program
	int BaseCase_aggregated(double argAgeStart=30.0, double argAgeEnd=86.0, double argAgeIncr=2.0);	// [Preferred] calculate age threshold inside this program
	int BaseCase_distributed(double argAgeStart,double argAgeEnd, double argAgeIncr);



	int ReduceDrugCost(double reductPerct);
	inline void SetInputScenarioFile(string v){FILE_scenario_input=v;}
	inline void SetOutputDetail(bool v){	_outputDetail=v; };

	vector<double> EstimateICER_age_fib(double argAge, stateType argFib);
	vector<double> EstimateICER_age(double argAge);	// for F0, F1, F2, F3, F4, F0-F2, F0-F3, F3-F4, respectively
	vector<double> GetCutoffAge(double argAgeStart, double argAgeEnd, double argAgeIncr, bool argOutput=false);

	int PSA_age_threshold(int argNumPSARunsPerBatch, int idx, double argAgeStart=30.0, double argAgeEnd=86.0, double argAgeIncr=2.0);
	int PSA_ICER_for_each_age(double argAge,int argBatchIdx=-1, long argSeed=-1);
	
	int DSA_OneWay_ICER_for_age(double argAge);

private:
	vector<double> CompareTxAndNoTx(string argTxArmName, const baseCohortType & argCohort);
	int DSA_OneWay_ChangeValue(string varName, double argVal);

	double GetWeight(string argTxExp, string argIFNTol, int argGenotype, char argGender);
	ofstream _outf, _outf_detail, _psaOutf_icer_conf, _psaOutf_threshold;
	bool _outputDetail;

	vector<double> _fibDistr;
	vector<double> _genoDistr;
	double _pMale, _pTE, _pIFN_intol;
	double _p8Week;

	
	string FILE_scenario_input; // default ="project_ageNoTx_input_scenarios.txt";
	modelParamType _modelParam, _modelParam_backup;
	double _drugCostDisc;
};



#endif