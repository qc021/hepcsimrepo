#ifndef __PROJECT_ACUTE_H__
#define __PROJECT_ACUTE_H__

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


class Project_Acute{



public:
	Project_Acute();
	~Project_Acute(){};


	int CEA_BaseResults();
	int CEA_VaryTimeHorizon(int argIdx);
	int CEA_VaryAge(int argIdx);
	int CEA_VaryAcuteTxDuration_probSelfClearance(int argIdx);
	int CEA_VaryAcuteTxDuration_SVR(int argIdx);
	int CEA_VaryTxDuration_Acute_ChronicF0(int argIdx);
	int CEA_Vary_TxCost_SVRacute(int argIdx);
	int CEA_Vary_TxCost_ClearanceRate(int argIdx);
	int CEA_Vary_TxCost_SVRacute_ClrRate(int argIdx);
	int CEA_Vary_AccessProbability(int argIdx);

	int CEA_OneWay(int argCmpIdx, string argFileOnewayParam, int argIdxForBatch=0);
	int CEA_PSA(int argIdx, int argNBatches);

	int RunDifferentMCReplc();
	int PSA();

	void SetMaleRatio(double argV) { _mfDistr.push_back(argV); _mfDistr.push_back(1 - argV); }
	void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}
	void SetOutFile(string argStr) { _outFile = argStr; _outFileSummary.open(_outFile); _outFileSummary.close(); }; // added 20160830
	void SetInitialAge(double argD) { _age_initial = argD; }

	void SetTxDur_Chronic_Acute(int argChronic, int argAcute) { _nWks_f0_tx = argChronic; _nWks_acute_tx = argAcute; }

private:
	vector<double> Compare(string strArm1, string strArm2, int argGenotype);
	vector<double> EvaluateOneArm(string strArm, const baseCohortType & testCohort);
	vector<double> EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort );
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	int ReadPSASampledValues(ifstream & inf, modelParamType & argModelParam);
	vector<double> GetAggregatedResults();

	double _drugCostReduction, _drugCostReduction_bak ; // e.g., -0.25 ==> reduce by 25%

	ofstream _outFileSummary;
	string _outFile; // added 20160830

	vector<string> _listCmp, _listArm1, _listArm2;
	vector<int> _listGenotype;

	modelParamType _modelParam, _baselineModelParam;

	double _age_initial;
	vector<double> _mfDistr;
	int _nWks_acute_tx, _nWks_f0_tx;

	string _note;
};




#endif