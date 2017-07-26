#ifndef __PROJECT_3WK_H__
#define __PROJECT_3WK_H__

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


class Project_3Wk{

public:
	Project_3Wk();
	~Project_3Wk(){};

	int BaseCase();
	int CEA_BaseResults();
	int CEA_OneWay(int argCmpIdx);
	int CEA_PSA(int argCmpIdx, int argBatch=0, int argBatchSize=5000);

	int RunDifferentMCReplc();
	int PSA();
	void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}

private:
	int Compare(string strArm1, string strArm2, int argGenotype);
	vector<double> EvaluateOneArm(string strArm, const baseCohortType & testCohort);
	vector<double> EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort );
	vector<double> EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort );
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	int ReadPSASampledValues(ifstream & inf, modelParamType & argModelParam);
	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%
	ofstream _outFileSummary;

	vector<string> _listCmp, _listArm1, _listArm2;
	vector<int> _listGenotype;

	modelParamType _modelParam, _baselineModelParam;

	
};




#endif