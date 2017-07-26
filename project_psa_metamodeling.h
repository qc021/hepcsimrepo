#ifndef __PROJECT_PSA_METAMODELING_H__
#define __PROJECT_PSA_METAMODELING_H__

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


class Project_PSA_MetaModeling{
public:
	Project_PSA_MetaModeling() { _doNotSync = false; _drugCostReduction = 0; _flagPrintParameterValues = false; };
	~Project_PSA_MetaModeling(){};
	
	int RunBaseCase();
	int RunPSA_Multivariate_OutputSingleMCRun(long argNumPSAIter, int argNumRplc);
	int RunParameterSmplingTest(int seed, int num);
	int RunModelWithDiffParamValue();

	inline void SetDrugCostReduction(double reductPerct) { _drugCostReduction = reductPerct; }
	inline void SetPrintParamValues(bool argV) { _flagPrintParameterValues = argV; }
	int ChangeValue(string varName, double argVal);

private:	

	int Compare(int argIdxMC=0, int argIdxPSA=0);

	vector<double> EvaluateOneArm(string strArm, const baseCohortType & testCohort);

	int ReadDistrForPSA(string argFileStr);

	void PrintHeader(ofstream & outf);

	vector<string> _distrPSA_varName;
	vector<string> _distrPSA_distrName;
	vector<double> _distrPSA_param1, _distrPSA_param2;

	mt19937 _seed_gnt;
	bool _doNotSync;

	ofstream _outf;

	ofstream _outFileSummary;
	double _drugCostReduction;
	modelParamType _modelParam, _baselineModelParam;
	int _seedOfOneMCRun;
	bool _flagPrintParameterValues;
};



#endif