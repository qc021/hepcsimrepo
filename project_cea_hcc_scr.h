#ifndef __PROJECT_CEA_HCC_SCR_H__
#define __PROJECT_CEA_HCC_SCR_H__

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<cassert> //header file for function "assert()"
#include"sim_screening.h"
using namespace std;


class Project_CEA_HCC_Screening{


public:
	Project_CEA_HCC_Screening();
	~Project_CEA_HCC_Screening(){};

	void SetDrugCostReduction(double reductPerct) { _drugCostReduction = reductPerct; }
	void SetOutFile(string argStr) { _outFile = argStr; _outFileSummary.open(_outFile); 
	if (_outFileSummary.fail()) {
		ExitWithMsg("Error @ cannot open file: " + _outFile);
	}
	_outFileSummary.close();
	};

	int Initialize();

	int CEA_BaseResults();	
	int CEA_OneWay(int argCmpIdx, string argFileOnewayParam, int argIdxForBatch=0);
	int CEA_PSA(int argIdx, int argNBatches);	
	

private:
	int ReadTable_InputParam();
	int ReadTable_BackgroundMortality(string argStr);
	double GetValue(string argVarName);
	vector<double> EvaluateScreeningPolicy(int argScrInterval);
	
	int SetInitialPopulation(type_scr_fib_state argFibState);
	int SetInitialPopulationMixed() { _modelParam._hccScreeningData._distr_initial_fib = _param_map_distrInitialFib; return 0.0; }

	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	
	
	double _drugCostReduction, _drugCostReduction_bak ; // e.g., -0.25 ==> reduce by 25%

	ofstream _outFileSummary;
	string _outFile; // added 20160830

	map<string, double> _param_map;
	map<string, vector<double>> _param_map_txDistr;
	map<type_scr_fib_state, double> _param_map_distrInitialFib;
	modelParamType _modelParam, _baselineModelParam;

	string _note;
};




#endif