#ifndef __PROJECT_HCVCALCULATOR_H__
#define __PROJECT_HCVCALCULATOR_H__

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


class Project_HCVCalculator{

public:
	Project_HCVCalculator();
	~Project_HCVCalculator(){};


	int SetInitialAge(double argAge){_initialAge=argAge; return 0;}
	int RunBaseCase();
	
	int RunPerturbations(string argPertValFile, int idxPertub);
	//int BaseCase();

	

	//int CEA_OneWay(int argCmpIdx);
	//int CEA_PSA(int argCmpIdx, int argBatch=0, int argBatchSize=5000);

	//int RunDifferentMCReplc();
	//int PSA();

	//void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}

private:
	int ReadPerturbedValues(string argPertValFile);
	
	int CompareAllCases();
	int Compare(string strArm1, string strArm2, int argGenotype);
	vector<double> EvaluateOneArm(string strArm, const baseCohortType & testCohort);
	vector<double> EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort );
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	double GetValue(string varName);


	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%
	ofstream _outFileSummary;

	vector<string> _listCmp, _listArm1, _listArm2;

	modelParamType _modelParam, _baselineModelParam;

	string _strOutputFile;
	string _strTxExperience;
	double _initialAge;
	
	vector<string> _vecVarOfInterest;
	vector<string> _vecPertVarName;
	vector<double> _vecPertVarVal;
};




#endif