#ifndef __PROJECT_INDIA_H__
#define __PROJECT_INDIA_H__

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


class Project_India{

public:
	Project_India();
	~Project_India(){};

	void SetOutFile(string argStr){_outFile=argStr; _outFileSummary.open(_outFile); _outFileSummary.close();};
	void SetInitialAge(double argD) { _age_initial = argD; }
	void SetOverridePerfectTreatment(bool argV) { _modelParam._transData._flag_override_pefct_tx = argV; }

	int CEA_BaseResults();
	int CEA_OneWay(int argCmpIdx, string argFileOnewayParam);
	int CEA_PSA(int argIdx, int argNBatches);
	int CEA_VaryTimeHorizon(int argIdx);
	int CEA_VaryAgeAndTimeHorizon(int argIdx);
	int CEA_VaryAge(int argIdx);
	int CEA_VaryDrugCost(int argIdx);

	int CEA_SA_BackgroundCost();
	int CEA_R1_LiverRelatedMortality();
	//int CEA_PSA(int argCmpIdx, int argBatch=0, int argBatchSize=5000);
	//int PSA();

	void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}
	int CEA_R1_LifeExpectancy_w_wo_HCV();
	int CEA_R1_discount();
	int CEA_R1_progRisk();
	int CEA_R1_progRisk_and_Horizon();
private:
	vector<double> GetAggregatedResults();
	vector<double> Compare(string strArm1, string strArm2, int argGenotype);
	vector<double> EvaluateOneArm(string strArm, const baseCohortType & testCohort);
	vector<double> EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort ); // added from HCVCaculator @ 7/9/2016
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);

	//vector<double> EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort );
	//int ReadPSASampledValues(ifstream & inf, modelParamType & argModelParam);

	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%
	ofstream _outFileSummary;
	string _outFile;

	vector<string> _listCmp, _listArm1, _listArm2;
	vector<int> _listGenotype;

	modelParamType _modelParam, _baselineModelParam;

	double _age_initial;
	


};




#endif