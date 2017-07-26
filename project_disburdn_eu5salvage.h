#ifndef __PROJECT_DISBURD_EU5_SAVALGE_H__
#define __PROJECT_DISBURD_EU5_SAVALGE_H__

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<map>
#include<cassert> //header file for function "assert()"

#include"sim_disburdn.h"

using namespace std;


class Project_SalvageEU5{



public:
	Project_SalvageEU5();
	~Project_SalvageEU5(){};

	inline void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}
	inline void SetScreeningOption(typeScreeningScenario argScr) { _modelParam._disBurdnData._option_screenScenario = argScr; }
	

	void SetOutFile(string argStr) { _outFile = argStr; _outFileSummary.open(_outFile); _outFileSummary.close(); _modelParam._disBurdnOutFile = argStr; }; // added 20160830
	int ReadTable_Param_by_Country(string argStr);

	int ReadTable_TreatmentCapacity(string argStr, string argCountry);
	int ReadTable_TransplantCapacity(string argStr, string argCountry);
	int ReadTable_BackgroundMortality(string argStr);
	int ReadTable_BackgroundMortality_5yAgeGroup_EU5(string argStr, string argCountry);
	int ReadTable_NS5AMarketShare(string argStr);
	int ReadTable_ExtendedUsePRandPI(string argStr, string argCountry);
	int ReadTable_Incidence_Num(string argStr, string argCountry);
	int ReadTable_Incidence_AgeDistr(string argStr);
	int ReadTable_Prob_Aware(string argStr);
	int Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap);

	int Initialize_Model_Parameter(string argStrCountry);

	int BaseCaseRun();

	//int CEA_OneWay(int argCmpIdx, string argFileOnewayParam);
	//int CEA_PSA(int argIdx, int argNBatches);





private:



	double GetValue_byCountry(string argVarName, string argStrCountry);
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	

	
	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%

	ofstream _outFileSummary;
	string _outFile; // added 20160830

	modelParamType _modelParam, _baselineModelParam;

	map<string, map<string, double>> _param_by_country;


};




#endif