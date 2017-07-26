#ifndef __PROJECT_DISBURD_PAKISTAN_H__
#define __PROJECT_DISBURD_PAKISTAN_H__

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

enum type_pakistan_incidence { INCI_CONSTANT, INCI_PROPORTIONAL, INCI_INCREASING };

class Project_DisBurd_Pakistan{

public:
	
	Project_DisBurd_Pakistan();
	~Project_DisBurd_Pakistan(){};

	inline void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}
	inline void SetScreeningOption(typeScreeningScenario argScr) { _modelParam._disBurdnData._option_screenScenario = argScr; }
	

	void SetOutFile(string argStr) { _outFile = argStr; _outFileSummary.open(_outFile); _outFileSummary.close(); _modelParam._disBurdnOutFile = argStr; }; // added 20160830
	int ReadTable_Param();
	double GetValue(string argVarName);

	//int ReadTable_TreatmentCapacity(string argStr, string argCountry);
	int ReadTable_BackgroundMortality(string argStr);
	int ReadTable_BackgroundMortality_5yAgeGroup(string argStr, string argCountry);

	//int ReadTable_NS5AMarketShare(string argStr);
	//int ReadTable_ExtendedUsePRandPI(string argStr, string argCountry);
	//int ReadTable_Incidence_Num(string argStr, string argCountry);
	//int ReadTable_Incidence_AgeDistr(string argStr);
	
	int ReadTable_Prob_Aware(string argStr);
	int Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap);

	int Initialize_Pakistan_Parameter(type_pakistan_incidence argInciMode = INCI_PROPORTIONAL);

	int BaseCaseRun();
	int TestWHOTarget(int argIdx=0);
	int CalibInitialSizeAndInciRate(int argIdx =0);

	//int CEA_OneWay(int argCmpIdx, string argFileOnewayParam);
	//int CEA_PSA(int argIdx, int argNBatches);



	void Set_param_pk_tx_cap_future(int val) { _paramPK_tx_cap_future = val; }

private:



	//double GetValue_byCountry(string argVarName, string argStrCountry);
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	

	
	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%

	ofstream _outFileSummary;
	string _outFile; // added 20160830

	modelParamType _modelParam, _baselineModelParam;

	map<string, map<string, double>> _param_by_country;


	// ---------- define additional parameters for pakistan ----------

	int _paramPK_tx_cap_future;
	double _paramPK_universal_screeningRate, _paramPK_universal_screeningCap;
	type_pakistan_incidence _paramPK_inci_mode;
	int _paramPK_inci_const;
	double _paramPK_inci_proportion;
	map<string, double> _paramPK_map;
};




#endif