#ifndef __PROJECT_DISBURD_NONNAHANES_H__
#define __PROJECT_DISBURD_NONNAHANES_H__

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

class Project_DisBurd_UScomprh{

	
public:
	
	Project_DisBurd_UScomprh();
	~Project_DisBurd_UScomprh(){};

	inline void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}
	inline void SetScreeningOption(typeScreeningScenario argScr) { _modelParam._disBurdnData._option_screenScenario = argScr; }
	

	void SetOutFile(string argStr) {
		_outFile = argStr; _outFileSummary.open(_outFile); 
		if (_outFileSummary.fail()) { ExitWithMsg("ERROR @ SetOutFile(): Unable to create file " + argStr); }
		_outFileSummary.close(); 
		_modelParam._disBurdnOutFile = argStr;
	}; // added 20160830

	int ReadTable_Param();
	double GetValue_byCohort(string argVarName, const  typeCohortDisBurdnModel & argCohort);

	int ReadTable_immigrants_hcv_cases(string argStr);
	int ReadTable_immigrant_age_distr(string argStr, vector<int> & argAgeCat, vector<double> & argAgeDistr);
	int ReadTable_TreatmentCapacity(string argStr);
	int ReadTable_BackgroundMortality(string argStr);
	int ReadTable_BackgroundMortality_5yAgeGroup_IndianResv(string argStr);

	int ReadTable_NS5AMarketShare(string argStr);
	//int ReadTable_ExtendedUsePRandPI(string argStr, string argCountry);
	int ReadTable_Incidence_Num(string argStr);
	int ReadTable_Incidence_AgeDistr(string argStr);
	
	int ReadTable_Prob_Aware(string argStr);
	int Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap);

	int Initialize_UScomprh_Parameter();

	int BaseCaseRun();




private:

	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	
	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%
	ofstream _outFileSummary;
	string _outFile; // added 20160830

	modelParamType _modelParam, _baselineModelParam;

};




#endif