#ifndef __PROJECT_GLOBAL_CEA_H__
#define __PROJECT_GLOBAL_CEA_H__

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


class Project_CEA_Global{

public:
	Project_CEA_Global();
	~Project_CEA_Global(){};

	void SetVersionNote(string argStr) { _strVerMark = argStr; }
	void SetOutFile(string argStr){_outFile=argStr; _outFileSummary.open(_outFile); _outFileSummary.close();};
	void SetInitialAge(double argD) { _age_initial = argD; }
	void SetOverridePerfectTreatment(bool argV) { _modelParam._transData._flag_override_pefct_tx = argV; }


	void RunCEAForOneCountry(string argCountry, int argIdx = 0);
	void RunCEAForListedCountries(int argIdx);





	void SetDrugCostReduction( double reductPerct ){_drugCostReduction=reductPerct;}

private:	
	int LoadParamByCountry(string argCountry);

	int CEA_BaseResults();
	int CEA_PSA(int argIdx, int argNBatches);
	int CEA_VaryAge(int argIdx);




	vector<double> GetAggregatedResults();
	vector<double> Compare(string strArm1, string strArm2, int argGenotype);
	vector<double> EvaluateOneArm(string strArm, const baseCohortType & testCohort);
	vector<double> EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort ); // added from HCVCaculator @ 7/9/2016
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);

	double _drugCostReduction; // e.g., -0.25 ==> reduce by 25%
	ofstream _outFileSummary;
	string _outFile;

	
	vector<int> _listGenotype;

	vector<string> _listCountries;
	modelParamType _modelParam, _baselineModelParam, _baselineModelParam_country;


	double _age_initial;
	int _id_batch;
	string _strVerMark;
	string _str_file_temp_life_table;
	map<string, map<string, double>> _map_countryParam;
	map<string, map<int, double>> _map_country_life_table_male, _map_country_life_table_female;
	map<int, double> _map_qol_male, _map_qol_female;

};




#endif