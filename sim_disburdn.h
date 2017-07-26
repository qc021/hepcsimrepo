#ifndef __SIM_DISBURDN_H__
#define __SIM_DISBURDN_H__

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<cassert>
#include<map>
#include<list>


#include"disburdn_config.h"

#include"data_type.h"
#include"rand_Gen.h"
#include"SysUtil.h"
#include"StatisticUtil.h"

#include"rndVarGen.h"	// Qiushi's random number generator

using namespace std;




class BurdenModelSim{
public:
	BurdenModelSim();
	~BurdenModelSim(){_tempOutf.close();};

	int Run(const modelParamType & argModelParam);
	//int OutputCounters(string argFile);
	int PrintOutputHeader(string argFile);
	int PrintOutputRow(int argModelCycle);
	int OutputAgeDistr(const map<double,int> & ageDistr, string argFile);


	int ReadSmpParamForPSA(ifstream & inf, modelParamType & argModelParam);
	int PSA_initializeSampler(const vector<string> & vecVarName, const vector<string> & vecDistrName, const vector<double> & vecP1, const vector<double> & vecP2);
	int PSA_initializeSampler(const psaDistrType & argPSADistr);
	int PSA_sampleModelParamValue(modelParamType & argModelParam);
	int PSA_sampleModelParamValue_Univariate(modelParamType & argModelParam, string argVarName);
	int PSA_sampleCohort(baseCohortType & argCohort);	// only change age/state/gender value

	disBurdnCounterType GetSimCounters() { return _simCounter; }
	void RedfineStartYear(int argV) { _start_year_redefined = argV; }

protected:
	int GeneratePatient(const modelParamType & argModelParam);
	int GeneratePatient_nonNHANES(const modelParamType & argModelParam);
	int GeneratePatient_immigrants_lpr(const modelParamType & argModelParam);

	int ChangeOfCohort(const modelParamType	 & argModelParam);
	int DetmInsrAndAwarenessStatus(patientType & argPt, const modelParamType & argModelParam, bool argIsForInit);
	int DetmRespsState(patientType & argPt, const modelParamType & argModelParam, bool argIsForInit);

	int Screening(const modelParamType & argModelParam);
	int Screening_UsualCare(const modelParamType & argModelParam);
	int Screening_BirthCohort(const modelParamType & argModelParam);
	int Screening_Universal_Rate(const modelParamType & argModelParam);
	int Screening_Universal_Capacity(const modelParamType & argModelParam);

	int Screening_MultipleCohort(const modelParamType & argModelParam);
	int Screening_Universal_Rate(const modelParamType & argModelParam, typeCohortDisBurdnModel argCohort);
	int Screening_BirthCohort_byGroup(const modelParamType & argModelParam, typeBirthCohortGroup argBCGroup, int argScrCap);


	int Treatment_DetermineElig(const modelParamType & argModelParam);
	bool Treatment_DetermineElig(const modelParamType & argModelParam, patientType * thePatientPointer);
	int Treatment_DetermineElig_w_ExtendedUsePRandPI(const modelParamType & argModelParam);
	

	int Treatment_Prioritize(const modelParamType & argModelParam);
	int Treatment_Prioritize_MultipleCohort(const modelParamType & argModelParam);
	int Treatment_Prioritize_byGroup(const modelParamType & argModelParam, typeTxCohortGroup argTxGroup);
	int Treatment_Treat(const modelParamType & argModelParam);

	int RunNaturalHistory_indvPt(patientType & argPt, const modelParamType & argModelParam);
	int UpdateInsurance_indvPt(patientType & argPt, const modelParamType & argModelParam);
	int UpdateCounter_Treamtent(const patientType & argPt, const modelParamType & argModelParam);
	int UpdateCounter_HealthState(const patientType & argPt, const modelParamType & argModelParam);
	int UpdateCounter_Insurance(const patientType & argPt);
	int UpdateCounter_DeathEvent(const patientType & argPt, const modelParamType & argModelParam);
	int UpdateCounter_TxEligibility();
	

	template<typename T> int PrintOutputRow_MapTemplate(const map<T, vector<int>> & argMap, int argIdx);
	template<typename T> int SaveOutputRow_MapTemplate(const map<T, vector<int>> & argMap, int argIdx, vector<double> & argVec, bool argSmooth = false);
	template<typename T> int SaveOutputRow_MapTemplate(const map<T, vector<double> > & argMap, int argIdx, vector<double>& argVec, bool argSmooth=false);

	double CalcDALY_GetYLD(double argAge, typeStateDisBurdnModel argState, double argL, const DALYType & argDALYData);
	double CalcDALY_GetYLL(double argAge, char argGender, const DALYType & argDALYData);
	double CalcDALY_GetDisutilityWeight_DisBurdenModel(typeStateDisBurdnModel argState, const DALYType & argDALYData);

	double CalcCost_Treatment(const patientType & argPt, const costType & argCostData, int argYear);
	double CalcCost_Screening(const costType & argCostData);
	double CalcCost_NaturalHistory(const patientType & argPt, const costType & argCostData);

	
	inline int PrintOutputHeader_MapStringIndex(const map<string, vector<int>> & argMap) {
		for (map<string, vector<int>>::const_iterator it = argMap.begin(); it != argMap.end(); it++) { _outf << it->first << "\t"; }
		return 0;
	};

	inline int PrintOutputHeader_MapStringIndex(const map<string, vector<double>> & argMap) {
		for (map<string, vector<double>>::const_iterator it = argMap.begin(); it != argMap.end(); it++) { _outf << it->first << "\t"; }
		return 0;
	};

	int ResetCounters();
	
	long _randomSeed;	// =0 by default

	disBurdnCounterType _simCounter;
	double _sim_QALY, _sim_cost, _sim_cost_tx;

	double _sim_DALY_YLD, _sim_DALY_YLL;	// added for DALY results @ 8/26/2016

	ofstream _tempOutf;

	map<string, CRndVarGen> _psa_sampler;
	CRndVarGen _psa_sampler_age, _psa_sampler_gender, _psa_sampler_state;
	mt19937 _psa_seed;
	CRndVarGen _rnd_vector_shuffler;
	CRndVarGen _rnd_pt, _rnd_sim;

	list<patientType*> _listAllPts;
	vector<patientType*> _vecPt_BirthCohortScr;
	vector<patientType*> _vecPtTxElig_F0F2, _vecPtTxElig_F3F4, _vecPtTxElig_all, _vecTxCandidates;
	vector<patientType*> _vecPtTxElig_F0F2_nonDefaultCohort, _vecPtTxElig_F3F4_nonDefaultCohort,_vecPtTxElig_all_nonDefaultCohort; // added for nonNHANES, 5/24/2017

	int _start_year_redefined;
	int _curModelCycle;
	int _curYear;
	typeTxCategory _curTxWave;

	long long _cum_num_pt_generated;
	int _num_aware_by_usual_care, _num_aware_by_screening;
	int _num_screening;
	map<typeBirthCohortGroup, int> _scrCap_birthCohort_byGroup;


	vector<int> _counter_livTr_DC, _counter_livTr_HCC;

	string _outf_str;
	ofstream _outf;

	vector<double> _outputRowData, _outputRowData_previousCycle;
};


template<typename T>
int BurdenModelSim::SaveOutputRow_MapTemplate(const map<T, vector<int> > & argMap, int argIdx, vector<double>& argVec, bool argSmooth)
{

		for (typename map<T, vector<int>>::const_iterator it = argMap.begin(); it != argMap.end(); it++) {
			if (argIdx == 0) {
				argVec.push_back(it->second[argIdx]);
			}
			else {
				if (argSmooth) {
					argVec.push_back((it->second[argIdx] + it->second[argIdx - 1]) / 2.0);
				}
				else {
					argVec.push_back(it->second[argIdx]);
				}
			}
		}

		
	return 0;
}

template<typename T>
int BurdenModelSim::SaveOutputRow_MapTemplate(const map<T, vector<double> > & argMap, int argIdx, vector<double>& argVec, bool argSmooth)
{

	for (typename map<T, vector<double>>::const_iterator it = argMap.begin(); it != argMap.end(); it++) {
		if (argIdx == 0) {
			argVec.push_back(it->second[argIdx]);
		}
		else {
			if (argSmooth) {
				argVec.push_back((it->second[argIdx] + it->second[argIdx - 1]) / 2.0);
			}
			else {
				argVec.push_back(it->second[argIdx]);
			}
		}
	}


	return 0;
}


template<typename T>
int BurdenModelSim::PrintOutputRow_MapTemplate(const map<T, vector<int> > & argMap, int argIdx) {
	for (typename map<T, vector<int>>::const_iterator it = argMap.begin(); it != argMap.end(); it++) { _outf << it->second[argIdx] << "\t"; }
	return 0;
};


#endif
