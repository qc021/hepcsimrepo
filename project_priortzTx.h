#ifndef __PROJECT_PRIORTZTX_H__
#define __PROJECT_PRIORTZTX_H__

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


class Project_priortzTx{
public:
	Project_priortzTx();
	~Project_priortzTx(){};

	// return: [0]=Q-treat-now, [1]=Q-delay-1-year, [2]=difference
	vector<double> GetLossInQALYs(double argAge, stateType argFib,txHistryTolType argTxHistrTol,char argGender, string argArmName);	
	vector<double> GetLossInQALYs_Combine8and12Weeks(double argAge, stateType argFib,txHistryTolType argTxHistrTol,char argGender);	

	int BaseCase();
	int OneWaySA(double argAge, int argDelayedWeeks, txHistryTolType argHistryTol, stateType argFib);
private:
	
	string GetArmName(stateType argFib, txHistryTolType argTxHistr);
	int ChangeValue(string varName, double argVal, modelParamType & _modelParam);
	
	ofstream _outf, _outf_detail;

	vector<double> _fibDistr;
	vector<double> _genoDistr;
	double _pMale, _pTE, _pIFN_intol, _p8Week;
	modelParamType _modelParam, _baselineModelParam;

	int _delayedWeeks;

	vector<double> rs_debug0;
	vector<double> rs_debug1;
	vector<double> rs_debug2;
};



#endif