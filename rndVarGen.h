#ifndef _RNDVARGEN_H
#define _RNDVARGEN_H

#include <iostream>
#include <random>
#include "SysUtil.h"

using namespace std;

/* Sample code for demonstration of this random variate generator

inline int TestRndVarGnt() {
	mt19937 mt(2);
	//for(int k=0; k<20; k++) cout<<mt()<<endl;

	//cout<<endl<<endl;
	CRndVarGen rvGen1,rvGen2;
	seed_type seed1 = mt();
	seed_type seed2 = mt();
	rvGen1.Initialize(seed1,UNIFREAL,0,1);
	rvGen2.Initialize(seed2, UNIFREAL, 0, 1);
	cout << seed1 << "\t" << seed2 << endl << endl;

	for (int k = 0; k < 100; k++)
		cout << rvGen1.GetU01() << "\t" << rvGen2.GetU01() << endl;

	vector<string> label;
	vector<int> count;
	for(int k=0; k<10; k++){
		label.push_back(basicToStr(k*2));
		count.push_back(k*3);
	}



	cout<<endl;
	//for(int k=0; k<100; k++) cout<<DiscreteDistrSampler(label,count,rvGen)<<endl;
	//for(int k=0; k<20; k++) cout<<rvGen.GetRV()<<endl;
	cout<<endl;

	return 0;
}
*/


typedef unsigned long seed_type;
enum TYPE_RANDOM_DISTR {U01, NORMAL, EXP, GAMMA, BETA, UNIFREAL, EMPCDF, CONSTANT, DISCRETE, UNIFINT};

class CRndVarGen {
public:
	CRndVarGen(){
		for(int k=0;k<3;k++) _distrParam.push_back(0);
		//ctor
	};

	~CRndVarGen(){};

    int Initialize(seed_type argSeed, TYPE_RANDOM_DISTR argDist, double arg1=0, double arg2=0, double arg3=0);

    inline void SetSeed(seed_type  arg){_eng.seed(arg);}

    inline double GetU01(){return distrU01(_eng);}
	inline int GetUnifInt(int argL, int argU) {
		uniform_int_distribution<int> distrInt(argL, argU);
		return distrInt(_eng);		
	}

	inline int GetUnifInt_Rounded(int argL, int argU) {
		int val = (int)floor(distrU01(_eng) * (argU - argL + 1.0 - EPSILON) + (double)argL);
		if (val > argU) {
			ExitWithMsg("Error @ GetUnifInt()");
		}
		return val;
	}


	inline double GetExp(double lambda) {
		exponential_distribution<double> distrExp(lambda);
		return distrExp(_eng);
	}

	inline int GetExp_Integer(double lambda) {
		exponential_distribution<double> distrExp(lambda);
		return (int)distrExp(_eng);
	}


	template <class RandomAccessIterator>
	void ShuffleVector(RandomAccessIterator first, RandomAccessIterator last);

    double GetRV();

protected:

private:
    mt19937 _eng; // or any other generate, see http://www.cplusplus.com/reference/random/
    uniform_real_distribution<double> distrU01;
    TYPE_RANDOM_DISTR _distrType;
    vector<double> _distrParam;

};

template <class RandomAccessIterator>
void CRndVarGen::ShuffleVector(RandomAccessIterator first, RandomAccessIterator last)
{
	// see reference here: http://www.cplusplus.com/reference/algorithm/shuffle/	
	for (int i = (last - first) - 1; i>0; --i) {
		int val = (int)floor(distrU01(_eng) * (i + 1.0 - EPSILON));
		if (val > i) {
			ExitWithMsg("Error @ Shuffle vector()");
		}
		std::swap(first[i], first[val]);
	}
}


// T2 can only be int/double
template<typename T1, typename T2>
inline T1 DiscreteDistrSampler(const vector<T1> & argListLabel, const vector<T2> & argListCount, CRndVarGen & argRnd){
	double u=argRnd.GetU01();
	// rescale
	vector<double> tempvec;
	double s=(double)(Sum(argListCount));
	double cumsum=0;
	int smpIdx=argListCount.size();
	for(int k=0; k<argListCount.size(); k++){
		cumsum=cumsum+((double)argListCount[k])/s;

		if(cumsum>u){
			smpIdx=k;
			break;
		}
	}
	if (smpIdx == argListCount.size()) {
		ExitWithMsg("Error @ DiscreteDistrSampler(): smpIdx == listCount.size()");
	}

	return argListLabel[smpIdx];

	//return 0;
}


// T2 can only be int/double
template<typename T2>
inline int DiscreteDistrSampler(const vector<T2> & argListCount, CRndVarGen & argRnd) {
	double u = argRnd.GetU01();
	// rescale
	vector<double> tempvec;
	double s = (double)(Sum(argListCount));
	double cumsum = 0;
	int smpIdx = argListCount.size();
	for (int k = 0; k<argListCount.size(); k++) {
		cumsum = cumsum + ((double)argListCount[k]) / s;

		if (cumsum>u) {
			smpIdx = k;
			break;
		}
	}
	assert(smpIdx != argListCount.size());

	return smpIdx;

	//return 0;
}

// T2 can only be int/double
template<typename T1>
T1 DiscreteDistrSampler_DistrTable(const map<T1, double> & argTable, CRndVarGen & argRnd) 
{
	double u = argRnd.GetU01();
	
	T1 smpIdx;
	double cumsum = 0;
	
	for (typename map<T1, double>::const_iterator it = argTable.begin(); it != argTable.end(); it++) {
		cumsum = cumsum + it->second;

		if (cumsum>u) {
			smpIdx = it->first;
			break;
		}
	}
	//assert(it != argTable.end());

	return smpIdx;

	//return 0;
}

#endif
