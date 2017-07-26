#include "rndVarGen.h"


int CRndVarGen::Initialize(seed_type argSeed, TYPE_RANDOM_DISTR argDist, double arg1, double arg2, double arg3){
	_eng.seed(argSeed);
	_distrType = argDist;
	_distrParam[0] = arg1;
	_distrParam[1] = arg2;
	_distrParam[2] = arg3;
	return 0;
}

double CRndVarGen::GetRV() {
	switch (_distrType) {
		//enum TYPE_RANDOM_DISTR {U01, NORMAL, EXP, GAMMA, BETA, UNIFREAL, EMPCDF, CONSTANT, DISCRETE, UNIFINT};
	case EXP: {
				  exponential_distribution<double> distrExp(_distrParam[0]);
				  return distrExp(_eng);
				  break;
	}
	case NORMAL: {
					 normal_distribution<double> distrNormal(_distrParam[0], _distrParam[1]);
					 return distrNormal(_eng);
					 break;
	}
	case GAMMA: {
					gamma_distribution<double> distrGamma(_distrParam[0], _distrParam[1]);
					return distrGamma(_eng);
					break;
	}
	case BETA: {
				   gamma_distribution<double> distrGamma1(_distrParam[0], 1);
				   gamma_distribution<double> distrGamma2(_distrParam[1], 1);
				   double x = distrGamma1(_eng);
				   double y = distrGamma2(_eng);
				   return x / (x + y);
				   break;
	}case UNIFINT: {
				   uniform_int_distribution<int> distrInt((int)_distrParam[0], (int)_distrParam[1]);
				   return distrInt(_eng);
				   break;
	}case UNIFREAL:{
		uniform_real_distribution<double> distrUnifReal(_distrParam[0], _distrParam[1]);
		return distrUnifReal(_eng);
		break;
	}
	case CONSTANT: {
		return _distrParam[0];
		break;
	}default: {
		ExitWithMsg("Unknown type of random variable distribution...");
		break;
	}
	}
	return 0;
}
