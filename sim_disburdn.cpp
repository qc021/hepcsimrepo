#include "sim_disburdn.h"

int BurdenModelSim::ReadSmpParamForPSA(ifstream & inf, modelParamType & argModelParam)
{

	inf >> argModelParam._cohortData.baseState >> argModelParam._cohortData.baseGender >> argModelParam._cohortData.baseAge
		>> argModelParam._qolData.q_F0 >> argModelParam._qolData.q_F1 >> argModelParam._qolData.q_F2 >> argModelParam._qolData.q_F3 >> argModelParam._qolData.q_CoCirr
		>> argModelParam._qolData.q_DeCirr >> argModelParam._qolData.q_HCC >> argModelParam._qolData.q_LivTr >> argModelParam._qolData.q_PostLivTr
		>> argModelParam._qolData.q_SVR >> argModelParam._qolData.q_Dec_Anemia >> argModelParam._qolData.q_TX_oSOC >> argModelParam._qolData.q_TX_DAA
		>> argModelParam._costData.c_F0 >> argModelParam._costData.c_F1 >> argModelParam._costData.c_F2 >> argModelParam._costData.c_F3 >> argModelParam._costData.c_CoCirr
		>> argModelParam._costData.c_DeCirr >> argModelParam._costData.c_DeCirr1yrPlus >> argModelParam._costData.c_HCC >> argModelParam._costData.c_LivTr >> argModelParam._costData.c_PostLivTr
		>> argModelParam._transData.pr_F0_F1 >> argModelParam._transData.pr_F1_F2 >> argModelParam._transData.pr_F2_F3 >> argModelParam._transData.pr_F3_CoCirr
		>> argModelParam._transData.pr_CoCirr_DeCirr >> argModelParam._transData.pr_CoCirr_HCC >> argModelParam._transData.pr_DeCirr_HCC
		>> argModelParam._transData.pr_DeCirr_LivTr >> argModelParam._transData.pr_DeCirr_DeathLiv >> argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv
		>> argModelParam._transData.pr_HCC_LivTr >> argModelParam._transData.pr_HCC_DeathLiv >> argModelParam._transData.pr_LivTr_DeathLiv
		>> argModelParam._transData.pr_LivTr1yrPlus_DeathLiv >> argModelParam._transData.pr_SVR_CoCirr_DeCirr >> argModelParam._transData.pr_SVR_CoCirr_HCC
		>> argModelParam._transData.pr_SVR_Delta_oSOC >> argModelParam._transData.pr_SVR_Delta_DAA;
	return 0;
}

int BurdenModelSim::PSA_initializeSampler(const vector<string> & vecVarName, const vector<string> & vecDistrName, const vector<double> & vecP1, const vector<double> & vecP2)
{
	seed_type s = _psa_seed();
	for (int k = 0; k < vecVarName.size(); k++) {
		CRndVarGen rvGen;
		seed_type theSeed = PSA_SAME_SEED_FOR_ALL_PARAM ? s : (_psa_seed());
		TYPE_RANDOM_DISTR  theType;
		if ("Beta" == vecDistrName[k])
			theType = BETA;
		else if ("Gamma" == vecDistrName[k])
			theType = GAMMA;
		else if ("Uniform" == vecDistrName[k])
			theType = UNIFREAL;
		else
			ExitWithMsg("[Error]: PSA_initializeSampler: unknown distr type in the input file: " + vecDistrName[k]);
		rvGen.Initialize(theSeed, theType, vecP1[k], vecP2[k]);
		_psa_sampler[vecVarName[k]] = rvGen;

	}
	return 0;


}

int BurdenModelSim::PSA_initializeSampler(const psaDistrType & argPSADistr)
{
	return PSA_initializeSampler(argPSADistr._distrPSA_varName, argPSADistr._distrPSA_distrName, argPSADistr._distrPSA_param1, argPSADistr._distrPSA_param2);
}

int BurdenModelSim::PSA_sampleModelParamValue(modelParamType & argModelParam)
{
	for (map<string, CRndVarGen>::const_iterator it = _psa_sampler.begin(); it != _psa_sampler.end(); it++) {
		string varName = it->first;
		if ("pF0_F1" == varName) {
			argModelParam._transData.pr_F0_F1 = _psa_sampler[varName].GetRV();
		}
		else if ("pF1_F2_SA" == varName) {
			argModelParam._transData.pr_F1_F2 = _psa_sampler[varName].GetRV();
		}
		else if ("pF2_F3_SA" == varName) {
			argModelParam._transData.pr_F2_F3 = _psa_sampler[varName].GetRV();
		}
		else if ("pF3_F4_SA" == varName) {
			argModelParam._transData.pr_F3_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("pF4_DC_SA" == varName) {
			argModelParam._transData.pr_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("pF4_HCC_SA" == varName) {
			argModelParam._transData.pr_CoCirr_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("pDC_HCC_SA" == varName) {
			argModelParam._transData.pr_DeCirr_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("pDC_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_DeCirr_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_dc_cyc_1_SA" == varName) {
			argModelParam._transData.pr_DeCirr_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_dc_cyc_2_SA" == varName) {
			argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pHCC_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_HCC_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_hcc_cyc_SA" == varName) {
			argModelParam._transData.pr_HCC_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_LivTr_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_Post_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pr_SVR_CoCirr_DeCirr" == varName) {
			argModelParam._transData.pr_SVR_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("pr_SVR_CoCirr_HCC" == varName) {
			argModelParam._transData.pr_SVR_CoCirr_HCC = _psa_sampler[varName].GetRV();

		}
		else if ("qolType::q_Acute" == varName) {
			argModelParam._qolData.q_acute = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F0" == varName) {
			argModelParam._qolData.q_F0 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F1" == varName) {
			argModelParam._qolData.q_F1 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F2" == varName) {
			argModelParam._qolData.q_F2 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F3" == varName) {
			argModelParam._qolData.q_F3 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_CoCirr" == varName) {
			argModelParam._qolData.q_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_DeCirr" == varName) {
			argModelParam._qolData.q_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_HCC" == varName) {
			argModelParam._qolData.q_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_LivTr" == varName) {
			argModelParam._qolData.q_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_PostLivTr" == varName) {
			argModelParam._qolData.q_PostLivTr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_SVR" == varName) {
			argModelParam._qolData.q_SVR = _psa_sampler[varName].GetRV();

		}
		else if ("Anemia_multiplier" == varName) {
			argModelParam._qolData.q_Dec_Anemia = _psa_sampler[varName].GetRV();
		}
		else if ("Therapy-related_multiplier-oSOC" == varName) {
			argModelParam._qolData.q_TX_oSOC = _psa_sampler[varName].GetRV();
		}
		else if ("Therapy-related_multiplier-DAA" == varName) {
			argModelParam._qolData.q_TX_DAA = _psa_sampler[varName].GetRV();

		}
		else if ("costType::c_Acute" == varName) {
			argModelParam._costData.c_acute = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F0" == varName) {
			argModelParam._costData.c_F0 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F1" == varName) {
			argModelParam._costData.c_F1 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F2" == varName) {
			argModelParam._costData.c_F2 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F3" == varName) {
			argModelParam._costData.c_F3 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_CoCirr" == varName) {
			argModelParam._costData.c_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_DeCirr" == varName) {
			argModelParam._costData.c_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_DeCirr1yrPlus" == varName) {
			argModelParam._costData.c_DeCirr1yrPlus = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_HCC" == varName) {
			argModelParam._costData.c_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_LivTr" == varName) {
			argModelParam._costData.c_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_PostLivTr" == varName) {
			argModelParam._costData.c_PostLivTr = _psa_sampler[varName].GetRV();

		}
		else if ("pSVR_Delta_DAA" == varName) {
			argModelParam._transData.pr_SVR_Delta_DAA = _psa_sampler[varName].GetRV();
		}
		else if ("pSVR_Delta_oSOC" == varName) {
			argModelParam._transData.pr_SVR_Delta_oSOC = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F0" == varName) {
			argModelParam._dalyData._dw_f0 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F1" == varName) {
			argModelParam._dalyData._dw_f1 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F2" == varName) {
			argModelParam._dalyData._dw_f2 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F3" == varName) {
			argModelParam._dalyData._dw_f3 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_CoCirr" == varName) {
			argModelParam._dalyData._dw_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_DeCirr" == varName) {
			double v = _psa_sampler[varName].GetRV();
			argModelParam._dalyData._dw_DeCirr = v;
			argModelParam._dalyData._dw_DeCirr1yrPlus = v;
		}
		else if ("dalyType::dw_HCC" == varName) {
			argModelParam._dalyData._dw_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_LivTr" == varName) {
			argModelParam._dalyData._dw_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_PostLivTr" == varName) {
			argModelParam._dalyData._dw_LivTr1yrPlus = _psa_sampler[varName].GetRV();
		}

		else if ("costType::c_testing_preTx" == varName) {
			argModelParam._costData.c_testing_preTx = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_testing_postTx" == varName) {
			argModelParam._costData.c_testing_postTx = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_testing_RNA" == varName) {
			argModelParam._costData.c_testing_RNA = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_testing_genotype" == varName) {
			argModelParam._costData.c_testing_genotype = _psa_sampler[varName].GetRV();
		}
		else if ("pAcuteSelfClearance" == varName) {
			argModelParam._transData._pr_acute_spont_clearance = _psa_sampler[varName].GetRV();
		}
		else if ("pr_SVR" == varName) {
			argModelParam._transData._pr_svr_for_SA = _psa_sampler[varName].GetRV();
		}



		else {
			ExitWithMsg("[Error] SampleModelParamValue(): Unknown parameter " + varName);
		}
	}

	if (PSA_SAME_QOL_F0_F3) {
		argModelParam._qolData.q_F1 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F2 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F3 = argModelParam._qolData.q_F0;


		argModelParam._costData.c_F1 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F2 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F3 = argModelParam._costData.c_F0;

		argModelParam._dalyData._dw_f1 = argModelParam._dalyData._dw_f0;
		argModelParam._dalyData._dw_f2 = argModelParam._dalyData._dw_f0;
		argModelParam._dalyData._dw_f3 = argModelParam._dalyData._dw_f0;
	}

	return 0;
}


int BurdenModelSim::PSA_sampleModelParamValue_Univariate(modelParamType & argModelParam, string argVarName)
{
	map<string, CRndVarGen>::const_iterator it = _psa_sampler.find(argVarName);
	string varName = it->first;
	if ("pF0_F1" == varName) {
		argModelParam._transData.pr_F0_F1 = _psa_sampler[varName].GetRV();
	}
	else if ("pF1_F2_SA" == varName) {
		argModelParam._transData.pr_F1_F2 = _psa_sampler[varName].GetRV();
	}
	else if ("pF2_F3_SA" == varName) {
		argModelParam._transData.pr_F2_F3 = _psa_sampler[varName].GetRV();
	}
	else if ("pF3_F4_SA" == varName) {
		argModelParam._transData.pr_F3_CoCirr = _psa_sampler[varName].GetRV();
	}
	else if ("pF4_DC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("pF4_HCC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("pDC_HCC_SA" == varName) {
		argModelParam._transData.pr_DeCirr_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("pDC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_DeCirr_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_dc_cyc_1_SA" == varName) {
		argModelParam._transData.pr_DeCirr_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_dc_cyc_2_SA" == varName) {
		argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pHCC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_HCC_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_hcc_cyc_SA" == varName) {
		argModelParam._transData.pr_HCC_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_Post_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pr_SVR_CoCirr_DeCirr" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("pr_SVR_CoCirr_HCC" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_HCC = _psa_sampler[varName].GetRV();

	}
	else if ("qolType::q_F0" == varName) {
		argModelParam._qolData.q_F0 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_F1" == varName) {
		argModelParam._qolData.q_F1 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_F2" == varName) {
		argModelParam._qolData.q_F2 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_F3" == varName) {
		argModelParam._qolData.q_F3 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_CoCirr" == varName) {
		argModelParam._qolData.q_CoCirr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_DeCirr" == varName) {
		argModelParam._qolData.q_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_HCC" == varName) {
		argModelParam._qolData.q_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_LivTr" == varName) {
		argModelParam._qolData.q_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_PostLivTr" == varName) {
		argModelParam._qolData.q_PostLivTr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_SVR" == varName) {
		argModelParam._qolData.q_SVR = _psa_sampler[varName].GetRV();

	}
	else if ("Anemia_multiplier" == varName) {
		argModelParam._qolData.q_Dec_Anemia = _psa_sampler[varName].GetRV();
	}
	else if ("Therapy-related_multiplier-oSOC" == varName) {
		argModelParam._qolData.q_TX_oSOC = _psa_sampler[varName].GetRV();
	}
	else if ("Therapy-related_multiplier-DAA" == varName) {
		argModelParam._qolData.q_TX_DAA = _psa_sampler[varName].GetRV();

	}
	else if ("costType::c_F0" == varName) {
		argModelParam._costData.c_F0 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_F1" == varName) {
		argModelParam._costData.c_F1 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_F2" == varName) {
		argModelParam._costData.c_F2 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_F3" == varName) {
		argModelParam._costData.c_F3 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_CoCirr" == varName) {
		argModelParam._costData.c_CoCirr = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_DeCirr" == varName) {
		argModelParam._costData.c_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_DeCirr1yrPlus" == varName) {
		argModelParam._costData.c_DeCirr1yrPlus = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_HCC" == varName) {
		argModelParam._costData.c_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_LivTr" == varName) {
		argModelParam._costData.c_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_PostLivTr" == varName) {
		argModelParam._costData.c_PostLivTr = _psa_sampler[varName].GetRV();

	}
	else if ("pSVR_Delta_DAA" == varName) {
		argModelParam._transData.pr_SVR_Delta_DAA = _psa_sampler[varName].GetRV();
	}
	else if ("pSVR_Delta_oSOC" == varName) {
		argModelParam._transData.pr_SVR_Delta_oSOC = _psa_sampler[varName].GetRV();


	}
	else {
		ExitWithMsg("[Error] SampleModelParamValue(): Unknown parameter " + varName);
	}


	if (PSA_SAME_QOL_F0_F3) {
		argModelParam._qolData.q_F1 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F2 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F3 = argModelParam._qolData.q_F0;


		argModelParam._costData.c_F1 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F2 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F3 = argModelParam._costData.c_F0;

	}

	return 0;
}

int BurdenModelSim::PSA_sampleCohort(baseCohortType & argCohort)
{
	// sample gender
	if (_psa_sampler_gender.GetU01() < 0.64) {
		argCohort.baseGender = 'M';
	}
	else {
		argCohort.baseGender = 'F';
	}

	// sample initial state
	vector<stateType> label_state;
	vector<double> count_state;
	label_state.push_back(s_F0);		count_state.push_back(0.1);	//F0
	label_state.push_back(s_F1);		count_state.push_back(0.25);	//F1
	label_state.push_back(s_F2);		count_state.push_back(0.2);	//F2
	label_state.push_back(s_F3);		count_state.push_back(0.2);		//F3
	label_state.push_back(s_CoCirr);		count_state.push_back(0.25);	//F4
	argCohort.baseState = DiscreteDistrSampler(label_state, count_state, _psa_sampler_state);

	// sample initial age

	return 0;
}


BurdenModelSim::BurdenModelSim()
{
	_cum_num_pt_generated = 0;
	_randomSeed = 0;
	// initialize
	_tempOutf.open("output_temp.txt", ios::app);
	_tempOutf << fixed << setprecision(9);

	_psa_seed.seed(PSA_SEED);
	_psa_sampler_age.Initialize(_psa_seed(), UNIFREAL);
	_psa_sampler_state.Initialize(_psa_seed(), UNIFREAL);
	_psa_sampler_gender.Initialize(_psa_seed(), UNIFREAL);

	_rnd_pt.Initialize(SIM_SEED, UNIFREAL);
	_rnd_sim.Initialize(SIM_SEED + 1, UNIFREAL);
	_rnd_vector_shuffler.Initialize(SIM_SEED + 2, UNIFREAL);
	
	

	_num_aware_by_screening = 0;
	_num_aware_by_usual_care = 0;
	_num_screening = 0;


	vector<int> emptyVec;
	for (int k = 0; k <= END_YR - START_YR; k++) {
		_outputRowData.push_back(0);
		_outputRowData_previousCycle.push_back(0);
	}

	_start_year_redefined = START_YR;
}


int BurdenModelSim::Run(const modelParamType & argModelParam)
{

	clock_t startclock = clock();

	ResetCounters();
	PrintOutputHeader(argModelParam._disBurdnOutFile);

	_curModelCycle = 0;
	_curYear = START_YR;
	_listAllPts.clear();



	while (START_YR + _curModelCycle <= END_YR) {
		if (_curYear < _start_year_redefined) {
			PrintOutputRow(_curModelCycle);
			_curYear++;
			_curModelCycle++;
			continue;
		}

		if (_curYear == _start_year_redefined) {
			// ------------ [1] initialize the population ----------------
			GeneratePatient(argModelParam);
		}

		// ------------ [2] simulate each year the population ----------------

		int nDeath = 0;
		cout << "============ Year " << _curYear << " ==============" << endl;

		cout << "    Change of cohort groups"<<endl;
		ChangeOfCohort(argModelParam); // Prison release, all hospitalized change to NHANES next year

		cout << "    Add new incidences" << endl;
		GeneratePatient(argModelParam);
		if(argModelParam._disBurdnData._flag_multipleCohort){
			GeneratePatient_nonNHANES(argModelParam);

			if (argModelParam._disBurdnData._flag_include_immigrants_lpr) {
				GeneratePatient_immigrants_lpr(argModelParam);
			}
		}


		// update since cycle 1 (not cycle 0, according to the original code)
		if (_curModelCycle >= 1) {
			// --- update insurance information --------------
			cout << "    Update insurance status" << endl;
			for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
				UpdateInsurance_indvPt(*(*ptIter), argModelParam);
			}



			// --- [2.1] Screen patients for treatment
			cout << "    Screening" << endl;
			if (argModelParam._disBurdnData._flag_multipleCohort) {
				Screening_MultipleCohort(argModelParam);
			}
			else {
				Screening(argModelParam);
			}
			

			// --- [2.2] treatment
			cout << "    Treatment: determine eligibility;" << endl;	
			if (_curYear >= TX_START_YR) {

				//Treatment_DetermineElig_OrignalVersion(argModelParam);
				if (argModelParam._disBurdnData._flag_extended_use_PR_PI) {
					Treatment_DetermineElig_w_ExtendedUsePRandPI(argModelParam);
				}
				else {
					// identify treatment eligible patients
					Treatment_DetermineElig(argModelParam);
				}
			}

			UpdateCounter_TxEligibility();

			cout << "    Treatment: prioritize treatment;" << endl;
			
			if (argModelParam._disBurdnData._flag_multipleCohort) {
				Treatment_Prioritize_MultipleCohort(argModelParam);
			}
			else {
				Treatment_Prioritize(argModelParam);
			}
			cout << "    Treatment: treat patients;" << endl;
			Treatment_Treat(argModelParam);


	

			// --- [2.3] run natural history/progression
			cout << "    Run natural history" << endl;
			for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ) {
				if ((*ptIter)->_state == state_NA) { ExitWithMsg("Error @ Natural history: patient state = state_NA!"); }

				// update treatment failure counter first!
				UpdateCounter_Treamtent(*(*ptIter), argModelParam);				
				RunNaturalHistory_indvPt(*(*ptIter), argModelParam);				
				ptIter++;
			}

		}// for if (cycleNum>=1)

		 // --- remove dead patients and update counters for alive patients --------------
		cout << "    Remove dead patients and update counters for alive patients;" << endl;
		for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ) {
			if ((*ptIter)->_state == state_NA) { ExitWithMsg("Error @ Natural history: patient state = state_NA!"); }


			// remove dead patient
			if ((*ptIter)->_state == state_Death || (*ptIter)->_state == state_DeathLR) {				
				UpdateCounter_DeathEvent(*(*ptIter), argModelParam);

				FreeMem(*ptIter); // free the memory of the patient object
				ptIter = _listAllPts.erase(ptIter);				
				
				nDeath++;
			}
			else {
				
				// update counters for all alive patients					
				UpdateCounter_HealthState(*(*ptIter), argModelParam);				
				UpdateCounter_Insurance(*(*ptIter));
				ptIter++;
			}
		}// end for each patient


		PrintOutputRow(_curModelCycle);
		

#if defined(__PAKISTAN_SETTING__)
		if (_curYear == YEAR_RECORD_AGE_DISTR) {OutputAgeDistr(_simCounter._counter_ageDistr_xsectional_given_year, 
			"output_project_pakistan/output_age_distr_prev_"+basicToStr(YEAR_RECORD_AGE_DISTR)+".txt"); }
		if (_curYear == _simCounter._record_ageDistr_year_lrd) {
			OutputAgeDistr(_simCounter._counter_LRD_ageDistr, 
			"output_project_pakistan/output_age_distr_lrd_"
			+basicToStr(_simCounter._year_start_lrd_ageDistr)+"_"
			+basicToStr(_simCounter._record_ageDistr_year_lrd)+".txt"); 
		}
#endif


		if (_curYear == argModelParam._disBurdnData._record_ageDistr_year) {
			OutputAgeDistr(_simCounter._counter_ageDistr,
				argModelParam._disBurdnData._record_ageDistr_outputFolder + "/output_age_distr_"
				+ basicToStr(argModelParam._disBurdnData._record_ageDistr_year) + ".txt");
		}
		cout << fixed << setprecision(3)
			<< "    SUMMARY of year " << _curYear << ":  total # pts = " << _listAllPts.size() << "; \t"
			<< "+ " << _simCounter._counter_new_hcv_incidence[_curModelCycle] << "\t"
			<< "- " << nDeath << "\t"
			<< setprecision(0) << "[" << Time(startclock) << " sec]"
			<< endl << endl;

		_curModelCycle++;
		_curYear++;
	}// end for year

	// release memory for all patients
	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		FreeMem(*ptIter);
	}

	return 0;

}



int BurdenModelSim::GeneratePatient(const modelParamType & argModelParam)
{
	// ---- list of patient attributes -----
	//- gender;
	//- initial age, current Age;
	//- genotype;

	//- bool _flagReceivable;
	//- bool _flagTrR;
	//- bool _flagBirthCohort_45_65_yr;

	//- typeInsurance _insurStatus;
	//- bool _flagAware;

	//- Tx Naive  vs. tx experienced
	//- typeStateDisBurdnModel _state;
	//- typeTxResponseStatus _curTxResps;

	// ----------------------------------------
	bool isForInitialization = _listAllPts.empty();

	// modified on 4/23/2017
	//long numNewPt = isForInitialization ? argModelParam._disBurdnData._num_initial_population : argModelParam._disBurdnData._table_incidence.find(_curYear)->second;
	long numNewPt;
	if(isForInitialization){
		numNewPt = argModelParam._disBurdnData._num_initial_population;
	}else{
		if(argModelParam._disBurdnData._flag_inci_proportional_to_prevalence
			|| _curYear<argModelParam._disBurdnData._year_start_constant_new_incidence){
			// For Pakistan analysis, and State-X analysis:
			// calculate incidence based on the prevalence			
			int curHCVPrev = _curModelCycle == 0? argModelParam._disBurdnData._num_initial_population : _simCounter._counterHCV[_curModelCycle-1];
			numNewPt = (argModelParam._disBurdnData._table_incidence_proportion.find(_curYear)->second) * curHCVPrev;
		}else{
			// use preset table for the incidence
			if (argModelParam._disBurdnData._flag_multipleCohort) {
				// ------------------- updated @ 7/20/2017 ------------------------------------------------------------------------------------------
				// Notes: allocate the incidences to the default NHANES population 
				//        based on the % of prevalence among other populations
				//        in the previous year. 
				double pNHANES = 1.0;
				if (_curYear == _start_year_redefined) {
					int totalInitialPopulation = argModelParam._disBurdnData._num_initial_population;
					for (map<string, typeCohortDisBurdnModel>::const_iterator it = argModelParam._nonNHANESPplData._mapNonNHANESCohortName.begin();
						it != argModelParam._nonNHANESPplData._mapNonNHANESCohortName.end(); it++) {
						if (argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_year", it->second) == _start_year_redefined) {
							totalInitialPopulation += argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_num", it->second);
						}
					}
					pNHANES = argModelParam._disBurdnData._num_initial_population / (double) totalInitialPopulation;
				}
				else if(_curYear > _start_year_redefined) {
					pNHANES = _simCounter._counter_byCohort_HCV[cohort_defaultSingleCohort][_curYear -1] / (double)_simCounter._counterHCV[_curYear - 1];
				}
				else {
					ExitWithMsg("ERROR @ GenratePatients(): This function should not be excuted when _curYear < _start_year_redfined...");
				}

				numNewPt = pNHANES* argModelParam._disBurdnData._table_incidence.find(_curYear)->second;
			}
			else {
				// ------------ base case -----------------------------------------------------
				numNewPt = argModelParam._disBurdnData._table_incidence.find(_curYear)->second;
				// ------------ base case -----------------------------------------------------
			}
		}
	}
	_simCounter._counter_new_hcv_incidence[_curModelCycle] = numNewPt;
	

	for (int k = 0; k < numNewPt; k++) {
		if (isForInitialization && k % (numNewPt / 10) == 0) {
			cout << fixed << "- Year " << _curYear << ": Generating initial population ... " << (int)(double(k) / double(numNewPt) * 100) << " % completed..." << endl;
		}
		patientType * newPt = new patientType;
		newPt->_ptID = _cum_num_pt_generated;
		newPt->gender = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_male ? 'M' : 'F';
		newPt->_flagReceivable = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_receivable ? true : false;
		newPt->genotype = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_genotype, _rnd_pt);
		newPt->_flagTrR = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_TrR ? true : false;
		newPt->_state = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_chronic ? state_NA : state_curedAcute;
		//newPt->_curTxResps = newPt->_state == state_curedAcute ? trStatus_svr : trStatus_unknown;
		newPt->_curTxResps = trStatus_unknown;
		newPt->_yrMostRecentTx = -1;
		newPt->_yrMostRecentScr = -1;
		newPt->_yrDiagnosedByUsualCare = -1;	// added for cost burden model @ 6/5/2017
		newPt->_flagInterferonTol = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_IFNintol ? false : true; // added for cost burden model @ 6/5/2017

		// ------------ Sample initial age ------------------------------------------------------------------
		// for initialization ONLY
		vector<int> ageCategory = isForInitialization
			? argModelParam._disBurdnData._vec_age_category_forInitialization
			: argModelParam._disBurdnData._vec_age_category_forNewIncidence;

		vector<double> ageDistr = isForInitialization
			? argModelParam._disBurdnData._table_inci_age_distr.find(0)->second
			: argModelParam._disBurdnData._table_inci_age_distr.lower_bound(_curYear)->second;

		if (ageDistr.size() != ageCategory.size()) {
			ExitWithMsg("Error @ GeneratePatient(): ageCategory and ageDistr do not match...");
		}

		int idxAgeCat = DiscreteDistrSampler(ageDistr, _rnd_pt);



		if (idxAgeCat >= (int)ageDistr.size()) {
			ExitWithMsg("Error @ GeneratePatient(): index out of range for ageDistr");
		}

		if (idxAgeCat != ageCategory.size() - 1) {
			newPt->initialAge = FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION
				? _rnd_pt.GetUnifInt_Rounded(ageCategory[idxAgeCat], ageCategory[idxAgeCat + 1] - 1)
				: _rnd_pt.GetUnifInt(ageCategory[idxAgeCat], ageCategory[idxAgeCat + 1] - 1);
		}
		else {
				newPt->initialAge = min(100, ageCategory[idxAgeCat] + _rnd_pt.GetExp_Integer(1.0 / 7.0)); // temporary upper bound
				//newPt->initialAge = _rnd_pt.GetUnifInt(ageCategory[idxAgeCat], 85); // temporary upper bound
		}
		newPt->currentAge = newPt->initialAge;
		newPt->_yrBirth = _curYear - (int)newPt->currentAge;
		newPt->_flagBirthCohortForScr = (newPt->_yrBirth >= SCREEN_BIRTH1_YR && newPt->_yrBirth <= SCREEN_BIRTH2_YR);

		// ------------ Determine awareness and insurance types ---------------------------------------
		DetmInsrAndAwarenessStatus(*newPt, argModelParam, isForInitialization);

		// ------------ sample response status state (treatment naive/experienced/...) ----------------
		if (newPt->_flagAware && (newPt->_state != state_curedAcute)) {
			DetmRespsState(*newPt, argModelParam, isForInitialization);
		}


		// ------------ sample fibrosis states --------------------------------------------------------
		if (newPt->_state != state_curedAcute) {
			if (isForInitialization) {
				newPt->_state = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_fib, _rnd_pt);
			}
			else {
				newPt->_state = state_F0;
			}
		}

		newPt->_state_previous = newPt->_state;
		newPt->_stateAtTx;


		// ---------- treatment failure flags -----

		newPt->_tfaflag_eligible_for_tx = true;
		newPt->_tfaflag_num_fail_pegrbv = 0;
		newPt->_tfaflag_num_fail_pi = 0;
		newPt->_tfaflag_current_tx_type = txCat_NoTx;
		newPt->_tfaflag_previous_failed_tx_type = txCat_NoTx;
		newPt->_tfaflag_total_tx_times = 0;

		newPt->_tfaflag_ever_failed_nonns5a_before_2015 = false;
		newPt->_tfaflag_ever_failed_ns5a = false;
		newPt->_tfaflag_ever_failed_nonns5a_after_2015 = false;

		// ---------- 


		_listAllPts.push_back(newPt);
		_cum_num_pt_generated++;

	}

	return 0;
}


int BurdenModelSim::GeneratePatient_immigrants_lpr(const modelParamType & argModelParam)
{
	if (_curYear < argModelParam._nonNHANESPplData._immigrants_year_start) { return 0; }

	map<int, map<typeStateDisBurdnModel, map<int, int>>> numMap = argModelParam._nonNHANESPplData._immigrants_hcv_num_by_year_fib_genotype;
	typeStateDisBurdnModel listState[] = { state_F0,state_F1,state_F2,state_F3,state_CoCirr };
	int listGenotype[] = { 1,2,3,456 };
	for (int k = 0; k < sizeof(listState) / sizeof(typeStateDisBurdnModel); k++) {
		for (int j = 0; j < sizeof(listGenotype) / sizeof(int); j++) {
			


			bool isForInitialization = false;
			int numNewPt = numMap[_curYear][listState[k]][listGenotype[j]];

			// ---- reload the distribution for immigrants
			vector<int> ageCategory = argModelParam._nonNHANESPplData._immigrants_age_category;
			vector<double> ageDistr = argModelParam._nonNHANESPplData._immigrant_distr_age;
			for (int k = 0; k < numNewPt; k++) {
				// create new patient with "default" cohort type
				patientType * newPt = new patientType;

				// ---- assign the cohort type and birth cohort screening groups ------
				newPt->_cohort = cohort_defaultSingleCohort;
				newPt->_bcScr_group = bc_group_general;
				newPt->_tx_group = tx_group_general;

				newPt->_ptID = _cum_num_pt_generated;
				newPt->gender = _rnd_pt.GetU01() < argModelParam._nonNHANESPplData._immigrants_male_pct ? 'M' : 'F';  // replace with male ratio by cohort
				newPt->_flagReceivable = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_receivable ? true : false;
				newPt->genotype = listGenotype[j];
				newPt->_flagTrR = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_TrR ? true : false;
				newPt->_state = listState[k];
				newPt->_curTxResps = trStatus_unknown;
				newPt->_yrMostRecentTx = -1;
				newPt->_yrMostRecentScr = -1;
				newPt->_yrDiagnosedByUsualCare = -1; // added for cost burden model @ 6/5/2017
				newPt->_flagInterferonTol = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_IFNintol ? false : true; // added for cost burden model @ 6/5/2017

				if (ageDistr.size() != ageCategory.size()) {
					ExitWithMsg("Error @ GeneratePatient(): ageCategory and ageDistr do not match...");
				}

				int idxAgeCat = DiscreteDistrSampler(ageDistr, _rnd_pt);

				if (idxAgeCat >= (int)ageDistr.size()) {
					ExitWithMsg("Error @ GeneratePatient(): index out of range for ageDistr");
				}

				if (idxAgeCat != ageCategory.size() - 1) {
					newPt->initialAge = FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION
						? _rnd_pt.GetUnifInt_Rounded(ageCategory[idxAgeCat], ageCategory[idxAgeCat + 1] - 1)
						: _rnd_pt.GetUnifInt(ageCategory[idxAgeCat], ageCategory[idxAgeCat + 1] - 1);
				}
				else {
					newPt->initialAge = min(100, ageCategory[idxAgeCat] + _rnd_pt.GetExp_Integer(1.0 / 7.0)); // temporary upper bound
				}
				newPt->currentAge = newPt->initialAge;
				newPt->_yrBirth = _curYear - (int)newPt->currentAge;
				newPt->_flagBirthCohortForScr = (newPt->_yrBirth >= SCREEN_BIRTH1_YR && newPt->_yrBirth <= SCREEN_BIRTH2_YR);

				// ------------ Determine awareness and insurance types ---------------------------------------
				// Note @ 05/25/2017: some changes to prison population and Indian Reservation population
				DetmInsrAndAwarenessStatus(*newPt, argModelParam, isForInitialization);

				// ------------ sample response status state (treatment naive/experienced/...) ----------------
				// Note @ 05/25/2017: No change
				if (newPt->_flagAware && (newPt->_state != state_curedAcute)) {
					DetmRespsState(*newPt, argModelParam, isForInitialization);
				}

				// ------------ sample fibrosis states --------------------------------------------------------
				// ------------ [no need for immigrants] ------------------------------------------------------
				//// Note @ 05/25/2017: only changed the distribution of fibrosis states
				//if (newPt->_state != state_curedAcute) {
				//	if (isForInitialization) {
				//		newPt->_state = DiscreteDistrSampler_DistrTable(state_distr, _rnd_pt);
				//	}
				//	else {
				//		newPt->_state = state_F0;
				//	}
				//}

				newPt->_state_previous = newPt->_state;
				newPt->_stateAtTx;


				// ---------- treatment failure flags -----
				// Note @ 05/25/2017: No change
				newPt->_tfaflag_eligible_for_tx = true;
				newPt->_tfaflag_num_fail_pegrbv = 0;
				newPt->_tfaflag_num_fail_pi = 0;
				newPt->_tfaflag_current_tx_type = txCat_NoTx;
				newPt->_tfaflag_previous_failed_tx_type = txCat_NoTx;
				newPt->_tfaflag_total_tx_times = 0;

				newPt->_tfaflag_ever_failed_nonns5a_before_2015 = false;
				newPt->_tfaflag_ever_failed_ns5a = false;
				newPt->_tfaflag_ever_failed_nonns5a_after_2015 = false;

				// ---------- 

				_listAllPts.push_back(newPt);
				_cum_num_pt_generated++;
			}

			cout << "    Add ## immigrants with fib state " << listState[k] << " genotype " <<listGenotype[j]<<" = "<< numNewPt << endl;

		}
	}
	return 0;
}



int BurdenModelSim::GeneratePatient_nonNHANES(const modelParamType & argModelParam)
{
	
	for (map<string, typeCohortDisBurdnModel>::const_iterator it = argModelParam._nonNHANESPplData._mapNonNHANESCohortName.begin();
		it != argModelParam._nonNHANESPplData._mapNonNHANESCohortName.end(); it++) {
		
		bool isForInitialization = false;

		// ---- calculate number of new patients to be added
		int numNewPt = 0;
		int numNewPt_initialPopulation = 0;
		int yearOfIncludingInitialPopulation = (int)argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_year", it->second);

		
		if (_curYear == yearOfIncludingInitialPopulation) {
			// ***** ADDING initial population *******
			numNewPt_initialPopulation = (int)argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_num", it->second);
			numNewPt = numNewPt_initialPopulation;
		}
		else {
		}


		if (_curYear >= yearOfIncludingInitialPopulation) {
			// ***** ADDING NEW incidences *******
			// ***** ONLY for prison and indianreservation populations
			// NOTES: modified by Qiushi @ 7/20/2017
			// Calculate the new incidences for [Prison] and [Indian reservation] population
			// Distribute some incidences from the total estimates (from CDC) based on the previous-year prevalence. 
			if (it->second == cohort_incarcerated || it->second == cohort_indianReservation) {

				double prop_new_CDC_inci = _simCounter._counter_byCohort_HCV[it->second][_curYear - 1] / (double)_simCounter._counterHCV[_curYear - 1];
				long new_CDC_inci = argModelParam._disBurdnData._table_incidence.find(_curYear)->second;
				numNewPt += (int)prop_new_CDC_inci * new_CDC_inci;


			}

		} else {
			ExitWithMsg("ERROR @ GenratePatients_nonHNAMES(): Error in conditions of _curYear vs. yearOfIncludingInitialPopulation");
		}


		if (numNewPt == 0) {
			continue;
		}

		// ---- we include incidences into the counter ------
		_simCounter._counter_new_hcv_incidence[_curModelCycle] += numNewPt;
		



		// ---- reload the distribution for specific nonNHANES cohort
		map<int, double> geno_distr = argModelParam._nonNHANESPplData._distr_genotype_by_cohort.find(it->second)->second;
		map<typeStateDisBurdnModel, double> state_distr = argModelParam._nonNHANESPplData._distr_fib_by_cohort.find(it->second)->second;

		vector<int> age_category_initial = argModelParam._nonNHANESPplData._distr_age_initial_category_by_cohort.find(it->second)->second;
		vector<double> age_distr_initial = argModelParam._nonNHANESPplData._distr_age_initial_by_cohort.find(it->second)->second;

		// ---- reload age distribution for annual new incidences only for Indian Reservation & Prison population
		vector<int> age_category_new_inci;
		vector<double> age_distr_new_inci;
		if (it->second == cohort_indianReservation || it->second == cohort_incarcerated) {			
			age_category_new_inci = argModelParam._nonNHANESPplData._distr_age_new_inci_category_by_cohort.find(it->second)->second;
			age_distr_new_inci = argModelParam._nonNHANESPplData._distr_age_new_inci_by_cohort.find(it->second)->second;
		}

		
		for (int k = 0; k < numNewPt; k++) {
			// initial population: for the first *numNewPt_initialpopulation* people
			// new incidence: for the remaining (numNewPt - numNewPt_initialPopulation), which may be 0.
			isForInitialization = (k < numNewPt_initialPopulation)? true:false;


			// create new patient with "default" cohort type
			patientType * newPt = new patientType;

			// ---- assign the cohort type and birth cohort screening groups ------
			newPt->_cohort = it->second;
			if (it->second == cohort_homeless || it->second == cohort_hospitalized || it->second == cohort_nursingHome || it->second == cohort_defaultSingleCohort) {
				newPt->_bcScr_group = bc_group_general;
				newPt->_tx_group = tx_group_general;
			}
			else if (it->second == cohort_indianReservation) {
				newPt->_bcScr_group = bc_group_indianRsv;
				newPt->_tx_group = tx_group_indianRsv;
			}
			else if (it->second == cohort_incarcerated) {
				newPt->_bcScr_group = bc_group_noBirthCohortScr;
				newPt->_tx_group = tx_group_incarcerated;
			}
			else {
				ExitWithMsg("Error @ GeneratePatient_nonNHANES(): unknown cohort type = "+ it->first + " in tx/scr group assignment");
			}

			newPt->_ptID = _cum_num_pt_generated;			
			newPt->gender = _rnd_pt.GetU01() < argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_male",it->second) ? 'M' : 'F';  // replace with male ratio by cohort
			newPt->_flagReceivable = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_receivable ? true : false;
			newPt->genotype = DiscreteDistrSampler_DistrTable(geno_distr, _rnd_pt); // replace with genotype distribution by cohort
			newPt->_flagTrR = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_TrR ? true : false;
			newPt->_state = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_chronic ? state_NA : state_curedAcute;
			newPt->_curTxResps = trStatus_unknown;
			newPt->_yrMostRecentTx = -1;
			newPt->_yrMostRecentScr = -1;
			newPt->_yrDiagnosedByUsualCare = -1; // added for cost burden model @ 6/5/2017
			newPt->_flagInterferonTol = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_IFNintol ? false : true; // added for cost burden model @ 6/5/2017

			// ------------ Sample initial age ------------------------------------------------------------------
			// for initialization ONLY
			vector<int> ageCategory = isForInitialization
				? age_category_initial
				: age_category_new_inci;
				//? argModelParam._disBurdnData._vec_age_category_forInitialization
				//: argModelParam._disBurdnData._vec_age_category_forNewIncidence;

			vector<double> ageDistr = isForInitialization
				? age_distr_initial
				: age_distr_new_inci;

			if (ageDistr.size() != ageCategory.size()) {
				ExitWithMsg("Error @ GeneratePatient(): ageCategory and ageDistr do not match...");
			}

			int idxAgeCat = DiscreteDistrSampler(ageDistr, _rnd_pt);
			
			if (idxAgeCat >= (int)ageDistr.size()) {
				ExitWithMsg("Error @ GeneratePatient(): index out of range for ageDistr");
			}

			if (idxAgeCat != ageCategory.size() - 1) {
				newPt->initialAge = FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION
					? _rnd_pt.GetUnifInt_Rounded(ageCategory[idxAgeCat], ageCategory[idxAgeCat + 1] - 1)
					: _rnd_pt.GetUnifInt(ageCategory[idxAgeCat], ageCategory[idxAgeCat + 1] - 1);
			}
			else {
				newPt->initialAge = min(100, ageCategory[idxAgeCat] + _rnd_pt.GetExp_Integer(1.0 / 7.0)); // temporary upper bound
			}
			newPt->currentAge = newPt->initialAge;
			newPt->_yrBirth = _curYear - (int)newPt->currentAge;
			newPt->_flagBirthCohortForScr = (newPt->_yrBirth >= SCREEN_BIRTH1_YR && newPt->_yrBirth <= SCREEN_BIRTH2_YR);


			// ------------ Determine awareness and insurance types ---------------------------------------
			// Note @ 05/25/2017: some changes to prison population and Indian Reservation population
			DetmInsrAndAwarenessStatus(*newPt, argModelParam, isForInitialization);

			// ------------ sample response status state (treatment naive/experienced/...) ----------------
			// Note @ 05/25/2017: No change
			if (newPt->_flagAware && (newPt->_state != state_curedAcute)) {
				DetmRespsState(*newPt, argModelParam, isForInitialization);
			}


			// ------------ sample fibrosis states --------------------------------------------------------
			// Note @ 05/25/2017: only changed the distribution of fibrosis states
			if (newPt->_state != state_curedAcute) {
				if (isForInitialization) {
					newPt->_state = DiscreteDistrSampler_DistrTable(state_distr, _rnd_pt);
				}
				else {
					newPt->_state = state_F0;
				}
			}

			newPt->_state_previous = newPt->_state;
			newPt->_stateAtTx;


			// ---------- treatment failure flags -----
			// Note @ 05/25/2017: No change
			newPt->_tfaflag_eligible_for_tx = true;
			newPt->_tfaflag_num_fail_pegrbv = 0;
			newPt->_tfaflag_num_fail_pi = 0;
			newPt->_tfaflag_current_tx_type = txCat_NoTx;
			newPt->_tfaflag_previous_failed_tx_type = txCat_NoTx;
			newPt->_tfaflag_total_tx_times = 0;

			newPt->_tfaflag_ever_failed_nonns5a_before_2015 = false;
			newPt->_tfaflag_ever_failed_ns5a = false;
			newPt->_tfaflag_ever_failed_nonns5a_after_2015 = false;

			// ---------- 

			_listAllPts.push_back(newPt);
			_cum_num_pt_generated++;
		}
		cout << "    Add ## Non-NHANES ## " << it->first << ": " << numNewPt << endl;

	}// end for each nonNHANES population




	return 0;
}



int BurdenModelSim::DetmInsrAndAwarenessStatus(patientType & argPt, const modelParamType & argModelParam, bool argIsForInit)
{
	// ***********************************************
	// ------- determine insurance status ------------
	// set the value of "argPt._insurStatus"
	// ***********************************************
	//  ------ For the default population ------------------------------------------------------------------------
	if (argPt._cohort == cohort_defaultSingleCohort) {

		// first determine medicare or not
		if (argPt.currentAge >= AGE_CUTOFF_MEDICARE) {
			argPt._insurStatus = insr_medicare;
		}
		// for younger patients
		else {
			// before 2014
			if (argIsForInit || _curYear < ACA_START_YR) {
				argPt._insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_insurance, _rnd_pt);
			}
			// since 2014 (ACA starts)
			else {
				argPt._insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_insurance_ACA_by_year.find(_curYear)->second, _rnd_pt);
			}
		}

	}
	//  ------ For the non-NHANES population [incacerated, homeless, nursing home, hospitalized, indian reservation] -------------------------------------
	else {
		argPt._insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._nonNHANESPplData._distr_insur_by_cohort.find(argPt._cohort)->second, _rnd_pt);
		// checkpoints:
		if (argPt._cohort == cohort_incarcerated) {
			if(argPt._insurStatus != insr_incarcerated)
			ExitWithMsg("ERROR @ DetmInsrAndAwarenessStatus(): incacerated population have wrong insurance types = " + basicToStr((int)argPt._insurStatus));
		}
		if (argPt._cohort == cohort_indianReservation) {
			if (argPt._insurStatus != insr_indianRsv)
				ExitWithMsg("ERROR @ DetmInsrAndAwarenessStatus(): indian reservation population have wrong insurance types = " + basicToStr((int)argPt._insurStatus));
		}
	}


	// ***********************************************
	// ------- determine awareness status ------------
	// set the value of "argPt._flagAware"
	// set the value of "argPt._flagTrInsur"
	// ***********************************************

	// -----------------------
	// for initial population
	// -----------------------
	if (argIsForInit) {

		// ========================= IF overrided by certain fixed input value =========================
		if (argModelParam._disBurdnData._flag_override_awareness_for_initial_population) {
			// modified by Qiushi @ 6/11/2017
			// for initial population, assign a particular awareness rate
			argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._awareness_for_initial_population_if_overrided);

			if (argPt._insurStatus == insr_medicare) {
				argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
			}
			else {
				argPt._flagTrInsur = true;
			}
		}
		// ========================= IF NOT overrided by certain fixed input value =========================
		// awareness is determined by insurance type
		else {
			double age = argPt.currentAge;
			// --------------------------------------- default cohort, indian reservation --------------------------------------------------------------------
			if (argPt._cohort == cohort_defaultSingleCohort || argPt._cohort == cohort_indianReservation) {
				// uninsured:
				if (argPt._insurStatus == insr_uninsured) {
					argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._table_aware_prob_initial_population_uninsured.lower_bound(age)->second); 
					argPt._flagTrInsur = false;
				}
				// insured:
				else {
					argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._table_aware_prob_initial_population_insured.lower_bound(age)->second);  
					if (argPt._insurStatus == insr_medicare) {
						argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
					}
					else {
						argPt._flagTrInsur = true;
					}
				}
			}
			// --------------------------------------- incacerated --------------------------------------------------------------------------------------------
			else if(argPt._cohort == cohort_incarcerated){
				argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_awareness", cohort_incarcerated));
				argPt._flagTrInsur = true;
			}
			// --------------------------------------- homeless, nusring home, hospitalized--------------------------------------------------------------------
			else {
				double pAware = _simCounter._counter_aware[_curYear - 1] /(double)_simCounter._counterHCV[_curYear - 1];
				if (argPt._insurStatus == insr_uninsured) {
					argPt._flagTrInsur = false;
				}else if (argPt._insurStatus == insr_medicare) {
					argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
				}
				else {
					argPt._flagTrInsur = true;
				}
			}// end cohort type
		}// end override part



	}
	// -----------------------
	// for new incidences
	// -----------------------
	else {
		double age = argPt.currentAge;
		// --------------------------------------- default cohort, indian reservation, and incacerated --------------------------------------------------------------------
		// now awareness% among new incidences = 7%, regardlesss of their insured/uninsured status and age
		// we still keep the old code structure, and only change the values of "_table_aware_prob_among_new_incidence_uninsured" to be 7% for all categories.
		if (argPt._cohort == cohort_defaultSingleCohort || argPt._cohort == cohort_indianReservation || argPt._cohort == cohort_incarcerated) {
			// uninsured:
			if (argPt._insurStatus == insr_uninsured) {
				argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured.lower_bound(age)->second);
				argPt._flagTrInsur = false;
			}
			// insured:
			else {
				argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._table_aware_prob_among_new_incidence_insured.lower_bound(age)->second);
				if (argPt._insurStatus == insr_medicare) {
					argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
				}
				else {
					argPt._flagTrInsur = true;
				}
			}
		}
		// --------------------------------------- homeless, hospitalized, nursing home population    --------------------------------------------------------------------
		else {
			ExitWithMsg("ERROR @ DetmInsrAndAwarenessStatus(): Should not generate new incidences for homeless/hospitalized/nursing home population: argPt._cohort = "+basicToStr((int)argPt._cohort));
		}
	} // end if initial or new-incidences



	return 0;


	// ******************************************************************************************************************************************
	// ********** below are the old version *****************************************************************************************************
	// ******************************************************************************************************************************************
	// ----- EXCEPTIONS: prison population -----------
	if (argPt._cohort == cohort_incarcerated) {
		argPt._insurStatus = insr_incarcerated;
		argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_awareness", cohort_incarcerated));
		argPt._flagTrInsur = true;
		return 0;
	}
	
	
	
	// -----------------------------------------------
	// ------- determine insurance status ------------
	// -----------------------------------------------
	// first determine medicare or not
	if (argPt.currentAge >= AGE_CUTOFF_MEDICARE) {
		argPt._insurStatus = insr_medicare;
	}
	// for younger patients
	else {
		// before 2014
		if (argIsForInit || _curYear < ACA_START_YR) {
			argPt._insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_insurance, _rnd_pt);
		}
		// since 2014 (ACA starts)
		else {
			argPt._insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_insurance_ACA_by_year.find(_curYear)->second, _rnd_pt);
		}

		// ----- OVERRIDE & EXCEPTIONS: Indian reservation population -----------
		if (argPt._cohort == cohort_indianReservation) {
			argPt._insurStatus = insr_indianRsv;
		}
		// ----------------------------------------------------------------------
	}


	double age = argPt.currentAge;

	// this recovers the WRONG code in the original version.
	if (RECOVER_SALVAGE_POSTER_VERSION) {
		if (argIsForInit) age = age + 10;
	}

	// -------------------------------------------------------------------------------
	// ------- determine awareness based on [age] and [insured/uninsured] ------------
	// -------------------------------------------------------------------------------

	// ******************************************** OVERRIDE: if the awareness in initial population is assigned with a specific value **********************
	if (argModelParam._disBurdnData._flag_override_awareness_for_initial_population && argIsForInit) {		
		// modified by Qiushi @ 6/11/2017
		// for initial population, assign a particular awareness rate
		argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._awareness_for_initial_population_if_overrided);
		
		if (argPt._insurStatus == insr_medicare) {
			argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
		}
		else {
			argPt._flagTrInsur = true;
		}
	}
	// ******************************************** If not for OVERRIDE ***************************************************************************************
	else {	

		// awareness is determined by insurance type

		// uninsured:
		if (argPt._insurStatus == insr_uninsured) {
			//argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured.lower_bound(age)->second); // commented on 7/16/2017
			
			double pAware = argIsForInit ?
				argModelParam._disBurdnData._table_aware_prob_initial_population_uninsured.lower_bound(age)->second :
				argModelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured.lower_bound(age)->second;

			argPt._flagAware = (_rnd_pt.GetU01() <= pAware);


			argPt._flagTrInsur = false;
		}
		// insured:
		else {
			//argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._disBurdnData._table_aware_prob_among_new_incidence_insured.lower_bound(age)->second);  // commented on 7/16/2017

			double pAware = argIsForInit ?
				argModelParam._disBurdnData._table_aware_prob_initial_population_insured.lower_bound(age)->second :
				argModelParam._disBurdnData._table_aware_prob_among_new_incidence_insured.lower_bound(age)->second;
			argPt._flagAware = (_rnd_pt.GetU01() <= pAware);

			if (argPt._insurStatus == insr_medicare) {
				argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
			}
			else {
				argPt._flagTrInsur = true;
			}
		}
	}


	//// ----- OVERRIDE & EXCEPTIONS: Indian reservation population -----------
	//if (argIsForInit) {
	//	double pCurAwarenessRate = ((double)_simCounter._counter_aware[_curModelCycle - 1]) / ((double)_simCounter._counterHCV[_curModelCycle - 1]);
	//	if (argPt._cohort == cohort_homeless) {
	//		argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_awareness", cohort_homeless));
	//	}
	//	else if (argPt._cohort == cohort_nursingHome) {
	//		//argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_awareness", cohort_nursingHome));
	//		argPt._flagAware = (_rnd_pt.GetU01() <= pCurAwarenessRate);
	//	}
	//	else if (argPt._cohort == cohort_hospitalized) {
	//		//argPt._flagAware = (_rnd_pt.GetU01() <= argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_awareness", cohort_hospitalized));
	//		argPt._flagAware = (_rnd_pt.GetU01() <= pCurAwarenessRate);
	//	}
	//	else {
	//		//do nothing
	//	}
	//}
	//
	//// ----------------------------------------------------------------------

	return 0;

}

int BurdenModelSim::DetmRespsState(patientType & argPt, const modelParamType & argModelParam, bool argIsForInit)
{
	// for initialization
	if (argIsForInit) {
		if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_naive_initialization) {
			argPt._curTxResps = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_respStatus_distr_naive, _rnd_pt);			
		}
		else {
			if (argPt.genotype == 1) {
				argPt._curTxResps = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_respStatus_distr_expr_G1, _rnd_pt);
			}
			else {
				argPt._curTxResps = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_respStatus_distr_expr_G234, _rnd_pt);
			}
			argPt._cycleOnDrug_PEGRBV = 1;

		}
	}
	// for new incidence added to the model
	else {
		argPt._curTxResps = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_respStatus_distr_naive, _rnd_pt);

	}	

	return 0;
}



int BurdenModelSim::Screening_MultipleCohort(const modelParamType & argModelParam)
{
	// ============ 3 types of screening policies ======================================================================================
	// universal screening rate <Incarcerated population>
	// usual care, depending on fibrosis, insured/uninsured, and age <default, homeless, hospitalized, nursing home, Indian Reservation>	
	// birth cohort screening: <default, hospitalized, nursing home> <Indian Reservation> 
	// =================================================================================================================================
	_num_aware_by_screening = 0;
	_num_aware_by_usual_care = 0;
	_num_screening = 0;


	// ----------- incacerated population ------------
	Screening_Universal_Rate(argModelParam, cohort_incarcerated);


	// ----------- all other population ------------
	// --- Change: skipped usual care for incacerated population
	Screening_UsualCare(argModelParam);


	// ----------- birth cohort screening ---------------------------------------------------------------------------------------------------------
	if (_curYear == SCREEN_BIRTHCOHORT_START_YR) {
		int spanYr = SCREEN_BIRTHCOHORT_END_YR - SCREEN_BIRTHCOHORT_START_YR + 1;
		_scrCap_birthCohort_byGroup[bc_group_general] 		= 
			(_simCounter._counter_byCohort_HCV[cohort_homeless][_curModelCycle - 1] - _simCounter._counter_byCohort_awareness[cohort_homeless][_curModelCycle - 1]
			+ _simCounter._counter_byCohort_HCV[cohort_hospitalized][_curModelCycle - 1] - _simCounter._counter_byCohort_awareness[cohort_hospitalized][_curModelCycle - 1]
			+ _simCounter._counter_byCohort_HCV[cohort_nursingHome][_curModelCycle - 1] - _simCounter._counter_byCohort_awareness[cohort_nursingHome][_curModelCycle - 1]
			+ _simCounter._counter_byCohort_HCV[cohort_defaultSingleCohort][_curModelCycle - 1] - _simCounter._counter_byCohort_awareness[cohort_defaultSingleCohort][_curModelCycle - 1]) / spanYr;
		_scrCap_birthCohort_byGroup[bc_group_general] = (int)(_scrCap_birthCohort_byGroup[bc_group_general] * argModelParam._nonNHANESPplData._bcscr_coverage);

		_scrCap_birthCohort_byGroup[bc_group_indianRsv] = (_simCounter._counter_byCohort_HCV[cohort_indianReservation][_curModelCycle - 1]- _simCounter._counter_byCohort_awareness[cohort_indianReservation][_curModelCycle - 1]) / spanYr;
		_scrCap_birthCohort_byGroup[bc_group_indianRsv] = (int)(_scrCap_birthCohort_byGroup[bc_group_indianRsv] * argModelParam._nonNHANESPplData._bcscr_coverage);

	}
	// -- in homeless, hospitalized, nusring home, and default populations --------
	Screening_BirthCohort_byGroup(argModelParam, bc_group_general, _scrCap_birthCohort_byGroup[bc_group_general]);

	// -- birth cohort screening in Indian reservation populations -----------------
	Screening_BirthCohort_byGroup(argModelParam, bc_group_indianRsv, _scrCap_birthCohort_byGroup[bc_group_indianRsv]);



	cout << "    Year " << _curYear << " - Screening - #becomingAware=" << _num_aware_by_usual_care << ", #screened=" << _num_aware_by_screening << endl;

	return 0;
}


int BurdenModelSim::Screening(const modelParamType & argModelParam)
{
	_num_aware_by_screening = 0;
	_num_aware_by_usual_care = 0;
	_num_screening = 0;

	Screening_UsualCare(argModelParam);

	switch (argModelParam._disBurdnData._option_screenScenario)
	{
	case screen_no_extra:
		break;
	case screen_birthCohort:
		Screening_BirthCohort(argModelParam);
		break;
	case screen_universal_rate:
		Screening_Universal_Rate(argModelParam);
		break;
	case screen_universal_capacity:
		Screening_Universal_Capacity(argModelParam);
		break;
	default:
		ExitWithMsg("Error @ BurdenModelSim::Screening: Unknown value of _option_screenScenario = " + basicToStr(argModelParam._disBurdnData._option_screenScenario));
		break;
	}
	cout << "    Year " << _curYear << " - Screening - #becomingAware=" << _num_aware_by_usual_care << ", #screened=" << _num_aware_by_screening << endl;

	return 0;
}

int BurdenModelSim::Screening_UsualCare(const modelParamType & argModelParam)
{
	//
	// transition to awareness with probability by fibrosis, insurance, and age category
	//

	cout << "    -- Screening: usual care" << endl;

	//map<typeStateDisBurdnModel, map<bool, map<double, double>>> mapProbAware = argModelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age;

	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		if (!(*ptIter)->_flagAware && (*ptIter)->_state <= state_CoCirr) {
			if ((*ptIter)->_flagCured) { ExitWithMsg("Screening_UsualCare(): _flagCured is true"); }
			// if the patient is not aware of infection
			// he becomes aware with certain probability
			double pAware = 0;
			if ((*ptIter)->_cohort == cohort_defaultSingleCohort ||
				(*ptIter)->_cohort == cohort_nursingHome ||
				(*ptIter)->_cohort == cohort_hospitalized ||
				(*ptIter)->_cohort == cohort_homeless ||
				(*ptIter)->_cohort == cohort_indianReservation) {
				// argModelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age[ptIter->_state][ptIter->_insurStatus != insr_uninsured].lower_bound(ptIter->currentAge)->second;

				// modified @ 4/29/2017 by Qiushi
				map<typeStateDisBurdnModel, map<bool, map<double, double>>>::const_iterator  it1 = argModelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.find((*ptIter)->_state);
				if (it1 == argModelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.end()) {
					ExitWithMsg("Error @ Screening_UsualCare(): Can't find state " + basicToStr((*ptIter)->_state) + " in the table _disBurdnData._table_prob_becoming_aware_by_fib_insur_age");
				}
				map<bool, map<double, double>>::const_iterator it2 = (it1->second).find((*ptIter)->_insurStatus != insr_uninsured);
				if (it2 == it1->second.end()) {
					ExitWithMsg("Error @ Screening_UsualCare(): Can't identify the insurance status" + basicToStr((*ptIter)->_insurStatus) + " in the table _disBurdnData._table_prob_becoming_aware_by_fib_insur_age");
				}
				pAware = (it2->second).lower_bound((*ptIter)->currentAge)->second;
			}
			else if ((*ptIter)->_cohort == cohort_incarcerated) {
				continue; // skipped incacerated population
			}
			else {
				ExitWithMsg("ERROR @ Screening_UsualCare(): Unknown cohort type = " + basicToStr((int)((*ptIter)->_cohort)));
			}// end if cohort types

			if (_rnd_pt.GetU01() < pAware) {
				(*ptIter)->_flagAware = true;
				(*ptIter)->_flagAwareByUsualCare = true;
				(*ptIter)->_curTxResps = trStatus_naive;
				(*ptIter)->_yrDiagnosedByUsualCare = _curYear;
				_num_aware_by_usual_care++;

			}
		}
	}
	return 0;
}


int BurdenModelSim::Screening_BirthCohort_byGroup(const modelParamType & argModelParam, typeBirthCohortGroup argBCGroup, int argScrCap)
{
	int num_scr_birthCohort = 0;
	if (_curYear < SCREEN_BIRTHCOHORT_START_YR || _curYear > SCREEN_BIRTHCOHORT_END_YR) {
		return 0;
	}

	cout << "    -- Screening: Birth cohort" << endl;
	// only for 2013-2018
	// find the target population
	vector<patientType*> vecPt;

	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		if ((*ptIter)->_bcScr_group != argBCGroup) continue;

		if ((!(*ptIter)->_flagAware) && (*ptIter)->_flagBirthCohortForScr) {
			if ((*ptIter)->_flagCured) {
				ExitWithMsg("Screening_BirthCohort_byGroup(): _flagCured is true");
			}
			vecPt.push_back(*ptIter);
		}
	}
	// ----------------- screening capacity ---------------------------------
	// random sample the patients for birth-cohort screening	
	

	if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
		_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
	}
	else {
		std::random_shuffle(vecPt.begin(), vecPt.end());
	}

	for (int k = 0; (k < (int)vecPt.size() - 1) && (num_scr_birthCohort <= argScrCap); k++) {
		num_scr_birthCohort++; // moved from the inner loop, it represents # of screening, not # of diagnosed patients.

		double rnd = _rnd_pt.GetU01(); 
		if ((vecPt[k]->_insurStatus != insr_uninsured && rnd < argModelParam._disBurdnData._prob_birthCohortScr_insured) // insured
			|| (vecPt[k]->_insurStatus == insr_uninsured && rnd < argModelParam._disBurdnData._prob_birthCohortScr_uninsured)) {//uninsured			
			if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_screen_accept) { // percentage of patients who accepted screening and received correct results
				vecPt[k]->_flagAware = true;
				vecPt[k]->_flagAwareByScreening = true;
				vecPt[k]->_curTxResps = trStatus_naive;
				vecPt[k]->_yrMostRecentScr = _curYear;
				_num_aware_by_screening++;
				_num_screening++;
				//num_scr_birthCohort++;
			}
		}
	} // for vecPt



	return 0;
}

int BurdenModelSim::Screening_BirthCohort(const modelParamType & argModelParam)
{
	cout<<"    -- Screening: Birth cohort"<<endl;
	_num_screening = 0;

	if (_curYear >= SCREEN_BIRTHCOHORT_START_YR && _curYear <= SCREEN_BIRTHCOHORT_END_YR) { // only for 2013-2018

		// find the target population
		_vecPt_BirthCohortScr.clear();

		for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
			if ((!(*ptIter)->_flagAware) && (*ptIter)->_flagBirthCohortForScr) {
				if ((*ptIter)->_flagCured) {
					ExitWithMsg("Screening_BirthCohort(): _flagCured is true");
				}
				_vecPt_BirthCohortScr.push_back(*ptIter);
			}
		}

		if (FIXED_CAP_BIRTH_COHORT_SCREENING) {
			// ----------------- screening capacity ---------------------------------
			// random sample the patients for birth-cohort screening	
			vector<patientType*> vecPt = _vecPt_BirthCohortScr;
			
			if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
				_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
			}
			else {
				std::random_shuffle(vecPt.begin(), vecPt.end());
			}
						
			//std::shuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
			

			for (int k = 0; (k < (int)vecPt.size() - 1) && (_num_screening <= NUM_SCREEN_BIRTH_COHORT_EACH_YEAR); k++) {

				double rnd = _rnd_pt.GetU01();
				if ((vecPt[k]->_insurStatus != insr_uninsured && rnd < argModelParam._disBurdnData._prob_birthCohortScr_insured) // insured
					|| (vecPt[k]->_insurStatus == insr_uninsured && rnd < argModelParam._disBurdnData._prob_birthCohortScr_uninsured)) {//uninsured
					if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_screen_accept) { // percentage of patients who accepted screening and received correct results
						vecPt[k]->_flagAware = true;
						vecPt[k]->_flagAwareByScreening = true;
						vecPt[k]->_curTxResps = trStatus_naive;
						vecPt[k]->_yrMostRecentScr = _curYear;
						_num_aware_by_screening++;
						_num_screening++;
					}
				}
			} // for vecPt

		}
		else {
			// ----------------- screening rate ---------------------------------
			for (int k = 0; k < (int)_vecPt_BirthCohortScr.size(); k++) {
				double rnd = _rnd_pt.GetU01();

				if ((_vecPt_BirthCohortScr[k]->_insurStatus != insr_uninsured && rnd < argModelParam._disBurdnData._prob_birthCohortScr_insured / (double)(SCREEN_BIRTHCOHORT_END_YR - SCREEN_BIRTHCOHORT_START_YR +1.0)) // insured
					|| (_vecPt_BirthCohortScr[k]->_insurStatus == insr_uninsured && rnd < argModelParam._disBurdnData._prob_birthCohortScr_uninsured / (double)(SCREEN_BIRTHCOHORT_END_YR - SCREEN_BIRTHCOHORT_START_YR + 1.0))) {//uninsured

					if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_screen_accept) { // percentage of patients who accepted screening and received correct results
						_vecPt_BirthCohortScr[k]->_flagAware = true;
						_vecPt_BirthCohortScr[k]->_flagAwareByScreening = true;
						_vecPt_BirthCohortScr[k]->_curTxResps = trStatus_naive;
						_vecPt_BirthCohortScr[k]->_yrMostRecentScr = _curYear;
						_num_aware_by_screening++;
						_num_screening++;
					}
				}
			} // for _vecPtBirthCohortScr
		}// end if FIXED_CAP

	}// end if SCREENING PERIODS
	return 0;
}

int BurdenModelSim::Screening_Universal_Rate(const modelParamType & argModelParam)
{
	double pScreening = argModelParam._disBurdnData._table_universal_scr_rate.lower_bound(_curYear - START_YR_SCREEN)->second;
	map<typeCohortDisBurdnModel, double> pScreeningByCohort = argModelParam._nonNHANESPplData._scrRate_by_year_cohort.lower_bound(_curYear)->second;

	if(pScreening>EPSILON)cout<<setprecision(3)<<"    -- Screening: Universal screening at year "<<_curYear<<" - with screening rate = "<<pScreening<<" ********"<<endl;

	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		if (!(*ptIter)->_flagAware && (*ptIter)->_state <= state_CoCirr) {
			if ((*ptIter)->_flagCured) { ExitWithMsg("Screening_UsualCare(): _flagCured is true"); }

			if(argModelParam._disBurdnData._flag_multipleCohort){
				pScreening = pScreeningByCohort[(*ptIter)->_cohort];
			}
			// corrected on 4/23/2017 by Qiushi
			//if (_rnd_pt.GetU01() < argModelParam._disBurdnData._table_scr_ratio.lower_bound(_curYear - START_YR_SCREEN)->second) {
			if (_rnd_pt.GetU01() < pScreening) {
				_num_screening++;

				if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_screen_accept) {
					(*ptIter)->_flagAware = true;
					(*ptIter)->_flagAwareByScreening = true;
					(*ptIter)->_curTxResps = trStatus_naive;
					(*ptIter)->_yrMostRecentScr = _curYear;
					_num_aware_by_screening++;
				}
			}

		}
	}
	return 0;
}


int BurdenModelSim::Screening_Universal_Rate(const modelParamType & argModelParam, typeCohortDisBurdnModel argCohort)
{
	
	map<typeCohortDisBurdnModel, double> pScreeningByCohort = argModelParam._nonNHANESPplData._scrRate_by_year_cohort.lower_bound(_curYear)->second;
	double pScreening = pScreeningByCohort[argCohort];

	if (pScreening>EPSILON)cout << setprecision(3) << "    -- Screening: Universal screening at year " << _curYear << " - with screening rate = " << pScreening << " ********" << endl;

	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		if ((*ptIter)->_cohort != argCohort) continue;

		if (!(*ptIter)->_flagAware && (*ptIter)->_state <= state_CoCirr) {
			if ((*ptIter)->_flagCured) { ExitWithMsg("Screening_UsualCare(): _flagCured is true"); }

			if (_rnd_pt.GetU01() < pScreening) {
				_num_screening++;

				if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_screen_accept) {
					(*ptIter)->_flagAware = true;
					(*ptIter)->_flagAwareByScreening = true;
					(*ptIter)->_curTxResps = trStatus_naive;
					(*ptIter)->_yrMostRecentScr = _curYear;
					_num_aware_by_screening++;
				}
			}
		}
	}
	return 0;
}

int BurdenModelSim::Screening_Universal_Capacity(const modelParamType & argModelParam)
{
	if(argModelParam._disBurdnData._flag_multipleCohort){
		ExitWithMsg("Error @ Screening_universal_capacity(): not implemented for multiple cohorts yet...[5/21/2017]");
	}

	int nScreening = argModelParam._disBurdnData._table_universal_scr_cap.lower_bound(_curYear - START_YR_SCREEN)->second;
	cout << setprecision(3) << "    -- Screening: Universal screening at year " << _curYear << " - with screening capacity = " << nScreening << ", among " ;
	_num_screening = 0;

	vector<patientType *> vecPt;
	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		if (!(*ptIter)->_flagAware && (*ptIter)->_state <= state_CoCirr) {
			if ((*ptIter)->_flagCured) { ExitWithMsg("Screening_UsualCare(): _flagCured is true"); }
			vecPt.push_back(*ptIter);
		}
	}

	cout << vecPt.size() << " screening candidates" << endl;
	if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
		_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
	}
	else {
		std::random_shuffle(vecPt.begin(), vecPt.end());
	}

	for (int k = 0; (k < (int)vecPt.size() - 1) && (_num_screening < nScreening); k++) {
			if (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_screen_accept) { // percentage of patients who accepted screening and received correct results
				vecPt[k]->_flagAware = true;
				vecPt[k]->_flagAwareByScreening = true;
				vecPt[k]->_curTxResps = trStatus_naive;
				vecPt[k]->_yrMostRecentScr = _curYear;
				_num_aware_by_screening++;
				_num_screening++;
			}

	} // for vecPt

	return 0;
}



int BurdenModelSim::Treatment_DetermineElig(const modelParamType & argModelParam)
{
	// FUNCTION: determine eligbility and what treatment to take


	// identify treatment eligible patients
	_vecPtTxElig_F0F2.clear();
	_vecPtTxElig_F3F4.clear();
	_vecPtTxElig_all.clear();
	// added 2017/5/21 for other cohorts except default cohorts
	_vecPtTxElig_F0F2_nonDefaultCohort.clear();
	_vecPtTxElig_F3F4_nonDefaultCohort.clear();
	_vecPtTxElig_all_nonDefaultCohort.clear();


	if (_curYear < TX_START_YR) {
		return 0;
	}

	int ptCounter = 0;
	int listSize = (int)_listAllPts.size();
	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ++ptIter) {
		if (ptCounter % (listSize / 5) == 0) {
			cout << fixed << "        - Determine treatment eligibility for " << (int)(double(ptCounter) / double(_listAllPts.size()) * 100) << " % patients..." << endl;
		}

		ptCounter++;
		
		patientType * thePatientPointer = *ptIter;
		// -------- set default value ----------------------------
		thePatientPointer->_flagTxEligible = false;
		thePatientPointer->_txCat = txCat_NoTx;
		thePatientPointer->_flagDelayTx = false;
		// --------------------------------------------------------

		// --- check ---
		if (thePatientPointer->genotype != 1 && thePatientPointer->_curTxResps == trStatus_failed_PI1) {
			ExitWithMsg("Error @ BurdenModelSim::Treatment_DetermineElig(): G2-6 patients cannot have trStatus_failed_PI1 tx response [genotype =  "
				+ basicToStr(thePatientPointer->genotype) + ", tx response = " + basicToStr(thePatientPointer->_curTxResps));
		}

		/******************* eligibility: <=F4, aware, not-cured, >=Tx_START_YR, insured ******************/
		if (thePatientPointer->_state > state_CoCirr
			|| (!thePatientPointer->_flagAware)
			|| thePatientPointer->_flagCured
			|| thePatientPointer->_insurStatus == insr_uninsured
			|| (!thePatientPointer->_flagTrInsur)) {
			continue;
		}

		/******************* eligibility: # of previous treatment  ******************/
		
		if (RECOVER_SALVAGE_POSTER_VERSION) {
			if (thePatientPointer->_cycleOnDrug_PEGRBV + thePatientPointer->_cycleOnDrug_PI >= MAX_NUM_TX_PEGRBV + MAX_NUM_TX_PI
				|| thePatientPointer->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA) {
				continue;
			}
		}
		else {
			//if (thePatientPointer->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV
			//	|| thePatientPointer->_cycleOnDrug_PI >= MAX_NUM_TX_PI
			//	|| thePatientPointer->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA) {
			//	continue;
			//}

			// [QC comment] Modified on 03/09/2017
			if ((_curYear < START_YR_PI1 && thePatientPointer->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)
				|| ( _curYear >= START_YR_PI1 && _curYear < START_YR_DAA1 && thePatientPointer->_cycleOnDrug_PI >= MAX_NUM_TX_PI)
				|| ( _curYear >= START_YR_DAA1 && thePatientPointer->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA)) {
				continue;
			}
		}
		/******************************************************************************************************/
		/*******************      Prior to Wave 1, patients can delay treatment with prob     *****************/
		/******************************************************************************************************/
		
		if (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1 /*&& thePatientPointer->_cycleOnDrug_PEGRBV <= MAX_NUM_TX_ELIGIBLE_FOR_WAITING_PREWAVE1*/) {
			double rndDelay = _rnd_pt.GetU01();
			if (thePatientPointer->_state <= state_F2 && rndDelay < argModelParam._disBurdnData._prob_delay_tx_until_wave1_f0f2) {
				thePatientPointer->_flagDelayTx = true;
				continue;
			}
			else if (thePatientPointer->_state == state_F3 && rndDelay < argModelParam._disBurdnData._prob_delay_tx_until_wave1_f3) {
				thePatientPointer->_flagDelayTx = true;
				continue;
			}
		}
		

		if (_curYear >= START_YR_PEGRBV && _curYear < START_YR_PI1) {
			/******************************************************************************************************/
			/***********************************  before availability of triple therapy ***************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive) {
				thePatientPointer->_flagTxEligible = true;
			}
			else if (thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null) {
				thePatientPointer->_flagTxEligible = true;
			}
			else if (thePatientPointer->_curTxResps == trStatus_contraInd_mod) {
				//if (thePatientPointer->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PRETRIPLE_contraind_mod) {
				//	thePatientPointer->_flagTxEligible = false;
				//	continue;
				//}
				if (thePatientPointer->genotype == 1
					&& thePatientPointer->_state >= state_F0
					&& thePatientPointer->_state <= state_F2
					) {
					thePatientPointer->_flagTxEligible = false;
					continue;
				}
				else {
					thePatientPointer->_flagTxEligible = true;
				}
			}
			else {
				// remain ineligible for any other txRespType
			}
			thePatientPointer->_txCat = (thePatientPointer->_flagTxEligible ? txCat_PEGRBV : txCat_NoTx);

		}
		else if (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1) {
			/******************************************************************************************************/
			/****************************************** triple thearpy   ******************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null) {
				thePatientPointer->_flagTxEligible = true;
			}
			else if (thePatientPointer->_curTxResps == trStatus_contraInd_mod /*&& thePatientPointer->_cycleOnDrug_PEGRBV < MAX_NUM_TX_CONTRA_FOR_PI1*/) {
				if (thePatientPointer->genotype == 1
					&& thePatientPointer->_state >= state_F3
					&& thePatientPointer->_state <= state_CoCirr) {
					thePatientPointer->_flagTxEligible = true;
				}
				else {
					thePatientPointer->_flagTxEligible = false;
					continue;
				}
			}

			// if eligible for treatment, PI only for G1
			if (thePatientPointer->_flagTxEligible) {
				if (thePatientPointer->genotype == 1) {
					thePatientPointer->_txCat = txCat_PI1;
				}
				else {
					thePatientPointer->_txCat = txCat_PEGRBV; // G2, G3, G456 patients still recieve PEGRBV treatment
				}
			}
			else {
				thePatientPointer->_txCat = txCat_NoTx;
			}
		}

		else if (_curYear >= START_YR_DAA1 && _curYear < START_YR_DAA2) {
			/******************************************************************************************************/
			/****************************************** DAA 1 (non NS5A) ******************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null
				|| thePatientPointer->_curTxResps == trStatus_contraInd_mod
				|| thePatientPointer->_curTxResps == trStatus_contraInd_nonmod
				|| thePatientPointer->_curTxResps == trStatus_failed_PI1) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat =  txCat_DAA1_nonNS5A;
			}
			else {
				ExitWithMsg("Error @ Treatment_DetermineElig(): Wrong _curTxResps for DAA1 era.");
			}
			

		}

		else if (_curYear >= START_YR_DAA2 && _curYear < START_YR_DAA3) {
			/******************************************************************************************************/
			/****************************************** DAA 2 (NS5A + nonNS5A) ************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_contraInd_nonmod
				|| thePatientPointer->_curTxResps == trStatus_contraInd_mod		
				) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find(thePatientPointer->genotype)->second
					? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
			}
			else if ( thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null						
				|| thePatientPointer->_curTxResps == trStatus_failed_PI1
				) {
					thePatientPointer->_flagTxEligible = true;					
					if(RECOVER_SALVAGE_POSTER_VERSION){
						thePatientPointer->_txCat = txCat_DAA2_NS5A; // [poster version]
					}else{
						thePatientPointer->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find(thePatientPointer->genotype)->second
							? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A; 

					}
			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA1_nonNS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA2_nonNS5A) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = txCat_DAA2_NS5A;

			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA2_NS5A) {
				if (thePatientPointer->_state == state_CoCirr) {
					thePatientPointer->_flagTxEligible = true;
					thePatientPointer->_txCat = txCat_DAA2_nonNS5A; 					
				}
				else {
					// non-cirrhotic will wait until 2018 = DAA3 for the new NS5A reigmens
					thePatientPointer->_flagTxEligible = false;
					thePatientPointer->_txCat = txCat_NoTx;
				}
			}
		}
		else if (_curYear >= START_YR_DAA3) {
			/******************************************************************************************************/
			/****************************************** DAA 3 (NS5A + nonNS5A) ************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null
				|| thePatientPointer->_curTxResps == trStatus_contraInd_mod
				|| thePatientPointer->_curTxResps == trStatus_contraInd_nonmod
				|| thePatientPointer->_curTxResps == trStatus_failed_PI1
				) {
				thePatientPointer->_flagTxEligible = true;

				thePatientPointer->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find(thePatientPointer->genotype)->second
					? txCat_DAA3_NS5A : txCat_DAA3_nonNS5A;
			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA1_nonNS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA2_nonNS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA3_nonNS5A) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = txCat_DAA3_NS5A;
			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA2_NS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA3_NS5A) {
					thePatientPointer->_flagTxEligible = true;
					if (RECOVER_SALVAGE_POSTER_VERSION) {
						thePatientPointer->_txCat = txCat_DAA3_nonNS5A;
					}else{
						thePatientPointer->_txCat = txCat_DAA3_NS5A;
					}
			}

		}		
		else {
			ExitWithMsg("Error: TreatPatients(): uknown category of current year [=" + basicToStr(_curYear) + "].");
		}// end if identifying waves


		 /******************* ******************/

		if (thePatientPointer->_flagTxEligible && _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_receivable) {

			// --- check ---
			if (!thePatientPointer->_flagTrInsur) ExitWithMsg("Error @ Treatment_DetermineElig(): patient's _falgTrInsur = false");

			if (thePatientPointer->_state <= state_F2) {
				_vecPtTxElig_F0F2.push_back(*ptIter);
				_vecPtTxElig_all.push_back(*ptIter);
				if(thePatientPointer->_cohort != cohort_defaultSingleCohort){
					_vecPtTxElig_F0F2_nonDefaultCohort.push_back(*ptIter);
					_vecPtTxElig_all_nonDefaultCohort.push_back(*ptIter);
				}
			}
			else {
				_vecPtTxElig_F3F4.push_back(*ptIter);
				_vecPtTxElig_all.push_back(*ptIter);
				if(thePatientPointer->_cohort != cohort_defaultSingleCohort){
					_vecPtTxElig_F3F4_nonDefaultCohort.push_back(*ptIter);
					_vecPtTxElig_all_nonDefaultCohort.push_back(*ptIter);
				}
			}
		}



	}// end for all patients
	return 0;
}



bool BurdenModelSim::Treatment_DetermineElig(const modelParamType & argModelParam, patientType * thePatientPointer)
{
	// FUNCTION: determine eligbility and what treatment to take

		// -------- set default value ----------------------------
		thePatientPointer->_flagTxEligible = false;
		thePatientPointer->_txCat = txCat_NoTx;
		thePatientPointer->_flagDelayTx = false;
		// --------------------------------------------------------

		// --- check ---
		if (thePatientPointer->genotype != 1 && thePatientPointer->_curTxResps == trStatus_failed_PI1) {
			ExitWithMsg("Error @ BurdenModelSim::Treatment_DetermineElig(): G2-6 patients cannot have trStatus_failed_PI1 tx response [genotype =  "
				+ basicToStr(thePatientPointer->genotype) + ", tx response = " + basicToStr(thePatientPointer->_curTxResps));
		}

		/******************* eligibility: <=F4, aware, not-cured, >=Tx_START_YR, insured ******************/
		if (thePatientPointer->_state > state_CoCirr
			|| (!thePatientPointer->_flagAware)
			|| thePatientPointer->_flagCured
			|| thePatientPointer->_insurStatus == insr_uninsured
			|| (!thePatientPointer->_flagTrInsur)) {
			return false;
		}

		/******************* eligibility: # of previous treatment  ******************/

		if (RECOVER_SALVAGE_POSTER_VERSION) {
			if (thePatientPointer->_cycleOnDrug_PEGRBV + thePatientPointer->_cycleOnDrug_PI >= MAX_NUM_TX_PEGRBV + MAX_NUM_TX_PI
				|| thePatientPointer->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA) {
				return false;
			}
		}
		else {
			//if (thePatientPointer->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV
			//	|| thePatientPointer->_cycleOnDrug_PI >= MAX_NUM_TX_PI
			//	|| thePatientPointer->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA) {
			//	continue;
			//}

			// [QC comment] Modified on 03/09/2017
			if ((_curYear < START_YR_PI1 && thePatientPointer->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)
				|| (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1 && thePatientPointer->_cycleOnDrug_PI >= MAX_NUM_TX_PI)
				|| (_curYear >= START_YR_DAA1 && thePatientPointer->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA)) {
				return false;
			}
		}
		/******************************************************************************************************/
		/*******************      Prior to Wave 1, patients can delay treatment with prob     *****************/
		/******************************************************************************************************/

		if (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1 /*&& thePatientPointer->_cycleOnDrug_PEGRBV <= MAX_NUM_TX_ELIGIBLE_FOR_WAITING_PREWAVE1*/) {
			double rndDelay = _rnd_pt.GetU01();
			if (thePatientPointer->_state <= state_F2 && rndDelay < argModelParam._disBurdnData._prob_delay_tx_until_wave1_f0f2) {
				thePatientPointer->_flagDelayTx = true;
				return false;
			}
			else if (thePatientPointer->_state == state_F3 && rndDelay < argModelParam._disBurdnData._prob_delay_tx_until_wave1_f3) {
				thePatientPointer->_flagDelayTx = true;
				return false;
			}
		}


		if (_curYear >= START_YR_PEGRBV && _curYear < START_YR_PI1) {
			/******************************************************************************************************/
			/***********************************  before availability of triple therapy ***************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive) {
				thePatientPointer->_flagTxEligible = true;
			}
			else if (thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null) {
				thePatientPointer->_flagTxEligible = true;
			}
			else if (thePatientPointer->_curTxResps == trStatus_contraInd_mod) {
				//if (thePatientPointer->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PRETRIPLE_contraind_mod) {
				//	thePatientPointer->_flagTxEligible = false;
				//	continue;
				//}
				if (thePatientPointer->genotype == 1
					&& thePatientPointer->_state >= state_F0
					&& thePatientPointer->_state <= state_F2
					) {
					thePatientPointer->_flagTxEligible = false;
					return false;
				}
				else {
					thePatientPointer->_flagTxEligible = true;
				}
			}
			else {
				// remain ineligible for any other txRespType
			}
			thePatientPointer->_txCat = (thePatientPointer->_flagTxEligible ? txCat_PEGRBV : txCat_NoTx);

		}
		else if (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1) {
			/******************************************************************************************************/
			/****************************************** triple thearpy   ******************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null) {
				thePatientPointer->_flagTxEligible = true;
			}
			else if (thePatientPointer->_curTxResps == trStatus_contraInd_mod /*&& thePatientPointer->_cycleOnDrug_PEGRBV < MAX_NUM_TX_CONTRA_FOR_PI1*/) {
				if (thePatientPointer->genotype == 1
					&& thePatientPointer->_state >= state_F3
					&& thePatientPointer->_state <= state_CoCirr) {
					thePatientPointer->_flagTxEligible = true;
				}
				else {
					thePatientPointer->_flagTxEligible = false;
					return false;
				}
			}

			// if eligible for treatment, PI only for G1
			if (thePatientPointer->_flagTxEligible) {
				if (thePatientPointer->genotype == 1) {
					thePatientPointer->_txCat = txCat_PI1;
				}
				else {
					thePatientPointer->_txCat = txCat_PEGRBV; // G2, G3, G456 patients still recieve PEGRBV treatment
				}
			}
			else {
				thePatientPointer->_txCat = txCat_NoTx;
			}
		}

		else if (_curYear >= START_YR_DAA1 && _curYear < START_YR_DAA2) {
			/******************************************************************************************************/
			/****************************************** DAA 1 (non NS5A) ******************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null
				|| thePatientPointer->_curTxResps == trStatus_contraInd_mod
				|| thePatientPointer->_curTxResps == trStatus_contraInd_nonmod
				|| thePatientPointer->_curTxResps == trStatus_failed_PI1) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = txCat_DAA1_nonNS5A;
			}
			else {
				ExitWithMsg("Error @ Treatment_DetermineElig(): Wrong _curTxResps for DAA1 era.");
			}


		}

		else if (_curYear >= START_YR_DAA2 && _curYear < START_YR_DAA3) {
			/******************************************************************************************************/
			/****************************************** DAA 2 (NS5A + nonNS5A) ************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_contraInd_nonmod
				|| thePatientPointer->_curTxResps == trStatus_contraInd_mod
				) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find(thePatientPointer->genotype)->second
					? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
			}
			else if (thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null
				|| thePatientPointer->_curTxResps == trStatus_failed_PI1
				) {
				thePatientPointer->_flagTxEligible = true;
				if (RECOVER_SALVAGE_POSTER_VERSION) {
					thePatientPointer->_txCat = txCat_DAA2_NS5A; // [poster version]
				}
				else {
					thePatientPointer->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find(thePatientPointer->genotype)->second
						? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;

				}
			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA1_nonNS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA2_nonNS5A) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = txCat_DAA2_NS5A;

			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA2_NS5A) {
				if (thePatientPointer->_state == state_CoCirr) {
					thePatientPointer->_flagTxEligible = true;
					thePatientPointer->_txCat = txCat_DAA2_nonNS5A;
				}
				else {
					// non-cirrhotic will wait until 2018 = DAA3 for the new NS5A reigmens
					thePatientPointer->_flagTxEligible = false;
					thePatientPointer->_txCat = txCat_NoTx;
				}
			}
		}
		else if (_curYear >= START_YR_DAA3) {
			/******************************************************************************************************/
			/****************************************** DAA 3 (NS5A + nonNS5A) ************************************/
			/******************************************************************************************************/
			if (thePatientPointer->_curTxResps == trStatus_naive
				|| thePatientPointer->_curTxResps == trStatus_relap
				|| thePatientPointer->_curTxResps == trStatus_partial
				|| thePatientPointer->_curTxResps == trStatus_null
				|| thePatientPointer->_curTxResps == trStatus_contraInd_mod
				|| thePatientPointer->_curTxResps == trStatus_contraInd_nonmod
				|| thePatientPointer->_curTxResps == trStatus_failed_PI1
				) {
				thePatientPointer->_flagTxEligible = true;

				thePatientPointer->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find(thePatientPointer->genotype)->second
					? txCat_DAA3_NS5A : txCat_DAA3_nonNS5A;
			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA1_nonNS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA2_nonNS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA3_nonNS5A) {
				thePatientPointer->_flagTxEligible = true;
				thePatientPointer->_txCat = txCat_DAA3_NS5A;
			}
			else if (thePatientPointer->_curTxResps == trStatus_failed_DAA2_NS5A
				|| thePatientPointer->_curTxResps == trStatus_failed_DAA3_NS5A) {
				thePatientPointer->_flagTxEligible = true;
				if (RECOVER_SALVAGE_POSTER_VERSION) {
					thePatientPointer->_txCat = txCat_DAA3_nonNS5A;
				}
				else {
					thePatientPointer->_txCat = txCat_DAA3_NS5A;
				}
			}

		}
		else {
			ExitWithMsg("Error: TreatPatients(): uknown category of current year [=" + basicToStr(_curYear) + "].");
		}// end if identifying waves


		 /******************* ******************/


	return true;
}


int BurdenModelSim::Treatment_DetermineElig_w_ExtendedUsePRandPI(const modelParamType & argModelParam)
{
	// FUNCTION: determine eligbility and what treatment to take


	// identify treatment eligible patients
	_vecPtTxElig_F0F2.clear();
	_vecPtTxElig_F3F4.clear();
	_vecPtTxElig_all.clear();

	if (_curYear < TX_START_YR) {
		return 0;
	}

	int ptCounter = 0;
	int listSize = (int)_listAllPts.size();
	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		ptCounter++;
		if (ptCounter % ( listSize/ 5) == 0) {
			cout << fixed << "        - Determine treatment eligiblity for " << (int)(double(ptCounter) / double(_listAllPts.size()) * 100) << " % patients..." << endl;
		}

		// -------- set default value ----------------------------
		(*ptIter)->_flagTxEligible = false;
		(*ptIter)->_txCat = txCat_NoTx;
		(*ptIter)->_flagDelayTx = false;
		// --------------------------------------------------------

		// --- check ---
		if ((*ptIter)->genotype != 1 && (*ptIter)->_curTxResps == trStatus_failed_PI1) {
			ExitWithMsg("Error @ BurdenModelSim::Treatment_DetermineElig(): G2-6 patients cannot have trStatus_failed_PI1 tx response [genotype =  "
				+ basicToStr((*ptIter)->genotype) + ", tx response = " + basicToStr((*ptIter)->_curTxResps));
		}

		/******************* eligibility: <=F4, aware, not-cured, >=Tx_START_YR, insured ******************/
		if ((*ptIter)->_state > state_CoCirr
			|| (!(*ptIter)->_flagAware)
			|| (*ptIter)->_flagCured
			|| (*ptIter)->_insurStatus == insr_uninsured
			|| (!(*ptIter)->_flagTrInsur)) {
			continue;
		}

		/******************* eligibility: # of previous treatment  ******************/

		if (RECOVER_SALVAGE_POSTER_VERSION) {
			if ((*ptIter)->_cycleOnDrug_PEGRBV + (*ptIter)->_cycleOnDrug_PI >= MAX_NUM_TX_PEGRBV + MAX_NUM_TX_PI
				|| (*ptIter)->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA) {
				continue;
			}
		}
		else {
			if ((_curYear < START_YR_PI1 && (*ptIter)->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)
				|| (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1 && (*ptIter)->_cycleOnDrug_PI >= MAX_NUM_TX_PI)
				|| (_curYear >= START_YR_DAA1 && (*ptIter)->_cycleOnDrug_DAA >= MAX_NUM_TX_DAA)) {
				continue;
			}
		}
		/******************************************************************************************************/
		/*******************      Prior to Wave 1, patients can delay treatment with prob     *****************/
		/******************************************************************************************************/

		if (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1 /*&& (*ptIter)->_cycleOnDrug_PEGRBV <= MAX_NUM_TX_ELIGIBLE_FOR_WAITING_PREWAVE1*/) {
			double rndDelay = _rnd_pt.GetU01();
			if ((*ptIter)->_state <= state_F2 && rndDelay < argModelParam._disBurdnData._prob_delay_tx_until_wave1_f0f2) {
				(*ptIter)->_flagDelayTx = true;
				continue;
			}
			else if ((*ptIter)->_state == state_F3 && rndDelay < argModelParam._disBurdnData._prob_delay_tx_until_wave1_f3) {
				(*ptIter)->_flagDelayTx = true;
				continue;
			}
		}


		if (_curYear < START_YR_PI1) {
			/******************************************************************************************************/
			/***********************************  before availability of triple therapy ***************************/
			/******************************************************************************************************/
			if ((*ptIter)->_curTxResps == trStatus_naive) {
				(*ptIter)->_flagTxEligible = true;
			}
			else if ((*ptIter)->_curTxResps == trStatus_relap
				|| (*ptIter)->_curTxResps == trStatus_partial
				|| (*ptIter)->_curTxResps == trStatus_null) {
				(*ptIter)->_flagTxEligible = true;
			}
			else if ((*ptIter)->_curTxResps == trStatus_contraInd_mod) {
				//if ((*ptIter)->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PRETRIPLE_contraind_mod) {
				//	(*ptIter)->_flagTxEligible = false;
				//	continue;
				//}
				if ((*ptIter)->genotype == 1
					&& (*ptIter)->_state >= state_F0
					&& (*ptIter)->_state <= state_F2
					) {
					(*ptIter)->_flagTxEligible = false;
					continue;
				}
				else {
					(*ptIter)->_flagTxEligible = true;
				}
			}
			else {
				// remain ineligible for any other txRespType
			}
			(*ptIter)->_txCat = ((*ptIter)->_flagTxEligible ? txCat_PEGRBV : txCat_NoTx);

		}
		else if (_curYear >= START_YR_PI1 && _curYear < START_YR_DAA1) {
			/******************************************************************************************************/
			/****************************************** triple thearpy   ******************************************/
			/******************************************************************************************************/
			if ((*ptIter)->_curTxResps == trStatus_naive
				|| (*ptIter)->_curTxResps == trStatus_relap
				|| (*ptIter)->_curTxResps == trStatus_partial
				|| (*ptIter)->_curTxResps == trStatus_null) {
				(*ptIter)->_flagTxEligible = true;
			}
			else if ((*ptIter)->_curTxResps == trStatus_contraInd_mod /*&& (*ptIter)->_cycleOnDrug_PEGRBV < MAX_NUM_TX_CONTRA_FOR_PI1*/) {
				if ((*ptIter)->genotype == 1
					&& (*ptIter)->_state >= state_F3
					&& (*ptIter)->_state <= state_CoCirr) {
					(*ptIter)->_flagTxEligible = true;
				}
				else {
					(*ptIter)->_flagTxEligible = false;
					continue;
				}
			}

			// if eligible for treatment, PI only for G1
			if ((*ptIter)->_flagTxEligible) {
				if ((*ptIter)->genotype == 1) {

					//(*ptIter)->_txCat = txCat_PI1; // commented by Qiushi @ 03/25/2017
					// ---------------------------------------------
					// [QC modified 2017/3/25]: mix of PR+PI and PR
					// ---------------------------------------------
					double sharePRandPI = (argModelParam._disBurdnData._table_pr_extended_use_pr_pi.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;
					double sharePR = (argModelParam._disBurdnData._table_pr_extended_use_pr.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;
					if (_rnd_pt.GetU01() < sharePRandPI / (sharePRandPI + sharePR)) {
						(*ptIter)->_txCat = txCat_PI1;
					}
					else {
						(*ptIter)->_txCat = txCat_PEGRBV;
					}
				}
				else {
					(*ptIter)->_txCat = txCat_PEGRBV; // G2, G3, G456 patients still recieve PEGRBV treatment
				}
			}
			else {
				(*ptIter)->_txCat = txCat_NoTx;
			}
		}

		else if (_curYear >= START_YR_DAA1 && _curYear < START_YR_DAA2) {
			/******************************************************************************************************/
			/****************************************** DAA 1 (non NS5A) ******************************************/
			/******************************************************************************************************/
			if ((*ptIter)->_curTxResps == trStatus_naive
				|| (*ptIter)->_curTxResps == trStatus_relap
				|| (*ptIter)->_curTxResps == trStatus_partial
				|| (*ptIter)->_curTxResps == trStatus_null
				|| (*ptIter)->_curTxResps == trStatus_contraInd_mod
				|| (*ptIter)->_curTxResps == trStatus_contraInd_nonmod
				|| (*ptIter)->_curTxResps == trStatus_failed_PI1) {
				(*ptIter)->_flagTxEligible = true;

				// ---------------- [COMMENTS: MODIFICATION FOR EXTENDED USE OF OLD TREATMENT] ------------------------------
				double sharePRandPI = (argModelParam._disBurdnData._table_pr_extended_use_pr_pi.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;
				double sharePR = (argModelParam._disBurdnData._table_pr_extended_use_pr.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;

				if (((*ptIter)->genotype == 1 && (*ptIter)->_cycleOnDrug_PI >= MAX_NUM_TX_PI) 
				|| ((*ptIter)->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)) {
					//(*ptIter)->_txCat = txCat_DAA1_nonNS5A;
					(*ptIter)->_txCat = abs(sharePRandPI + sharePR - 1.0) < EPSILON ? txCat_NoTx : txCat_DAA1_nonNS5A;

					if ((*ptIter)->_txCat == txCat_NoTx) {
						(*ptIter)->_flagTxEligible = false;
					}
				

				}
				else {

					// sample from market share for those who do not have "constraints"
					// ---------------------------------------------
					// [QC modified 2017/3/25]: mix of PR+PI and PR
					// ---------------------------------------------
					double pExtUse = _rnd_pt.GetU01();
					if (pExtUse < sharePRandPI ) {
						(*ptIter)->_txCat = txCat_PI1;
					}
					else if (pExtUse >= sharePRandPI && pExtUse < sharePRandPI + sharePR) {
						(*ptIter)->_txCat = txCat_PEGRBV;
					}
					else {
						//(*ptIter)->_txCat = txCat_DAA1_nonNS5A;
						// ---------------------------- for most countries, ns5a share = 0% at 2014, except Germany ----------------------------
						(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
							? txCat_DAA2_NS5A : txCat_DAA1_nonNS5A;						
					}

					// override nonmod status
					if((*ptIter)->_curTxResps == trStatus_contraInd_nonmod
						|| ((*ptIter)->_curTxResps == trStatus_contraInd_mod && (*ptIter)->_state <= state_F2)){
						//(*ptIter)->_txCat = txCat_DAA1_nonNS5A;
						
						(*ptIter)->_txCat = abs(sharePRandPI + sharePR - 1.0) < EPSILON ? txCat_NoTx : txCat_DAA1_nonNS5A;

						if ((*ptIter)->_txCat == txCat_NoTx) {
							(*ptIter)->_flagTxEligible = false;
						}

					}
				}
				// -----------------------------------------------------------------------------------------------------------
				
			}


		}

		else if (_curYear >= START_YR_DAA2 && _curYear < START_YR_DAA3) {
			/******************************************************************************************************/
			/****************************************** DAA 2 (NS5A + nonNS5A) ************************************/
			/******************************************************************************************************/
			double sharePRandPI = (argModelParam._disBurdnData._table_pr_extended_use_pr_pi.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;
			double sharePR = (argModelParam._disBurdnData._table_pr_extended_use_pr.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;

			if ((*ptIter)->_curTxResps == trStatus_naive
				|| (*ptIter)->_curTxResps == trStatus_contraInd_nonmod
				|| (*ptIter)->_curTxResps == trStatus_contraInd_mod
				) {
				(*ptIter)->_flagTxEligible = true;

				// ---------------- [COMMENTS: MODIFICATION FOR EXTENDED USE OF OLD TREATMENT] ------------------------------
				if (((*ptIter)->genotype == 1 && (*ptIter)->_cycleOnDrug_PI >= MAX_NUM_TX_PI)
					||((*ptIter)->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)) {

						if(abs(sharePRandPI + sharePR - 1.0) < EPSILON){
							(*ptIter)->_txCat = txCat_NoTx;
						}else{
							(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
								? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
						}

						if ((*ptIter)->_txCat == txCat_NoTx) {
							(*ptIter)->_flagTxEligible = false;
						}
				}
				else {

					// sample from market share for those who do not have "constraints"
					double pExtUse = _rnd_pt.GetU01();
					if (pExtUse < sharePRandPI ) {
						(*ptIter)->_txCat = txCat_PI1;
					}
					else if (pExtUse >= sharePRandPI && pExtUse < sharePRandPI + sharePR) {
						(*ptIter)->_txCat = txCat_PEGRBV;
					}
					else {
						(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
							? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
					}
					
					// exceptions
					if((*ptIter)->_curTxResps == trStatus_contraInd_nonmod
						|| ((*ptIter)->_curTxResps == trStatus_contraInd_mod && (*ptIter)->_state <= state_F2)){

							if(abs(sharePRandPI + sharePR - 1.0) < EPSILON){
								(*ptIter)->_txCat = txCat_NoTx;
							}else{
								(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
									? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
							}
							if ((*ptIter)->_txCat == txCat_NoTx) {
								(*ptIter)->_flagTxEligible = false;
							}
					}
				}
				// -----------------------------------------------------------------------------------------------------------
			}
			else if ((*ptIter)->_curTxResps == trStatus_relap
				|| (*ptIter)->_curTxResps == trStatus_partial
				|| (*ptIter)->_curTxResps == trStatus_null
				|| (*ptIter)->_curTxResps == trStatus_failed_PI1
				) {
				(*ptIter)->_flagTxEligible = true;


					// ---------------- [COMMENTS: MODIFICATION FOR EXTENDED USE OF OLD TREATMENT] ------------------------------
					if (((*ptIter)->genotype == 1 && (*ptIter)->_cycleOnDrug_PI >= MAX_NUM_TX_PI)
						|| ((*ptIter)->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)) {

							if(abs(sharePRandPI + sharePR - 1.0) < EPSILON){
								(*ptIter)->_txCat = txCat_NoTx;
							}else{

								(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
									? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
							}
							
							if ((*ptIter)->_txCat == txCat_NoTx) {
								(*ptIter)->_flagTxEligible = false;
							}
					}
					else {

						// sample from market share for those who do not have "constraints"
						double pExtUse = _rnd_pt.GetU01();
						if (pExtUse < sharePRandPI ) {
							(*ptIter)->_txCat = txCat_PI1;
						}
						else if (pExtUse >= sharePRandPI && pExtUse < sharePRandPI + sharePR) {
							(*ptIter)->_txCat = txCat_PEGRBV;
						}
						else {
							(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
								? txCat_DAA2_NS5A : txCat_DAA2_nonNS5A;
						}

					}
					// -----------------------------------------------------------------------------------------------------------
				
			}
			else if ((*ptIter)->_curTxResps == trStatus_failed_DAA1_nonNS5A
				|| (*ptIter)->_curTxResps == trStatus_failed_DAA2_nonNS5A) {
				(*ptIter)->_flagTxEligible = true;
				(*ptIter)->_txCat = txCat_DAA2_NS5A;

			}
			else if ((*ptIter)->_curTxResps == trStatus_failed_DAA2_NS5A) {
				if ((*ptIter)->_state == state_CoCirr) {
					(*ptIter)->_flagTxEligible = true;
					(*ptIter)->_txCat = txCat_DAA2_nonNS5A;
				}
				else {
					// non-cirrhotic will wait until 2018 = DAA3 for the new NS5A reigmens
					(*ptIter)->_flagTxEligible = false;
					(*ptIter)->_txCat = txCat_NoTx;
				}

			}
		}
		else if (_curYear >= START_YR_DAA3) {
			/******************************************************************************************************/
			/****************************************** DAA 3 (NS5A + nonNS5A) ************************************/
			/******************************************************************************************************/
			double sharePRandPI = (argModelParam._disBurdnData._table_pr_extended_use_pr_pi.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;
			double sharePR = (argModelParam._disBurdnData._table_pr_extended_use_pr.lower_bound(_curYear)->second).find((*ptIter)->genotype)->second;

			if ((*ptIter)->_curTxResps == trStatus_naive
				|| (*ptIter)->_curTxResps == trStatus_relap
				|| (*ptIter)->_curTxResps == trStatus_partial
				|| (*ptIter)->_curTxResps == trStatus_null
				|| (*ptIter)->_curTxResps == trStatus_contraInd_mod
				|| (*ptIter)->_curTxResps == trStatus_contraInd_nonmod
				|| (*ptIter)->_curTxResps == trStatus_failed_PI1
				) {
				(*ptIter)->_flagTxEligible = true;

				// ---------------- [COMMENTS: MODIFICATION FOR EXTENDED USE OF OLD TREATMENT] ------------------------------
				if (((*ptIter)->genotype == 1 && (*ptIter)->_cycleOnDrug_PI >= MAX_NUM_TX_PI) 
					|| ((*ptIter)->genotype != 1 && (*ptIter)->_cycleOnDrug_PEGRBV >= MAX_NUM_TX_PEGRBV)) {
						if(abs(sharePRandPI + sharePR - 1.0) < EPSILON){
							(*ptIter)->_txCat = txCat_NoTx;
						}else{

							(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
								? txCat_DAA3_NS5A : txCat_DAA3_nonNS5A;
						}

						if ((*ptIter)->_txCat == txCat_NoTx) {
							(*ptIter)->_flagTxEligible = false;
						}
						
				}
				else {
					// sample from market share for those who do not have "constraints"
					double pExtUse = _rnd_pt.GetU01();
					if (pExtUse < sharePRandPI ) {
						(*ptIter)->_txCat = txCat_PI1;
					}
					else if (pExtUse >= sharePRandPI && pExtUse < sharePRandPI + sharePR) {
						(*ptIter)->_txCat = txCat_PEGRBV;
					}
					else {

						(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
							? txCat_DAA3_NS5A : txCat_DAA3_nonNS5A;
					}
					// override nonmod status
					if((*ptIter)->_curTxResps == trStatus_contraInd_nonmod
						|| ((*ptIter)->_curTxResps == trStatus_contraInd_mod && (*ptIter)->_state <= state_F2)){

							if(abs(sharePRandPI + sharePR - 1.0) < EPSILON){
								(*ptIter)->_txCat = txCat_NoTx;
							}else{
								(*ptIter)->_txCat = _rnd_pt.GetU01() < argModelParam._disBurdnData._table_ns5a_marketshare.find(_curYear)->second.find((*ptIter)->genotype)->second
									? txCat_DAA3_NS5A : txCat_DAA3_nonNS5A;
							}
							if ((*ptIter)->_txCat == txCat_NoTx) {
								(*ptIter)->_flagTxEligible = false;
							}
					}

				}
				// -----------------------------------------------------------------------------------------------------------
			}
			else if ((*ptIter)->_curTxResps == trStatus_failed_DAA1_nonNS5A
				|| (*ptIter)->_curTxResps == trStatus_failed_DAA2_nonNS5A
				|| (*ptIter)->_curTxResps == trStatus_failed_DAA3_nonNS5A) {
				(*ptIter)->_flagTxEligible = true;
				(*ptIter)->_txCat = txCat_DAA3_NS5A;
			}
			else if ((*ptIter)->_curTxResps == trStatus_failed_DAA2_NS5A
				|| (*ptIter)->_curTxResps == trStatus_failed_DAA3_NS5A) {
				(*ptIter)->_flagTxEligible = true;
				if (RECOVER_SALVAGE_POSTER_VERSION) {
					(*ptIter)->_txCat = txCat_DAA3_nonNS5A;
				}
				else {
					(*ptIter)->_txCat = txCat_DAA3_NS5A;
				}
			}

		}
		else {
			ExitWithMsg("Error: TreatPatients(): uknown category of current year [=" + basicToStr(_curYear) + "].");
		}// end if identifying waves


		 /******************* ******************/

		if ((*ptIter)->_flagTxEligible && _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_receivable) {

			// --- check ---
			if (!(*ptIter)->_flagTrInsur) ExitWithMsg("Error @ Treatment_DetermineElig(): patient's _falgTrInsur = false");

			if ((*ptIter)->_state <= state_F2) {
				_vecPtTxElig_F0F2.push_back(*ptIter);
				_vecPtTxElig_all.push_back(*ptIter);
			}
			else {
				_vecPtTxElig_F3F4.push_back(*ptIter);
				_vecPtTxElig_all.push_back(*ptIter);
			}
		}


	}
	return 0;
}

int BurdenModelSim::Treatment_Prioritize(const modelParamType & argModelParam)
{
	int txCap = argModelParam._disBurdnData._table_tx_capacity.find(_curYear)->second;

	_vecTxCandidates.clear();

		vector<patientType*> vecPtAll, vecPtF0F2, vecPtF3F4;

		for (vector<patientType*>::const_iterator ptIter = _vecPtTxElig_all.begin();
			ptIter != _vecPtTxElig_all.end(); ++ptIter) {
			vecPtAll.push_back(*ptIter);
			if ((*ptIter)->_state <= state_F2) {
				vecPtF0F2.push_back(*ptIter);
			}
			else {
				vecPtF3F4.push_back(*ptIter);
			}

		}

		if (_curYear < START_YR_TX_PRIORITIZATION || _curYear > END_YR_TX_PRIORITIZATION) {
			// random treatment before launch of PI


			//vector<patientType*> vecPt = _vecPtTxElig_all; // commented on 5/21/2017
			vector<patientType*> vecPt = vecPtAll;


			if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
				_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
			}
			else {
				std::random_shuffle(vecPt.begin(), vecPt.end());
			}

			//VectorShuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
			//std::shuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
			//std::random_shuffle(vecPt.begin(), vecPt.end());

			//if(_curYear == 1998){
			//	for(int k=0; k<20; k++){
			//		cout<<vecPt[k]->_ptID<< "\t"<<_eng_rnd_alloc()<< endl;
			//	}
			//	cout << "--- # of eligible patients = " << _vecPtTxElig_all.size() << endl;
			//	
			//}
			for (int k = 0; (k < (int)vecPt.size()) && (k < txCap); k++) {
				_vecTxCandidates.push_back(vecPt[k]);
			}

		}// <--- if not in the prioritization periods
		else {
			/****** Treatment Priorization ************/
			// First treat all F3-F4 patients
			// Then allocate the rest treatment capacity
			/******************************************/

			//vector<patientType*> vecPt = _vecPtTxElig_F3F4; // commented on 5/21/2017
			vector<patientType*> vecPt = vecPtF3F4;

			//std::random_shuffle(vecPt.begin(), vecPt.end());
			//std::shuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
			if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
				_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
			}
			else {
				std::random_shuffle(vecPt.begin(), vecPt.end());
			}


			int k = 0;
			for (; (k < vecPt.size()) && (k < txCap); k++) {
				_vecTxCandidates.push_back(vecPt[k]);
			}

			//vecPt = _vecPtTxElig_F0F2; // commented on 5/21/2017
			vecPt = vecPtF0F2;

			//std::random_shuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
			//std::shuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
			if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
				_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
			}
			else {
				std::random_shuffle(vecPt.begin(), vecPt.end());
			}

			int j = 0;
			for (; k < txCap && j < vecPt.size(); k++) {
				_vecTxCandidates.push_back(vecPt[j]);
				j++;
			}

		}// end if prioritization periods

	
	if (_vecTxCandidates.size() > txCap) {
		ExitWithMsg("Error @ Treatment_Prioritize(): Treatment list is longer than treatment capacity!");
	}
	return 0;
}


int BurdenModelSim::Treatment_Prioritize_MultipleCohort(const modelParamType & argModelParam) {
	_vecTxCandidates.clear();
	Treatment_Prioritize_byGroup(argModelParam, tx_group_general);
	Treatment_Prioritize_byGroup(argModelParam, tx_group_indianRsv);
	Treatment_Prioritize_byGroup(argModelParam, tx_group_incarcerated);
	return 0;
}

int BurdenModelSim::Treatment_Prioritize_byGroup(const modelParamType & argModelParam, typeTxCohortGroup argTxGroup)
{
	bool flag_override_no_prioritization = false;

	vector<patientType*> vecPtAll, vecPtF0F2, vecPtF3F4;
	vecPtAll.clear();
	vecPtF0F2.clear();
	vecPtF3F4.clear();
	for (vector<patientType*>::const_iterator ptIter = _vecPtTxElig_all.begin(); ptIter != _vecPtTxElig_all.end(); ++ptIter) {
		if ((*ptIter)->_tx_group == argTxGroup) {
			// skip the patients who are not in this given treatment group

			vecPtAll.push_back(*ptIter);
			if ((*ptIter)->_state <= state_F2) {
				vecPtF0F2.push_back(*ptIter);
			}
			else {
				vecPtF3F4.push_back(*ptIter);
			}
		}
	}


	int txCap = 0;
	if (argTxGroup == tx_group_general) {
		txCap = argModelParam._disBurdnData._table_tx_capacity.find(_curYear)->second;
	}
	else if (argTxGroup == tx_group_incarcerated) {
		txCap = argModelParam._nonNHANESPplData._txCap_by_year_cohort.lower_bound(_curYear)->second.find(cohort_incarcerated)->second;
		//flag_override_no_prioritization = true;
	}
	else if (argTxGroup == tx_group_indianRsv) {
		txCap = argModelParam._nonNHANESPplData._txCap_by_year_cohort.lower_bound(_curYear)->second.find(cohort_indianReservation)->second;
	}
	else {
		ExitWithMsg("Error @ Treatment_prioritization_byGroup(): unknown treatment group = " + basicToStr((int)argTxGroup));
	}
	cout << "        Prioritize " << txCap << " (" << vecPtAll.size() << ");";



	if (_curYear < START_YR_TX_PRIORITIZATION || _curYear > END_YR_TX_PRIORITIZATION || flag_override_no_prioritization) {
		// random treatment before launch of PI
		vector<patientType*> vecPt = vecPtAll;


		if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
			_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
		}
		else {
			std::random_shuffle(vecPt.begin(), vecPt.end());
		}

		for (int k = 0; (k < (int)vecPt.size()) && (k < txCap); k++) {
			_vecTxCandidates.push_back(vecPt[k]);
		}

	}// <--- if not in the prioritization periods
	else {
		/****** Treatment Priorization ************/
		// First treat all F3-F4 patients
		// Then allocate the rest treatment capacity
		/******************************************/

		//vector<patientType*> vecPt = _vecPtTxElig_F3F4; // commented on 5/21/2017
		vector<patientType*> vecPt = vecPtF3F4;

		//std::random_shuffle(vecPt.begin(), vecPt.end());
		//std::shuffle(vecPt.begin(), vecPt.end(), _eng_rnd_alloc);
		if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
			_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
		}
		else {
			std::random_shuffle(vecPt.begin(), vecPt.end());
		}


		int k = 0;
		for (; (k < vecPt.size()) && (k < txCap); k++) {
			_vecTxCandidates.push_back(vecPt[k]);
		}

		vecPt.clear();
		vecPt = vecPtF0F2;

		if (FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION) {
			_rnd_vector_shuffler.ShuffleVector(vecPt.begin(), vecPt.end());
		}
		else {
			std::random_shuffle(vecPt.begin(), vecPt.end());
		}

		int j = 0;
		for (; k < txCap && j < vecPt.size(); k++) {
			_vecTxCandidates.push_back(vecPt[j]);
			j++;
		}

	}// end if prioritization periods

	cout << "  total (cumulative) " << _vecTxCandidates.size() << " candidates waiting for treatment"<<endl;
	
	return 0;
}


int BurdenModelSim::Treatment_Treat(const modelParamType & argModelParam)
{
	map_svr_table table_SVR;
	//table_SVR = argModelParam._disBurdnData._table_SVR_old;
	table_SVR = argModelParam._disBurdnData._table_SVR;
	
	for (vector<patientType*>::iterator pat = _vecTxCandidates.begin(); pat != _vecTxCandidates.end(); pat++) {
		if (!((*pat)->_flagTxEligible && (*pat)->_txCat != txCat_NoTx)) { ExitWithMsg("ERROR @ Treatment_treat(): this patient must be eligible AND txCat != NoTx"); }

		// adding labels
		(*pat)->_yrMostRecentTx = _curYear;
		(*pat)->_stateAtTx = (*pat)->_state;

		// treat each patient
		int ptGenotype = (*pat)->genotype;
		typeTxCategory ptTxWave = (*pat)->_txCat;
		typeTxResponseStatus ptResponse = (*pat)->_curTxResps;
		typeStateDisBurdnModel ptState = (*pat)->_state;

		
		double pSVR = table_SVR[ptGenotype][ptTxWave][ptResponse][ptState];

		// added 04/10/2017
		if (ptTxWave == txCat_DAA2_NS5A || ptTxWave == txCat_DAA3_NS5A) { 
			//cout << fixed << setprecision(3) << pSVR << endl;
			(*pat)->_flagEverReceivedNS5A = true;
		}

		
		if (pSVR < EPSILON)
		{
			ExitWithMsg("ERROR @ Treatment_treat(): undefined SVR value " + basicToStr(pSVR) + " in the svr table for genotype "
				+ basicToStr(ptGenotype) + " txCategory " + basicToStr(ptTxWave) + " tx response " + basicToStr(ptResponse) + " and state " + basicToStr(ptState));
		}

		if (_rnd_pt.GetU01() < pSVR) {
			// ****** SVR achieved *****************
			if (ptState <= state_F2) {
				(*pat)->_state = state_curedF0F2;
			}
			else if (ptState == state_F3) {
				(*pat)->_state = state_curedF3;
			}
			else if (ptState == state_CoCirr) {
				(*pat)->_state = state_curedCoCirr;
			}
			else {
				ExitWithMsg("Error @ Treatment_Treat(): incorrect value of variable - ptState = " + basicToStr(ptState));
			}

			(*pat)->_curTxResps = trStatus_svr;
			(*pat)->_flagCured = true;
		}
		else {
			(*pat)->_flagCured = false;
			// ****** SVR not achieved *************
			switch (ptTxWave) {
			case txCat_PEGRBV:
				// sample relapse/partial/null
				if ((*pat)->genotype == 1) {
					(*pat)->_curTxResps = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_respStatus_distr_expr_G1, _rnd_pt);
				}
				else {
					(*pat)->_curTxResps = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_respStatus_distr_expr_G234, _rnd_pt);
				}

				break;
			case txCat_PI1:
				(*pat)->_curTxResps = trStatus_failed_PI1;
				break;
			case txCat_DAA1_nonNS5A:
				(*pat)->_curTxResps = trStatus_failed_DAA1_nonNS5A;
				break;
			case txCat_DAA2_nonNS5A:
				(*pat)->_curTxResps = trStatus_failed_DAA2_nonNS5A;
				break;
			case txCat_DAA2_NS5A:
				(*pat)->_curTxResps = trStatus_failed_DAA2_NS5A;
				(*pat)->_flagEverFailedNS5A = true;
				break;
			case txCat_DAA3_nonNS5A:
				(*pat)->_curTxResps = trStatus_failed_DAA3_nonNS5A;
				break;
			case txCat_DAA3_NS5A:
				(*pat)->_curTxResps = trStatus_failed_DAA3_NS5A;
				(*pat)->_flagEverFailedNS5A = true;
				break;
			default:
				ExitWithMsg("ERROR @ Treatment_Treat: Unknown treatment category (ptTxWave)");
				break;
			}
		}

		if (ptTxWave == txCat_PEGRBV) {
			(*pat)->_cycleOnDrug_PEGRBV++;
		}
		else if (ptTxWave == txCat_PI1) {
			(*pat)->_cycleOnDrug_PI++;
		}
		else if (ptTxWave == txCat_DAA1_nonNS5A || ptTxWave == txCat_DAA2_nonNS5A || ptTxWave == txCat_DAA2_NS5A
			|| ptTxWave == txCat_DAA3_nonNS5A || ptTxWave == txCat_DAA3_NS5A) {
			(*pat)->_cycleOnDrug_DAA++;

		}
		else {
			ExitWithMsg("Error @ Treatment_treat(): Unknown treatment category");
		}

		// ---- treatment cost, added @ 6/5/2017 -------


	} // end of for each patient


	return 0;
}


int BurdenModelSim::RunNaturalHistory_indvPt(patientType & argPt, const modelParamType & argModelParam)
{

	argPt._state_previous = argPt._state;

	double rndTrans = _rnd_pt.GetU01();
	
	// [QC comment] modified for life table based on 5-year age group in EU5 analysis
	// CHANGE: .find() to .lower_bound()
	// double bgMort = (argPt.gender == 'M') ? argModelParam._disBurdnData._bgMort_male.find((int)argPt.currentAge)->second
	//	: argModelParam._disBurdnData._bgMort_female.find((int)argPt.currentAge)->second;

	// ----------------------------- NHANES population only ---------------------------------------------------------------
	double bgMort = (argPt.gender == 'M') ? argModelParam._disBurdnData._bgMort_male.lower_bound((int)argPt.currentAge)->second
		: argModelParam._disBurdnData._bgMort_female.lower_bound((int)argPt.currentAge)->second;

	// ----------------------------- adjustment needed if running for multiple populations ---------------------------------------------------------------
	if (argModelParam._disBurdnData._flag_multipleCohort) {		
		if (argPt._cohort == cohort_defaultSingleCohort) {
			// no adjustment
		}
		else if (argPt._cohort == cohort_incarcerated
			|| argPt._cohort == cohort_homeless
			|| argPt._cohort == cohort_nursingHome
			|| argPt._cohort == cohort_hospitalized) {
			bgMort = 1.0 - pow(1.0 - bgMort, argModelParam._nonNHANESPplData.GetNonNHANESPplData("bgMortSMR", argPt._cohort));
		}
		else if (argPt._cohort == cohort_indianReservation) {
			bgMort = (argPt.gender == 'M') ? argModelParam._disBurdnData._bgMort_AI_male.lower_bound((int)argPt.currentAge)->second
				: argModelParam._disBurdnData._bgMort_AI_female.lower_bound((int)argPt.currentAge)->second;
		}
	}

	switch (argPt._state)
	{
	case state_F0:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort
			&& rndTrans < bgMort + argModelParam._transData.pr_F0_F1) {
			argPt._state = state_F1;
		}
		else {
			// no transition
		}
		break;

	case state_F1:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort
			&& rndTrans < bgMort + argModelParam._transData.pr_F1_F2) {
			argPt._state = state_F2;
		}
		else {
			// no transition
		}
		break;

	case state_F2:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort
			&& rndTrans < bgMort + argModelParam._transData.pr_F2_F3) {
			argPt._state = state_F3;
		}
		else {
			// no transition
		}
		break;

	case state_F3: {
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort
			&& rndTrans < bgMort + argModelParam._transData.pr_F3_CoCirr) {
			argPt._state = state_CoCirr;
		}
		else if (rndTrans >= bgMort + argModelParam._transData.pr_F3_CoCirr
			&& rndTrans < bgMort + argModelParam._transData.pr_F3_CoCirr + argModelParam._transData.pr_F3_HCC) {
			argPt._state = state_HCC;
		}
		else {
			// no transition
		}
		break;
	}
	case state_CoCirr: {
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort
			&& rndTrans < bgMort + argModelParam._transData.pr_CoCirr_DeCirr) {
			argPt._state = state_DeCirr;
		}
		else if (rndTrans >= bgMort + argModelParam._transData.pr_CoCirr_DeCirr
			&& rndTrans < bgMort + argModelParam._transData.pr_CoCirr_DeCirr + argModelParam._transData.pr_CoCirr_HCC) {
			argPt._state = state_HCC;
		}
		else {
			// no transition
		}
		break;
	}

	case state_DeCirr: {
		double cumProb1 = bgMort;
		double cumProb2 = cumProb1 + argModelParam._transData.pr_DeCirr_HCC;
		double cumProb3 = cumProb2 + argModelParam._transData.pr_DeCirr_LivTr;
		double cumProb4 = cumProb3 + argModelParam._transData.pr_DeCirr_DeathLiv;


		bool flag_elig_transp; // modified on 6/29/2017 by Qiushi, redefine time-varying treatment capacity
		if (argModelParam._disBurdnData._flag_use_constant_transplant_capacity) {
			flag_elig_transp = (_counter_livTr_DC[_curModelCycle] <= LIV_TRANSP_CAP_DC);
		}
		else {
			flag_elig_transp = ((_counter_livTr_DC[_curModelCycle] + _counter_livTr_HCC[_curModelCycle]) < argModelParam._disBurdnData._table_transplant_capacity.lower_bound(_curYear)->second);
		}

		if (rndTrans < cumProb1) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= cumProb1 && rndTrans < cumProb2) {
			argPt._state = state_HCC;
		}
		else if (rndTrans >= cumProb2 && rndTrans < cumProb3 && argPt.currentAge <= LIV_TRANSP_MAX_AGE && flag_elig_transp) {
			argPt._state = state_LivTr;
			_counter_livTr_DC[_curModelCycle]++;
		}
		else if (rndTrans > cumProb3 && rndTrans < cumProb4) {
			argPt._state = state_DeathLR;
		}
		else {
			argPt._state = state_DeCirr1yrPlus;
		}
		break;
	}
	case state_DeCirr1yrPlus: {
		double cumProb1 = bgMort;
		double cumProb2 = cumProb1 + argModelParam._transData.pr_DeCirr_HCC;
		double cumProb3 = cumProb2 + argModelParam._transData.pr_DeCirr_LivTr;
		double cumProb4 = cumProb3 + argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv;

		bool flag_elig_transp; // modified on 6/29/2017 by Qiushi, redefine time-varying treatment capacity
		if (argModelParam._disBurdnData._flag_use_constant_transplant_capacity) {
			flag_elig_transp = (_counter_livTr_DC[_curModelCycle] <= LIV_TRANSP_CAP_DC);
		}
		else {
			flag_elig_transp = ((_counter_livTr_DC[_curModelCycle] + _counter_livTr_HCC[_curModelCycle]) < argModelParam._disBurdnData._table_transplant_capacity.lower_bound(_curYear)->second);
		}

		if (rndTrans < cumProb1) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= cumProb1 && rndTrans < cumProb2) {
			argPt._state = state_HCC;
		}
		else if (rndTrans >= cumProb2 && rndTrans < cumProb3 && argPt.currentAge <= LIV_TRANSP_MAX_AGE && flag_elig_transp) {
			argPt._state = state_LivTr;
			_counter_livTr_DC[_curModelCycle]++;
		}
		else if (rndTrans > cumProb3 && rndTrans < cumProb4) {
			argPt._state = state_DeathLR;
		}
		else {
			// do nothing
		}
		break;
	}
	case state_LivTr:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort && rndTrans < bgMort + argModelParam._transData.pr_LivTr_DeathLiv) {
			argPt._state = state_DeathLR;
		}
		else {
			argPt._state = state_LivTr1yrPlus;
		}
		break;

	case state_LivTr1yrPlus:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort && rndTrans < bgMort + argModelParam._transData.pr_LivTr1yrPlus_DeathLiv) {
			argPt._state = state_DeathLR;
		}
		else {
			// do nothing
		}
		break;

	case state_HCC:

		bool flag_elig_transp; // modified on 6/29/2017 by Qiushi, redefine time-varying treatment capacity
		if (argModelParam._disBurdnData._flag_use_constant_transplant_capacity) {
			flag_elig_transp = (_counter_livTr_HCC[_curModelCycle] <= LIV_TRANSP_CAP_HCC);
		}
		else {
			flag_elig_transp = ((_counter_livTr_DC[_curModelCycle] + _counter_livTr_HCC[_curModelCycle]) < argModelParam._disBurdnData._table_transplant_capacity.lower_bound(_curYear)->second);
		}

		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort && rndTrans < bgMort + argModelParam._transData.pr_HCC_LivTr
			&& argPt.currentAge <= LIV_TRANSP_MAX_AGE && flag_elig_transp) {
			argPt._state = state_LivTr;
			_counter_livTr_HCC[_curModelCycle]++;
		}
		else if (rndTrans >= bgMort + argModelParam._transData.pr_HCC_LivTr && rndTrans < bgMort + argModelParam._transData.pr_HCC_LivTr + argModelParam._transData.pr_HCC_DeathLiv) {
			argPt._state = state_DeathLR;
		}
		else {
			// do nothing
		}
		break;

	case state_curedF0F2:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else {
			// do nothing
		}
		break;

	case state_curedAcute: // same as cured F0-F2, as normal people, only background mortality
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else {
			// do nothing
		}
		break;

	case state_curedF3:
		if (rndTrans < bgMort) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= bgMort
			&& rndTrans < bgMort + argModelParam._transData.pr_F3SVR_HCC) {
			argPt._state = state_HCC;
			argPt._flagCured = false;
		}
		else {
			// no transition
		}
		break;

	case state_curedCoCirr: {
		double pRegr = argPt._cycleSVRF4 == YEAR_REGRESSION_AFTER_SVR ? argModelParam._transData.pr_regression : 0;
		double cumProb1 = bgMort;
		double cumProb2 = cumProb1 + argModelParam._transData.pr_SVR_CoCirr_DeCirr;
		double cumProb3 = cumProb2 + argModelParam._transData.pr_SVR_CoCirr_HCC;
		double cumProb4 = cumProb3 + pRegr;
		if (rndTrans < cumProb1) {
			argPt._state = state_Death;
		}
		else if (rndTrans >= cumProb1 && rndTrans < cumProb2) {
			argPt._state = state_DeCirr;
			argPt._flagCured = false;
		}
		else if (rndTrans >= cumProb2 && rndTrans < cumProb3) {
			argPt._state = state_HCC;
			argPt._flagCured = false;
		}
		else if (rndTrans >= cumProb3 && rndTrans < cumProb4) {
			argPt._state = state_curedF3;
		}
		else {
			// do nothing
			argPt._cycleSVRF4++;
		}
		break;
	}
	default:
		ExitWithMsg("Error @ Natural_history_indiv: unknown patient state = " + basicToStr(argPt._state));
		break;
	}


	// increase age
	argPt.currentAge = argPt.currentAge + 1;

	return 0;
}

int BurdenModelSim::UpdateInsurance_indvPt(patientType & argPt, const modelParamType & argModelParam)
{
	if (argPt._state == state_Death || argPt._state == state_DeathLR)
		ExitWithMsg("Error @ UpdateInsurance_indvPt: identified dead patient here...");

	if (argPt._cohort == cohort_defaultSingleCohort
		|| argPt._cohort == cohort_homeless
		|| argPt._cohort == cohort_hospitalized
		|| argPt._cohort == cohort_nursingHome
		|| argPt._cohort == cohort_indianReservation
		) {
		// those to be 65 year old, transit to medicare
		//if (((int)argPt.currentAge) == AGE_CUTOFF_MEDICARE -1 ) {
		if (((int)argPt.currentAge) == AGE_CUTOFF_MEDICARE) {
			argPt._insurStatus = insr_medicare;
			argPt._flagTrInsur = (_rnd_pt.GetU01() < argModelParam._disBurdnData._prob_tx_coverage_by_medicarePartD);
		}

		// ACA starting from 2014, expand insurance coverage
		if (_curYear >= ACA_START_YR && argPt._insurStatus == insr_uninsured) {
			double r = _rnd_pt.GetU01();
			if (r < argModelParam._disBurdnData.ACAMedicaidchange.lower_bound(_curYear)->second) {
				argPt._insurStatus = insr_medicaid;
				argPt._flagTrInsur = true;
			}
			else if (r >= argModelParam._disBurdnData.ACAMedicaidchange.lower_bound(_curYear)->second
				&& r < argModelParam._disBurdnData.ACAMedicaidchange.lower_bound(_curYear)->second + argModelParam._disBurdnData.ACAPrivatechange.lower_bound(_curYear)->second) {
				argPt._insurStatus = insr_private;
				argPt._flagTrInsur = true;
			}
		}

	}
	else if (argPt._cohort == cohort_incarcerated ) {
		// no change
	}
	else {
		ExitWithMsg("ERROR @ UpdateInsurance_indvPt(): Unknown cohort type = " + basicToStr((int)argPt._cohort));
	}


	return 0;
}






int BurdenModelSim::ResetCounters()
{
	vector<int> emptyVecInt;
	vector<double> emptyVecDouble;

	for (int k = 0; k <= END_YR - START_YR; k++) {
		emptyVecInt.push_back(0);
		emptyVecDouble.push_back(0.0);
	}

	_simCounter._counterPpl = emptyVecInt;
	_simCounter._counterHCV = emptyVecInt;
	_simCounter._counterHCV_complement = emptyVecInt;

	_simCounter._counter_byCohort_HCV[cohort_defaultSingleCohort] = emptyVecInt;
	_simCounter._counter_byCohort_HCV[cohort_homeless] = emptyVecInt;
	_simCounter._counter_byCohort_HCV[cohort_nursingHome] = emptyVecInt;
	_simCounter._counter_byCohort_HCV[cohort_hospitalized] = emptyVecInt;
	_simCounter._counter_byCohort_HCV[cohort_incarcerated] = emptyVecInt;
	_simCounter._counter_byCohort_HCV[cohort_indianReservation] = emptyVecInt;

	_simCounter._counter_byCohort_awareness[cohort_defaultSingleCohort] = emptyVecInt;
	_simCounter._counter_byCohort_awareness[cohort_homeless] = emptyVecInt;
	_simCounter._counter_byCohort_awareness[cohort_nursingHome] = emptyVecInt;
	_simCounter._counter_byCohort_awareness[cohort_hospitalized] = emptyVecInt;
	_simCounter._counter_byCohort_awareness[cohort_incarcerated] = emptyVecInt;
	_simCounter._counter_byCohort_awareness[cohort_indianReservation] = emptyVecInt;

	_simCounter._counter_byCohort_screening[cohort_defaultSingleCohort] = emptyVecInt;
	_simCounter._counter_byCohort_screening[cohort_homeless] = emptyVecInt;
	_simCounter._counter_byCohort_screening[cohort_nursingHome] = emptyVecInt;
	_simCounter._counter_byCohort_screening[cohort_hospitalized] = emptyVecInt;
	_simCounter._counter_byCohort_screening[cohort_incarcerated] = emptyVecInt;
	_simCounter._counter_byCohort_screening[cohort_indianReservation] = emptyVecInt;


	_simCounter._counter_byCohort_treatment[cohort_defaultSingleCohort] = emptyVecInt;
	_simCounter._counter_byCohort_treatment[cohort_homeless] = emptyVecInt;
	_simCounter._counter_byCohort_treatment[cohort_nursingHome] = emptyVecInt;
	_simCounter._counter_byCohort_treatment[cohort_hospitalized] = emptyVecInt;
	_simCounter._counter_byCohort_treatment[cohort_incarcerated] = emptyVecInt;
	_simCounter._counter_byCohort_treatment[cohort_indianReservation] = emptyVecInt;

	_simCounter._counter_new_hcv_incidence = emptyVecInt;
	_simCounter._counter_liverRelatedDeath = emptyVecInt;
	_simCounter._counter_SVR = emptyVecInt;

	_simCounter._counter_failedSVR = emptyVecInt;
	_simCounter._counter_failedSVR_retreatable_everFailedNS5A = emptyVecInt;
	_simCounter._counter_failedSVR_retreatable_neverFailedNS5A = emptyVecInt;
	_simCounter._counter_failedSVR_unretreatable = emptyVecInt;
	_simCounter._counter_failedSVR_everReceivedNS5A = emptyVecInt;
	_simCounter._counter_failedSVR_neverReceivedNS5A = emptyVecInt;

	_simCounter._counter_aware = emptyVecInt;
	_simCounter._counter_aware_F0F4 = emptyVecInt;
	_simCounter._counter_aware_byFib[state_F0] = emptyVecInt;
	_simCounter._counter_aware_byFib[state_F1] = emptyVecInt;
	_simCounter._counter_aware_byFib[state_F2] = emptyVecInt;
	_simCounter._counter_aware_byFib[state_F3] = emptyVecInt;
	_simCounter._counter_aware_byFib[state_CoCirr] = emptyVecInt;
	_simCounter._counter_unaware = emptyVecInt;

	_simCounter._counter_incidence_DC = emptyVecInt;
	_simCounter._counter_incidence_HCC = emptyVecInt;

	_simCounter._counterHealthState[state_curedAcute] = emptyVecInt;
	_simCounter._counterHealthState[state_F0] = emptyVecInt;
	_simCounter._counterHealthState[state_F1] = emptyVecInt;
	_simCounter._counterHealthState[state_F2] = emptyVecInt;
	_simCounter._counterHealthState[state_F3] = emptyVecInt;
	_simCounter._counterHealthState[state_CoCirr] = emptyVecInt;
	_simCounter._counterHealthState[state_DeCirr] = emptyVecInt;
	_simCounter._counterHealthState[state_DeCirr1yrPlus] = emptyVecInt;
	_simCounter._counterHealthState[state_HCC] = emptyVecInt;
	_simCounter._counterHealthState[state_LivTr] = emptyVecInt;
	_simCounter._counterHealthState[state_LivTr1yrPlus] = emptyVecInt;

	////enum typeInsurance {insr_uninsured, insr_private, insr_medicare, insr_medicaid, insr_military};
	_simCounter._counterInsurance[insr_uninsured] = emptyVecInt;
	_simCounter._counterInsurance[insr_private] = emptyVecInt;
	_simCounter._counterInsurance[insr_medicare] = emptyVecInt;
	_simCounter._counterInsurance[insr_medicaid] = emptyVecInt;
	_simCounter._counterInsurance[insr_military] = emptyVecInt;
	// added two non-NHANES population
	_simCounter._counterInsurance[insr_incarcerated] = emptyVecInt;
	_simCounter._counterInsurance[insr_indianRsv] = emptyVecInt;


	_simCounter._counterTxFailure_numTx = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byNS5A["nonNS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byNS5A["NS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byNS5A["PR"] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byCirr["Cirr"] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byCirr["nonCirr"] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byGenotype[1] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byGenotype[2] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byGenotype[3] = emptyVecInt;
	_simCounter._counterTxFailure_numTx_byGenotype[456] = emptyVecInt;

	_simCounter._counterTxFailure_numFailure = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byGenotype[1] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byGenotype[2] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byGenotype[3] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byGenotype[456] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byNS5A["nonNS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byNS5A["NS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byNS5A["PR"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byCirr["Cirr"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byCirr["nonCirr"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byCirr_NS5A["Cirr_nonNS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byCirr_NS5A["Cirr_NS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byCirr_NS5A["nonCirr_nonNS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numFailure_byCirr_NS5A["nonCirr_NS5A"] = emptyVecInt;
	_simCounter._counterTxFailure_numRetxCandidates = emptyVecInt;

	_simCounter._counter_treated_alive = emptyVecInt;

	_counter_livTr_DC = emptyVecInt;
	_counter_livTr_HCC = emptyVecInt;

	_simCounter._counter_screening = emptyVecInt;
	_simCounter._counter_txElig_all = emptyVecInt;
	_simCounter._counter_txElig_F3F4 = emptyVecInt;

	_simCounter._counter_ageDistr_xsectional_given_year.clear();
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(0, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(5, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(20, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(30, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(40, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(50, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(60, 0));
	_simCounter._counter_ageDistr_xsectional_given_year.insert(pair<double, int>(70, 0));
	
	_simCounter._counter_birthcohort_infected = emptyVecInt;
	_simCounter._counter_birthcohort_aware = emptyVecInt;
	_simCounter._counter_birthcohort_unaware = emptyVecInt;
	_simCounter._counter_birthcohort_becomingAware = emptyVecInt;
	_simCounter._counter_birthcohort_cured = emptyVecInt;
	_simCounter._counter_birthcohort_gettingTx = emptyVecInt;
	_simCounter._counter_birthcohort_infected_insured = emptyVecInt;
	_simCounter._counter_birthcohort_unaware_insured = emptyVecInt;

	////// DALY counters /////////////
	_simCounter._counter_DALY_YLD = emptyVecDouble;
	_simCounter._counter_DALY_YLL = emptyVecDouble;
	_simCounter._counter_DALY_Total = emptyVecDouble;


	_simCounter._counter_LRD_ageDistr.clear();
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(0, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(5, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(20, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(30, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(40, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(50, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(60, 0));
	_simCounter._counter_LRD_ageDistr.insert(pair<double, int>(70, 0));


	_simCounter._counter_ageDistr.clear();
	_simCounter._counter_ageDistr.insert(pair<double, int>(0, 0));
	_simCounter._counter_ageDistr.insert(pair<double, int>(20, 0));
	_simCounter._counter_ageDistr.insert(pair<double, int>(30, 0));
	_simCounter._counter_ageDistr.insert(pair<double, int>(40, 0));
	_simCounter._counter_ageDistr.insert(pair<double, int>(50, 0));
	_simCounter._counter_ageDistr.insert(pair<double, int>(60, 0));
	_simCounter._counter_ageDistr.insert(pair<double, int>(70, 0));


	/////// cost counters ////////
	_simCounter._counterCost = emptyVecDouble;
	_simCounter._counterCost_byCategory["c_treatment"] = emptyVecDouble;
	_simCounter._counterCost_byCategory["c_diagnosis"] = emptyVecDouble;
	_simCounter._counterCost_byCategory["c_natr_hist"] = emptyVecDouble;
	_simCounter._counterCost_byCategory["c_screening"] = emptyVecDouble;
	
	////enum typeInsurance {insr_uninsured, insr_private, insr_medicare, insr_medicaid, insr_military};
	_simCounter._counterCost_byInsr[insr_uninsured] = emptyVecDouble;
	_simCounter._counterCost_byInsr[insr_private] = emptyVecDouble;
	_simCounter._counterCost_byInsr[insr_medicare] = emptyVecDouble;
	_simCounter._counterCost_byInsr[insr_medicaid] = emptyVecDouble;
	_simCounter._counterCost_byInsr[insr_military] = emptyVecDouble;
	_simCounter._counterCost_byInsr[insr_incarcerated] = emptyVecDouble;
	_simCounter._counterCost_byInsr[insr_indianRsv] = emptyVecDouble;

	_simCounter._counterCost_byLivDis["DC"] = emptyVecDouble;
	_simCounter._counterCost_byLivDis["HCC"] = emptyVecDouble;
	_simCounter._counterCost_byLivDis["LivTr"] = emptyVecDouble;


	_simCounter._counter_transplant_DC = emptyVecInt;
	_simCounter._counter_transplant_HCC = emptyVecInt;

	return 0;
}

int BurdenModelSim::UpdateCounter_Treamtent(const patientType & argPt, const modelParamType & argModelParam)
{

	if (argPt._yrMostRecentScr == _curYear) {
		_simCounter._counter_screening[_curModelCycle]++;
		map<typeCohortDisBurdnModel, vector<int >>::iterator it_scrNum = _simCounter._counter_byCohort_screening.find(argPt._cohort);
		it_scrNum->second[_curModelCycle] ++;
	}

	// ---------- count screening cost --------------
	if (argPt._yrMostRecentScr == _curYear || argPt._yrDiagnosedByUsualCare == _curYear) {		
		double c_scr = CalcCost_Screening(argModelParam._costData) / argModelParam._costData.scalingCoef;
		_simCounter._counterCost[_curModelCycle] += c_scr;
		_simCounter._counterCost_byCategory["c_diagnosis"][_curModelCycle] += c_scr;
		_simCounter._counterCost_byInsr[argPt._insurStatus][_curModelCycle] += c_scr;
		
		if (argPt._yrMostRecentScr == _curYear) {
			_simCounter._counterCost_byCategory["c_screening"][_curModelCycle] += c_scr;
		}
	}

	if (argPt._flagBirthCohortForScr) {
		if (argPt._yrMostRecentScr == _curYear) {
			_simCounter._counter_birthcohort_becomingAware[_curModelCycle]++;
		}
	}

	
	// skip if the patient does not get treatment this year
	if (argPt._yrMostRecentTx != _curYear) return 0;

	// ---------------------------
	// treatment cost, added @ 6/5/2017
	// ---------------------------
	double c_tx = CalcCost_Treatment(argPt, argModelParam._costData, _curYear) / argModelParam._costData.scalingCoef;
	_simCounter._counterCost[_curModelCycle] += c_tx;
	_simCounter._counterCost_byCategory["c_treatment"][_curModelCycle] += c_tx;
	_simCounter._counterCost_byInsr[argPt._insurStatus][_curModelCycle] += c_tx;



	// ---------------------------
	// treatment number
	// ---------------------------
	_simCounter._counterTxFailure_numTx[_curModelCycle]++;
	// ---- treatment by cohort ----
	map<typeCohortDisBurdnModel, vector<int >>::iterator it_txNum = _simCounter._counter_byCohort_treatment.find(argPt._cohort);
	if (it_txNum != _simCounter._counter_byCohort_treatment.end()) {
		it_txNum->second[_curModelCycle] ++;
	}
	else {
		ExitWithMsg("Error @ UpdateCounter_Treatment(): unknown cohort type");
	}


	if (argPt._flagBirthCohortForScr) {
		_simCounter._counter_birthcohort_gettingTx[_curModelCycle]++;
	}

	// by NS5A/nonNS5A categories
	if (argPt._txCat == txCat_PEGRBV) {
		_simCounter._counterTxFailure_numTx_byNS5A["PR"][_curModelCycle]++;
	}else if(argPt._txCat == txCat_PI1 ||
		argPt._txCat == txCat_DAA1_nonNS5A ||
		argPt._txCat == txCat_DAA2_nonNS5A ||
		argPt._txCat == txCat_DAA3_nonNS5A) {
		_simCounter._counterTxFailure_numTx_byNS5A["nonNS5A"][_curModelCycle]++;
	}
	else if (argPt._txCat == txCat_DAA2_NS5A || argPt._txCat == txCat_DAA3_NS5A) {
		_simCounter._counterTxFailure_numTx_byNS5A["NS5A"][_curModelCycle]++;
	}
	else {
		ExitWithMsg("ERROR @ UpdateCounter treatment failure: patient treatment category = NoTx (_txCat = " + basicToStr((int)argPt._txCat));
	}

	// by cirrhosis state
	if (argPt._stateAtTx < state_CoCirr) {
		_simCounter._counterTxFailure_numTx_byCirr["nonCirr"][_curModelCycle]++;
	}
	else if (argPt._stateAtTx == state_CoCirr) {
		_simCounter._counterTxFailure_numTx_byCirr["Cirr"][_curModelCycle]++;
	}
	else {
		ExitWithMsg("ERROR @ UpdateCounter treatment failure: patient state_at_treatment > F4 (_stateAtTx = " + basicToStr((int)argPt._stateAtTx));
	}

	// by genotype
	_simCounter._counterTxFailure_numTx_byGenotype[argPt.genotype][_curModelCycle]++;

	// ---------------------------
	// treatment failure
	// ---------------------------
	if (argPt._state != state_curedAcute && argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
		_simCounter._counterTxFailure_numFailure[_curModelCycle]++;

		if (argPt._curTxResps == trStatus_relap || argPt._curTxResps == trStatus_partial || argPt._curTxResps == trStatus_null
			//|| argPt._curTxResps == trStatus_contraInd_mod || argPt._curTxResps == trStatus_contraInd_nonmod  // [poster version]
			) {
			_simCounter._counterTxFailure_numFailure_byGenotype[argPt.genotype][_curModelCycle]++;
			_simCounter._counterTxFailure_numFailure_byNS5A["PR"][_curModelCycle]++;

			if (argPt._stateAtTx <= state_F3) {
				_simCounter._counterTxFailure_numFailure_byCirr["nonCirr"][_curModelCycle]++;
			}
			else if (argPt._stateAtTx == state_CoCirr) {
				_simCounter._counterTxFailure_numFailure_byCirr["Cirr"][_curModelCycle]++;
			}
			else {
				ExitWithMsg("Error @  UpdateCounter treatment failure: state at treatment = " + basicToStr((int)argPt._stateAtTx));
			}

		}else if(argPt._curTxResps == trStatus_failed_PI1
			|| argPt._curTxResps == trStatus_failed_DAA1_nonNS5A
			|| argPt._curTxResps == trStatus_failed_DAA2_nonNS5A
			|| argPt._curTxResps == trStatus_failed_DAA3_nonNS5A
			) {

			_simCounter._counterTxFailure_numFailure_byGenotype[argPt.genotype][_curModelCycle]++;
			_simCounter._counterTxFailure_numFailure_byNS5A["nonNS5A"][_curModelCycle]++;

			if (argPt._stateAtTx <= state_F3) {
				_simCounter._counterTxFailure_numFailure_byCirr_NS5A["nonCirr_nonNS5A"][_curModelCycle]++;
				_simCounter._counterTxFailure_numFailure_byCirr["nonCirr"][_curModelCycle]++;
			}
			else if (argPt._stateAtTx == state_CoCirr) {

				_simCounter._counterTxFailure_numFailure_byCirr_NS5A["Cirr_nonNS5A"][_curModelCycle]++;
				_simCounter._counterTxFailure_numFailure_byCirr["Cirr"][_curModelCycle]++;
			}
			else {
				ExitWithMsg("Error @  UpdateCounter treatment failure: state at treatment = " + basicToStr((int)argPt._stateAtTx));
			}
		}
		else if (argPt._curTxResps == trStatus_failed_DAA2_NS5A
			|| argPt._curTxResps == trStatus_failed_DAA3_NS5A) {

			_simCounter._counterTxFailure_numFailure_byGenotype[argPt.genotype][_curModelCycle]++;
			_simCounter._counterTxFailure_numFailure_byNS5A["NS5A"][_curModelCycle]++;
			if (argPt._stateAtTx <= state_F3) {
				_simCounter._counterTxFailure_numFailure_byCirr_NS5A["nonCirr_NS5A"][_curModelCycle]++;
				_simCounter._counterTxFailure_numFailure_byCirr["nonCirr"][_curModelCycle]++;
			}
			else if (argPt._stateAtTx == state_CoCirr) {

				_simCounter._counterTxFailure_numFailure_byCirr_NS5A["Cirr_NS5A"][_curModelCycle]++;
				_simCounter._counterTxFailure_numFailure_byCirr["Cirr"][_curModelCycle]++;
			}
			else {
				ExitWithMsg("Error @  UpdateCounter treatment failure: state at treatment = " + basicToStr((int)argPt._stateAtTx));
			}

		}
		else {
			// skip the treatment_naive, _svr, _unknown patients.
		}
	}
	return 0;
}

int BurdenModelSim::UpdateCounter_DeathEvent(const patientType & argPt, const modelParamType & argModelParam)
{

	if (argPt._state == state_DeathLR) {
		_simCounter._counter_liverRelatedDeath[_curModelCycle]++;
		_simCounter._counter_DALY_YLL[_curModelCycle] += CalcDALY_GetYLL(argPt.currentAge, argPt.gender, argModelParam._dalyData);
		if(_curYear >= _simCounter._year_start_lrd_ageDistr && _curYear <= _simCounter._record_ageDistr_year_lrd){
			map<double, int>::iterator it = _simCounter._counter_LRD_ageDistr.lower_bound(argPt.currentAge);
			if (it != _simCounter._counter_LRD_ageDistr.begin()) { it--; }
			it->second++;
		}
	}
	return 0;
}

int BurdenModelSim::UpdateCounter_HealthState(const patientType & argPt, const modelParamType & argModelParam)
{

	_simCounter._counterPpl[_curModelCycle]++;


	
	//if (argPt._curTxResps == trStatus_svr) {
	//	if (argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
	//		ExitWithMsg("Error @ UPdateCounter_HealthState (): response is SVR but state is not curedF0F2/curedF3/curedCoCirr (state = "+basicToStr((int)argPt._state));
	//	}
	//	_simCounter._counter_SVR[_curModelCycle]++;
	//}
	//if (argPt._flagCured) {
	//	if (argPt._curTxResps != trStatus_svr) {


	if (argPt._curTxResps == trStatus_svr ) { // changed on 04/10/2017
		if (argPt._state == state_HCC || argPt._state == state_DeCirr || argPt._state == state_DeCirr1yrPlus || argPt._state == state_LivTr || argPt._state == state_LivTr1yrPlus) {
			// do nothing
		}
		else {
			if (!argPt._flagCured) {
				ExitWithMsg("Error @ UPdateCounter_HealthState (): response is SVR but _flagCured = false (except DC, LivTr, HCC states)");
			}
			if (argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
				ExitWithMsg("Error @ UPdateCounter_HealthState (): response is SVR but state is not curedF0F2/curedF3/curedCoCirr (state = " + basicToStr((int)argPt._state));
			}
		}
		_simCounter._counter_SVR[_curModelCycle]++;
	}
	// --------------- added 4/6/2017 by Qiushi ---------------
	//trStatus_naive, trStatus_svr,
	//	trStatus_relap, trStatus_partial, trStatus_null, //treatment response/history status
	//	trStatus_contraInd_mod, trStatus_contraInd_nonmod,
	//	trStatus_failed_PI1,
	//	trStatus_failed_DAA1_nonNS5A,
	//	trStatus_failed_DAA2_nonNS5A, trStatus_failed_DAA2_NS5A,
	//	trStatus_failed_DAA3_nonNS5A, trStatus_failed_DAA3_NS5A,
	else if (argPt._curTxResps == trStatus_relap
		|| argPt._curTxResps == trStatus_partial
		|| argPt._curTxResps == trStatus_null
		|| argPt._curTxResps == trStatus_failed_PI1
		|| argPt._curTxResps == trStatus_failed_DAA1_nonNS5A
		|| argPt._curTxResps == trStatus_failed_DAA2_nonNS5A
		|| argPt._curTxResps == trStatus_failed_DAA2_NS5A
		|| argPt._curTxResps == trStatus_failed_DAA3_nonNS5A
		|| argPt._curTxResps == trStatus_failed_DAA3_NS5A) {
		_simCounter._counter_failedSVR[_curModelCycle]++;

		
		if(argPt._state <= state_CoCirr) {
			_simCounter._counterTxFailure_numRetxCandidates[_curModelCycle]++;

			if (argPt._flagEverFailedNS5A) {
				_simCounter._counter_failedSVR_retreatable_everFailedNS5A[_curModelCycle]++;
			}
			else {
				_simCounter._counter_failedSVR_retreatable_neverFailedNS5A[_curModelCycle]++;
			}
		}
		else {
			_simCounter._counter_failedSVR_unretreatable[_curModelCycle]++;
		}
	}
	else {

	}


	//// for those who have ever been treated (excluding treatment naive pts)
	//if (argPt._curTxResps != trStatus_naive && argPt._curTxResps != trStatus_contraInd_mod && argPt._curTxResps != trStatus_contraInd_nonmod && argPt._curTxResps != trStatus_unknown) {
	//	if (argPt._flagEverReceivedNS5A) {
	//		_simCounter._counter_failedSVR_everReceivedNS5A[_curModelCycle]++;
	//	}
	//	else {
	//		_simCounter._counter_failedSVR_neverReceivedNS5A[_curModelCycle]++;
	//	}
	//}

	// -----------------------------------------------------------


	// for HCV positive patients (excluding the cured)
	if (argPt._state != state_curedAcute && argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
		_simCounter._counterHCV[_curModelCycle]++;

		// count cost by fibrosis state
		double c_natHistory = CalcCost_NaturalHistory(argPt, argModelParam._costData) / argModelParam._costData.scalingCoef;
		_simCounter._counterCost[_curModelCycle] += c_natHistory;
		_simCounter._counterCost_byCategory["c_natr_hist"][_curModelCycle] += c_natHistory;
		_simCounter._counterCost_byInsr[argPt._insurStatus][_curModelCycle] += c_natHistory;
		if (argPt._state == state_DeCirr || argPt._state == state_DeCirr1yrPlus) {
			_simCounter._counterCost_byLivDis["DC"][_curModelCycle] += c_natHistory;
		}
		else if (argPt._state == state_LivTr || argPt._state == state_LivTr1yrPlus) {
			_simCounter._counterCost_byLivDis["LivTr"][_curModelCycle] += c_natHistory;
		}
		else if (argPt._state == state_HCC) {
			_simCounter._counterCost_byLivDis["HCC"][_curModelCycle] += c_natHistory;
		}
		else {
			// do nothing, skip other states.
		}

		// age distribution for a given year
		if (_curYear == YEAR_RECORD_AGE_DISTR) {
			map<double, int>::iterator it = _simCounter._counter_ageDistr_xsectional_given_year.lower_bound(argPt.currentAge);
			if (it != _simCounter._counter_ageDistr_xsectional_given_year.begin()) { it--; }
			it->second++;
		}


		if (argPt._flagAware) {
			_simCounter._counter_aware[_curModelCycle]++;
			if (argPt._state >= state_F0 && argPt._state <= state_CoCirr) {
				_simCounter._counter_aware_F0F4[_curModelCycle]++;
				_simCounter._counter_aware_byFib[argPt._state][_curModelCycle]++;
			}
		}
		else {
			_simCounter._counter_unaware[_curModelCycle]++;
		}

		// --------------------------------------------------------------------------------------------------------------------------
		// ------------------------------------------------- by cohort  -------------------------------------------------------------
		// --------------------------------------------------------------------------------------------------------------------------
		map<typeCohortDisBurdnModel, vector<int>>::iterator it = _simCounter._counter_byCohort_HCV.find(argPt._cohort);
		if (it != _simCounter._counter_byCohort_HCV.end()) {
			it->second[_curModelCycle]++;

			if (argPt._flagAware) {
				map<typeCohortDisBurdnModel, vector<int >>::iterator it2 = _simCounter._counter_byCohort_awareness.find(argPt._cohort);
				it2->second[_curModelCycle] ++;
			}
		}
		else {
			ExitWithMsg("Error @ BurdenModelSim::UpdateCounter_HealthState(const patientType & argPt): Can't find argPt._cohort = " + basicToStr((int)argPt._cohort));			
		}


	}
	else {
		_simCounter._counterHCV_complement[_curModelCycle]++;

	}


	// by health states
	map<typeStateDisBurdnModel, vector<int>>::iterator it = _simCounter._counterHealthState.find(argPt._state);
	if (it != _simCounter._counterHealthState.end()) {
		it->second[_curModelCycle]++;
	}
	else {
		if (argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
			ExitWithMsg("Error @ BurdenModelSim::UpdateCounter_HealthState(const patientType & argPt): Can't find state = " + basicToStr((int)argPt._state));
		}
	}


	// counter incidence of DC/HCC
	if (argPt._state == state_DeCirr && argPt._state_previous != state_DeCirr) {
		_simCounter._counter_incidence_DC[_curModelCycle]++;
	}
	if (argPt._state == state_HCC && argPt._state_previous != state_HCC) {
		_simCounter._counter_incidence_HCC[_curModelCycle]++;
	}

	// counter for birth cohort
	if (argPt._flagBirthCohortForScr) {
		if (argPt._state != state_curedAcute && argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
			_simCounter._counter_birthcohort_infected[_curModelCycle]++;
		}
		
		if (argPt._flagAware) {
			_simCounter._counter_birthcohort_aware[_curModelCycle]++;
		}
		else {
			if (argPt._state == state_curedF0F2 || argPt._state == state_curedF3 || argPt._state == state_curedCoCirr) {
				ExitWithMsg("ERROR @ UpdateCounter_HeatlhState: patients are cured but not awared");
			}
			_simCounter._counter_birthcohort_unaware[_curModelCycle]++;

			if (argPt._insurStatus != insr_uninsured) {
				_simCounter._counter_birthcohort_unaware_insured[_curModelCycle]++;
			}

		}


		if (argPt._insurStatus != insr_uninsured) {
			_simCounter._counter_birthcohort_infected_insured[_curModelCycle]++;
		}

		if (argPt._state == state_curedAcute || argPt._state == state_curedF0F2 || argPt._state == state_curedF3 || argPt._state == state_curedCoCirr) {
			_simCounter._counter_birthcohort_cured[_curModelCycle]++;
		}
	}




	///// Update DALYs - YLD ////////
	_simCounter._counter_DALY_YLD[_curModelCycle] += CalcDALY_GetYLD(argPt.currentAge, argPt._state, 1, argModelParam._dalyData);



	// ---- record age distribution ----
	if (_curYear == argModelParam._disBurdnData._record_ageDistr_year) {
		map<double, int>::iterator it = _simCounter._counter_ageDistr.lower_bound(argPt.currentAge);
		if (it != _simCounter._counter_ageDistr.begin()) { it--; }
		it->second++;
	}


	_simCounter._counter_transplant_DC = _counter_livTr_DC;
	_simCounter._counter_transplant_HCC = _counter_livTr_HCC;

	return 0;
}

int BurdenModelSim::UpdateCounter_Insurance(const patientType & argPt)
{
	if (argPt._state != state_curedAcute && argPt._state != state_curedF0F2 && argPt._state != state_curedF3 && argPt._state != state_curedCoCirr) {
		if (_simCounter._counterInsurance.find(argPt._insurStatus) != _simCounter._counterInsurance.end()) {
			_simCounter._counterInsurance[argPt._insurStatus][_curModelCycle]++;
		}
		else {
			ExitWithMsg("ERROR @ update insurance counter: unknown insurance type = " + basicToStr((int)argPt._insurStatus));
		}
	}
	return 0;
}

int BurdenModelSim::UpdateCounter_TxEligibility()
{
	
	_simCounter._counter_txElig_all[_curModelCycle] += _vecPtTxElig_all.size();
	_simCounter._counter_txElig_F3F4[_curModelCycle] += _vecPtTxElig_F3F4.size();

	// added 05/04/2017
	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ++ptIter) {
		if ((*ptIter)->_flagCured) {
		//if ((*ptIter)->_cycleOnDrug_PEGRBV > 0 || (*ptIter)->_cycleOnDrug_PI > 0 || (*ptIter)->_cycleOnDrug_DAA > 0) {
			_simCounter._counter_treated_alive[_curModelCycle]++;
		}

		
	}

	return 0;
}

double BurdenModelSim::CalcDALY_GetYLD(double argAge, typeStateDisBurdnModel argState, double argL, const DALYType & argDALYData)
{

	double dw = CalcDALY_GetDisutilityWeight_DisBurdenModel(argState, argDALYData);
	if (argDALYData._discountRate < EPSILON) {//no time discounting
		double partWithWeight = argDALYData._C * exp(-argDALYData._beta * argAge) / pow(argDALYData._beta, 2.0) * (exp(-argDALYData._beta * argL) * (-argDALYData._beta * (argL + argAge) - 1)
			- (-argDALYData._beta * argAge - 1));
		double partWithoutWeight = argL;
		return dw * (argDALYData._K * partWithWeight + (1 - argDALYData._K) * partWithoutWeight);

	}
	else {


		double rPlusBeta = argDALYData._discountRate + argDALYData._beta;
		double partWithWeight = argDALYData._C * exp(argDALYData._discountRate * argAge) / pow(rPlusBeta, 2.0) * (exp(-rPlusBeta * (argL + argAge)) * (-rPlusBeta * (argL + argAge) - 1)
			- exp(-rPlusBeta * argAge) * (-rPlusBeta * argAge - 1));
		double partWithoutWeight = (1 - exp(-argDALYData._discountRate * argL)) / argDALYData._discountRate;
		return dw * (argDALYData._K * partWithWeight + (1 - argDALYData._K) * partWithoutWeight);

	}
}

double BurdenModelSim::CalcDALY_GetYLL(double argAge, char argGender, const DALYType & argDALYData)
{

	double L;
	if (argGender == 'M') {
		L = Interpolate(argAge, argDALYData._lifeExpct_male);
	}
	else {
		L = Interpolate(argAge, argDALYData._lifeExpct_female);
	}


	if (argDALYData._discountRate < EPSILON) {//no time discounting

		double partWithWeight = argDALYData._C * exp(-argDALYData._beta * argAge) / pow(argDALYData._beta, 2.0) * (exp(-argDALYData._beta * L) * (-argDALYData._beta * (L + argAge) - 1)
			- (-argDALYData._beta * argAge - 1));
		double partWithoutWeight = L;
		return argDALYData._K * partWithWeight + (1 - argDALYData._K) * partWithoutWeight;

	}
	else { // with time discounting
		double rPlusBeta = argDALYData._discountRate + argDALYData._beta;
		double partWithWeight = argDALYData._C * exp(argDALYData._discountRate * argAge) / pow(rPlusBeta, 2.0) * (exp(-rPlusBeta * (L + argAge)) * (-rPlusBeta * (L + argAge) - 1)
			- exp(-rPlusBeta * argAge) * (-rPlusBeta * argAge - 1));
		double partWithoutWeight = (1 - exp(-argDALYData._discountRate * L)) / argDALYData._discountRate;
		return argDALYData._K * partWithWeight + (1 - argDALYData._K) * partWithoutWeight;

	}

}

double BurdenModelSim::CalcDALY_GetDisutilityWeight_DisBurdenModel(typeStateDisBurdnModel argState, const DALYType & argDALYData)
{

	switch (argState) {
	case state_F0:	return argDALYData._dw_f0;
	case state_F1:	return argDALYData._dw_f1;
	case state_F2:	return argDALYData._dw_f2;
	case state_F3:	return argDALYData._dw_f3;
	case state_CoCirr:	return argDALYData._dw_CoCirr;
	case state_DeCirr:	return argDALYData._dw_DeCirr;
	case state_DeCirr1yrPlus:	return argDALYData._dw_DeCirr1yrPlus;
	case state_HCC:	return argDALYData._dw_HCC;
	case state_LivTr:	return argDALYData._dw_LivTr;
	case state_LivTr1yrPlus:	return argDALYData._dw_LivTr1yrPlus;
	case state_curedF0F2:	return argDALYData._dw_SVR;
	case state_curedF3:	return argDALYData._dw_SVR;
	case state_curedCoCirr:	return argDALYData._dw_SVR;
	case state_curedAcute:	return argDALYData._dw_SVR;
	default:
		return ExitWithMsg("Incorrect input state = " + basicToStr(argState) + " to get DisutilityWeight @ class DALYType");
	}
}

double BurdenModelSim::CalcCost_Treatment(const patientType & argPt, const costType & argCostData, int argYear)
{
	//if (argYear >= argCostData.override_DAA_cost_since_when) {
	//	return argCostData.override_DAA_cost_value;
	//}
	//
	// double cTxDiscount = argCostData.map_drug_cost_discount.lower_bound(argYear)->second;
	
	double cTx = argCostData._table_tx_cost.find(argYear)->second;
	return cTx;
	
}

double BurdenModelSim::CalcCost_Screening(const costType & argCostData)
{
	return argCostData.c_screening;
}

double BurdenModelSim::CalcCost_NaturalHistory(const patientType & argPt, const costType & argCostData)
{

	switch (argPt._state) {
	case state_F0: return argCostData.c_F0;
	case state_F1: return argCostData.c_F1;
	case state_F2: return argCostData.c_F2;
	case state_F3: return argCostData.c_F3;
	case state_CoCirr: return argCostData.c_CoCirr;
	case state_DeCirr: return argCostData.c_DeCirr;
	case state_DeCirr1yrPlus: return argCostData.c_DeCirr1yrPlus;
	case state_HCC: return argCostData.c_HCC;
	case state_LivTr: return argCostData.c_LivTr;
	case state_LivTr1yrPlus: return argCostData.c_PostLivTr;
	default: ExitWithMsg("Error @ CalcCost_Natural_history: not applicable for state = " + basicToStr((int)argPt._state));

	}
	return 0.0;
}


int BurdenModelSim::PrintOutputHeader(string argFile)
{

	_outf.open(argFile);
	_outf << fixed << showpoint;

	// print table head
	_outf << "Year" << "\t" << "Total#" << "\t" << "HCV#" << "\t" << "SVR" << "\t" << "Aware" << "\t" << "Unaware" << "\t" 
		<<"NewHCVIncidence"<<"\t" <<"LRD#"<<"\t";
	for (map<typeStateDisBurdnModel, vector<int>>::const_iterator it = _simCounter._counterHealthState.begin();
		it != _simCounter._counterHealthState.end(); it++) {
		_outf << typeStateDisBurdnModel_label[it->first] << "\t";
	}
	for (map<typeInsurance, vector<int>>::const_iterator it = _simCounter._counterInsurance.begin();
		it != _simCounter._counterInsurance.end(); it++) {
		_outf << typeInsurance_label[it->first] << "\t";
	}


	PrintOutputHeader_MapStringIndex(_simCounter._counterTxFailure_numTx_byNS5A);
	PrintOutputHeader_MapStringIndex(_simCounter._counterTxFailure_numTx_byCirr);
	for (map<int, vector<int>>::const_iterator it = _simCounter._counterTxFailure_numTx_byGenotype.begin();
		it != _simCounter._counterTxFailure_numTx_byGenotype.end(); it++) {
		_outf << "GT" << it->first << "\t";
	}
	
	PrintOutputHeader_MapStringIndex(_simCounter._counterTxFailure_numFailure_byNS5A);
	for (map<int, vector<int>>::const_iterator it = _simCounter._counterTxFailure_numFailure_byGenotype.begin();
		it != _simCounter._counterTxFailure_numFailure_byGenotype.end(); it++) {
		_outf << "GT" << it->first << "\t";
	}
	PrintOutputHeader_MapStringIndex(_simCounter._counterTxFailure_numFailure_byCirr);
	PrintOutputHeader_MapStringIndex(_simCounter._counterTxFailure_numFailure_byCirr_NS5A);

	_outf << "Inci_DC\tInci_HCC\t#treatment\t#txFailure\t#screening\t#txEligAll\t#txEligF3F4\t";

	_outf << "birthCohort_HCV+\tbirthCohort_HCV+_insured\tbirthCohort_aware\tbirthCohort_unaware\tbirthCohort_unaware_insured\t"
		  << "BirthCohort_becomingAware\tbirthCohort_cured\tbirthCohort_GettingTx\t";

	_outf << "aliveWhoFailedSVR\taliveEverFailedNS5A\taliveNEVERfailedNS5A\taliveUnretreatable\taliveEverReceivedNS5A\taliveNEVERreceivedNS5A\t";

	_outf << "DALY_YLD\tDALY_YLL\tDALY_Total\t";


	for (map<typeCohortDisBurdnModel, vector<int>>::const_iterator it = _simCounter._counter_byCohort_HCV.begin();
		it != _simCounter._counter_byCohort_HCV.end(); it++) {
		_outf << "HCV_"<<typeCohort_label[(int) it->first] << "\t";
	}


	for (map<typeCohortDisBurdnModel, vector<int>>::const_iterator it = _simCounter._counter_byCohort_HCV.begin();
		it != _simCounter._counter_byCohort_HCV.end(); it++) {
		_outf << "awareness_" << typeCohort_label[(int)it->first] << "\t";
	}

	for (map<typeCohortDisBurdnModel, vector<int>>::const_iterator it = _simCounter._counter_byCohort_HCV.begin();
		it != _simCounter._counter_byCohort_HCV.end(); it++) {
		_outf << "screening_" << typeCohort_label[(int)it->first] << "\t";
	}

	for (map<typeCohortDisBurdnModel, vector<int>>::const_iterator it = _simCounter._counter_byCohort_HCV.begin();
		it != _simCounter._counter_byCohort_HCV.end(); it++) {
		_outf << "treatment_" << typeCohort_label[(int)it->first] << "\t";
	}

	_outf << "cost_total\t";
	PrintOutputHeader_MapStringIndex(_simCounter._counterCost_byCategory);

	for (map<typeInsurance, vector<int>>::const_iterator it = _simCounter._counterInsurance.begin();
		it != _simCounter._counterInsurance.end(); it++) {
		_outf << "c_"<<typeInsurance_label[it->first] << "\t";
	}

	PrintOutputHeader_MapStringIndex(_simCounter._counterCost_byLivDis);


	_outf << "numTranspl\tnumReTxCandidt\t#awareF0_F4\t";
	for (map<typeStateDisBurdnModel, vector<int>>::const_iterator it = _simCounter._counter_aware_byFib.begin();
		it != _simCounter._counter_aware_byFib.end(); it++) {
		_outf << "#aware_"<< typeStateDisBurdnModel_label[it->first] << "\t";
	}
	_outf << endl;


	return 0;
}

int BurdenModelSim::PrintOutputRow(int argModelCycle)
{
	int k = argModelCycle;
	_outputRowData.clear();
	_outputRowData.push_back(_simCounter._counterPpl[k]);
	_outputRowData.push_back(_simCounter._counterHCV[k]);
	_outputRowData.push_back(_simCounter._counter_SVR[k]);
	_outputRowData.push_back(_simCounter._counter_aware[k]);
	_outputRowData.push_back(_simCounter._counter_unaware[k]);
	_outputRowData.push_back(_simCounter._counter_new_hcv_incidence[k]);
	_outputRowData.push_back(_simCounter._counter_liverRelatedDeath[k]);
	SaveOutputRow_MapTemplate(_simCounter._counterHealthState, k, _outputRowData, FLAG_SMOOTH_RESULTS);	
	SaveOutputRow_MapTemplate(_simCounter._counterInsurance, k, _outputRowData, FLAG_SMOOTH_RESULTS);

	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numTx_byNS5A, k, _outputRowData);
	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numTx_byCirr, k, _outputRowData);
	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numTx_byGenotype, k, _outputRowData);

	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numFailure_byNS5A, k, _outputRowData);
	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numFailure_byGenotype, k, _outputRowData);
	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numFailure_byCirr, k, _outputRowData);
	SaveOutputRow_MapTemplate(_simCounter._counterTxFailure_numFailure_byCirr_NS5A, k, _outputRowData);

	_outputRowData.push_back(_simCounter._counter_incidence_DC[k]);
	_outputRowData.push_back(_simCounter._counter_incidence_HCC[k]);
	_outputRowData.push_back(_simCounter._counterTxFailure_numTx[k]);
	_outputRowData.push_back(_simCounter._counterTxFailure_numFailure[k]);
	_outputRowData.push_back(_simCounter._counter_screening[k]);
	//_outputRowData.push_back(_simCounter._counter_txElig_all[k] * 0.915);
	//_outputRowData.push_back(_simCounter._counter_txElig_F3F4[k] * 0.915);
	_outputRowData.push_back(_simCounter._counter_txElig_all[k]);
	_outputRowData.push_back(_simCounter._counter_txElig_F3F4[k]);

	//_outf << _vecPtTxElig_all.size() << "\t";
	//_outf << _vecPtTxElig_F3F4.size() << "\t";

	_outputRowData.push_back(_simCounter._counter_birthcohort_infected[k]);
	_outputRowData.push_back(_simCounter._counter_birthcohort_infected_insured[k]);
	_outputRowData.push_back(_simCounter._counter_birthcohort_aware[k]);
	_outputRowData.push_back(_simCounter._counter_birthcohort_unaware[k]);
	_outputRowData.push_back(_simCounter._counter_birthcohort_unaware_insured[k]);
	//_outputRowData.push_back(_simCounter._counter_birthcohort_becomingAware[k]);
	_outputRowData.push_back(_num_screening);
	_outputRowData.push_back(_simCounter._counter_birthcohort_cured[k]);
	_outputRowData.push_back(_simCounter._counter_birthcohort_gettingTx[k]);


	_outputRowData.push_back(_simCounter._counter_failedSVR[k]);
	_outputRowData.push_back(_simCounter._counter_failedSVR_retreatable_everFailedNS5A[k]);
	_outputRowData.push_back(_simCounter._counter_failedSVR_retreatable_neverFailedNS5A[k]);
	_outputRowData.push_back(_simCounter._counter_failedSVR_unretreatable[k]);
	_outputRowData.push_back(_simCounter._counter_failedSVR_everReceivedNS5A[k]);
	_outputRowData.push_back(_simCounter._counter_failedSVR_neverReceivedNS5A[k]);

	_outputRowData.push_back(_simCounter._counter_DALY_YLD[k]);
	_outputRowData.push_back(_simCounter._counter_DALY_YLL[k]);
	_simCounter._counter_DALY_Total[k] = _simCounter._counter_DALY_YLD[k] + _simCounter._counter_DALY_YLL[k];
	_outputRowData.push_back(_simCounter._counter_DALY_Total[k]);

	SaveOutputRow_MapTemplate(_simCounter._counter_byCohort_HCV, k, _outputRowData, FLAG_SMOOTH_RESULTS);
	SaveOutputRow_MapTemplate(_simCounter._counter_byCohort_awareness, k, _outputRowData, FLAG_SMOOTH_RESULTS);
	SaveOutputRow_MapTemplate(_simCounter._counter_byCohort_screening, k, _outputRowData, FLAG_SMOOTH_RESULTS);
	SaveOutputRow_MapTemplate(_simCounter._counter_byCohort_treatment, k, _outputRowData, FLAG_SMOOTH_RESULTS);

	_outputRowData.push_back(_simCounter._counterCost[k]);
	SaveOutputRow_MapTemplate(_simCounter._counterCost_byCategory, k, _outputRowData, FLAG_SMOOTH_RESULTS);
	SaveOutputRow_MapTemplate(_simCounter._counterCost_byInsr, k, _outputRowData, FLAG_SMOOTH_RESULTS);
	SaveOutputRow_MapTemplate(_simCounter._counterCost_byLivDis, k, _outputRowData, FLAG_SMOOTH_RESULTS);

	_outputRowData.push_back(_simCounter._counter_transplant_DC[k] + _simCounter._counter_transplant_HCC[k]);
	_outputRowData.push_back(_simCounter._counterTxFailure_numRetxCandidates[k]);

	_outputRowData.push_back(_simCounter._counter_aware_F0F4[k]);
	SaveOutputRow_MapTemplate(_simCounter._counter_aware_byFib, k, _outputRowData);

	// ----------- print out -------------------------------
	_outf << START_YR + k << "\t";
	_outf << fixed;
	for (int i = 0; i < _outputRowData.size(); i++) {		
			_outf << _outputRowData[i] << "\t";
	}
	_outf << endl;
	_outputRowData_previousCycle = _outputRowData;
	if (_curYear == END_YR) { _outf.close(); }
	return 0;
}

int BurdenModelSim::OutputAgeDistr(const map<double,int> & ageDistr, string argFile)
{
	ofstream outf(argFile);
	if (outf.fail()) {
		ExitWithMsg("Can't open/create file: " + argFile);
	}
	for (map<double, int>::const_iterator it = ageDistr.begin();
		it != ageDistr.end(); it++) {
		outf << it->first << "\t" << it->second << endl;
	}
	outf.close();



	return 0;
}


int BurdenModelSim::ChangeOfCohort(const modelParamType & argModelParam)
{
	for (list<patientType*>::iterator ptIter = _listAllPts.begin(); ptIter != _listAllPts.end(); ptIter++) {
		// -------------- general (NHANES) population -------------
		if ((*ptIter)->_cohort == cohort_defaultSingleCohort) {

		}
		// -------------- homeless --------------------------------
		else if ((*ptIter)->_cohort == cohort_homeless) {
		}
		// -------------- in prison/jail --------------------------
		else if ((*ptIter)->_cohort == cohort_incarcerated) {
			//double pRelease = _curYear == (int)argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_year", cohort_hospitalized) + 1
			//	? argModelParam._nonNHANESPplData._incarcerated_release_prob + argModelParam._nonNHANESPplData._incacerated_jail_pct
			//	: argModelParam._nonNHANESPplData._incarcerated_release_prob;
			//if (_rnd_pt.GetU01() < pRelease) {
			//	// change cohort 
			//	(*ptIter)->_cohort = cohort_defaultSingleCohort;
			//	(*ptIter)->_tx_group = tx_group_general;
			//	(*ptIter)->_bcScr_group = bc_group_general;

			//	// change insurance status, but no change in awareness status
			//	if ((*ptIter)->currentAge >= AGE_CUTOFF_MEDICARE) {
			//		(*ptIter)->_insurStatus = insr_medicare;
			//	}// for younger patients				
			//	else {// for older patients
			//		// before 2014
			//		if (_curYear < ACA_START_YR) {
			//			(*ptIter)->_insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_insurance, _rnd_pt);
			//		}
			//		// since 2014 (ACA starts)
			//		else {
			//			(*ptIter)->_insurStatus = DiscreteDistrSampler_DistrTable(argModelParam._disBurdnData._table_distr_insurance_ACA_by_year.find(_curYear)->second, _rnd_pt);
			//		}
			//	}// end age cut-off
			//}
		}
		// --------------  hospitalized ---------------------------
		else if ((*ptIter)->_cohort == cohort_hospitalized) {
			// come back to NHANES population next year
			// no change in insurance status (due to change of cohort)
			if (_curYear == (int)argModelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_year", cohort_hospitalized) + 1) {
				(*ptIter)->_cohort = cohort_defaultSingleCohort;
			}
		}
		// --------------  nursing home ---------------------------
		else if ((*ptIter)->_cohort == cohort_nursingHome) {

		}
		// --------------  indian reservation----------------------
		else if ((*ptIter)->_cohort == cohort_indianReservation) {
		}
		else {
			ExitWithMsg("Error @ ChangeOfCohort(): Unknown cohort type = " + basicToStr((int)(*ptIter)->_cohort));
		}
	}
	return 0;
}
