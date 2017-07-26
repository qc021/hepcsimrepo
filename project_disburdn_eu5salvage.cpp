#include "project_disburdn_eu5salvage.h"
using namespace std;

Project_SalvageEU5::Project_SalvageEU5()
{
	_drugCostReduction = 0;



	_modelParam._transData.pr_F3_HCC = ENABLE_HCC_FROM_F3 ? 0.008 : 0;
	_modelParam._transData.pr_F3SVR_HCC = ENABLE_HCC_FROM_F3SVR ? 0.001926 : 0;
	_modelParam._transData.pr_regression = ENABLE_REGRESSION ? 0.4615 : 0;


}





//
//int Project_Acute::CEA_OneWay(int argCmpIdx, string argFileOnewayParam)
//{
//
//
//
//	ofstream outf_low, outf_high;
//	outf_low.open("output_project_acute/SA_LOW_" + _listCmp[argCmpIdx] + ".txt");
//	outf_high.open("output_project_acute/SA_HIGH_" + _listCmp[argCmpIdx] + ".txt");
//	outf_low << fixed << showpoint;
//	outf_high << fixed << showpoint;
//
//	vector<string> sa_varName;
//	vector<SParamTriplet> sa_varVal;
//	ReadOneWaySAParamValueRange(argFileOnewayParam, sa_varName, sa_varVal);
//
//	int sa_counter = 1;
//	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
//	/********************************** Define trial arms **************************************************************/
//	// === scenarios ===
//	int nScenarios = 2;
//	stateType listState[] = { s_Acute,s_Acute };
//	char listGender[] = { 'M','F' };
//	double listAge[2];
//	for (int k = 0; k < nScenarios; k++) { listAge[k] = _age_initial; }
//
//
//	for (int k = 0; k < nScenarios; k++) {
//		// define cohort
//		// create patient profile:
//		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
//		baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', _listGenotype[argCmpIdx], 'N', 'A');
//
//		// ===== baseline ======
//
//		_modelParam = _baselineModelParam;// recover the baseline values
//
//		vector<double> r1_b = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
//		vector<double> r2_b = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);
//
//		assert(r1_b.size() == r2_b.size());
//		outf_low << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
//			<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
//			<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
//			<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
//		outf_high << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
//			<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
//			<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
//			<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
//		for (int i = 2; i < r1_b.size(); i++) {
//			outf_low << r1_b[i] << "\t" << r2_b[i] << "\t";
//			outf_high << r1_b[i] << "\t" << r2_b[i] << "\t";
//		}
//		outf_low << endl;
//		outf_high << endl;
//
//		sa_counter++;
//
//		for (int j = 0; j < sa_varName.size(); j++) {
//
//			// ===== low value ======
//			cout << endl<< " - change value of paramter: " << sa_varName[j] << " = " << sa_varVal[j]._lb<<"\t";
//			_modelParam = _baselineModelParam;// recover the baseline values
//			ChangeValue(sa_varName[j], sa_varVal[j]._lb, _modelParam);
//			vector<double> r1_low = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
//			vector<double> r2_low = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);
//
//			outf_low << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
//				<< r1_low[0] << "\t" << r2_low[0] << "\t" << r2_low[0] - r1_low[0] << "\t"
//				<< r1_low[1] << "\t" << r2_low[1] << "\t" << r2_low[1] - r1_low[1] << "\t"
//				<< (r2_low[1] - r1_low[1]) / (r2_low[0] - r1_low[0]) << "\t";
//
//			for (int i = 2; i < r1_low.size(); i++) {
//				outf_low << r1_low[i] << "\t" << r2_low[i] << "\t";
//			}
//			outf_low << endl;
//
//			// ===== high value ======
//			cout << endl<< " - change value of paramter: " << sa_varName[j] << " = " << sa_varVal[j]._ub << "\t";
//			_modelParam = _baselineModelParam;// recover the baseline values
//			ChangeValue(sa_varName[j], sa_varVal[j]._ub, _modelParam);
//			vector<double> r1_high = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
//			vector<double> r2_high = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);
//
//			outf_high << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
//				<< r1_high[0] << "\t" << r2_high[0] << "\t" << r2_high[0] - r1_high[0] << "\t"
//				<< r1_high[1] << "\t" << r2_high[1] << "\t" << r2_high[1] - r1_high[1] << "\t"
//				<< (r2_high[1] - r1_high[1]) / (r2_high[0] - r1_high[0]) << "\t";
//			for (int i = 2; i < r1_high.size(); i++) {
//				outf_high << r1_high[i] << "\t" << r2_high[i] << "\t";
//			}
//			outf_high << endl;
//
//			sa_counter++;
//		}
//
//	}
//
//
//	_modelParam = _baselineModelParam;//restore the nominal values
//
//
//
//	outf_low.close(); outf_high.close();
//	return 0;
//}
//
//
//int Project_Acute::CEA_PSA(int argIdx, int argNBatches) {
//	ofstream outf_psa;
//	outf_psa.open("output_project_acute/PSA_" + basicToStr(argIdx) + ".txt");
//	outf_psa << fixed << showpoint;
//
//	psaDistrType psaDistr;
//	psaDistr.ReadDistrForPSA("project_acute_input_PSA_distr.txt");
//
//	HepCSim mySim_psaSampler;		//used for sampling
//	mySim_psaSampler.PSA_initializeSampler(psaDistr);
//
//	int nIterOfThisBatch = NUM_RUNS_PSA / argNBatches;
//	for (int n = 0; n < nIterOfThisBatch; n++) {
//		if ((n + 1) % (nIterOfThisBatch / 10) == 0) {
//			cout << (n + 1) << " PSA iterations finished...(" << (n + 1)*1.0 / (1.0 * nIterOfThisBatch) * 100 << "%)" << endl;
//		}
//		// sample
//		mySim_psaSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
//		// solve
//
//
//		//if (argIdx == 0 && n == 0) {
//		//	// print table head
//		//	outf_psa << "iter\t"
//		//		<< "QALY1\tCost1\t"
//		//		<< "QALY2\tCost2\t"
//		//		<< endl;
//		//}
//
//		vector<double> vecResult;
//		outf_psa << argIdx * nIterOfThisBatch + n << "\t";
//		vecResult = GetAggregatedResults();
//		for (int k = 0; k < vecResult.size(); k++) {
//			outf_psa << vecResult[k] << "\t";
//		}
//		outf_psa << endl;
//
//
//
//	}
//	outf_psa.close();
//
//	return 0;
//
//
//}



double Project_SalvageEU5::GetValue_byCountry(string argVarName, string argStrCountry)
{
	map<string, map<string, double>>::const_iterator it = _param_by_country.find(argVarName);
	if (it == _param_by_country.end()) {
		ExitWithMsg("Error @ Project_SalvageEU5::GetValue_byCountry: Can't find parameter - " + argVarName);
	}
	map<string, double > ::const_iterator it2 = it->second.find(argStrCountry);
	if (it2 == it->second.end()) {
		ExitWithMsg("Error @ Project_SalvageEU5::GetValue_byCountry: Can't find country code - " + argStrCountry);
	}
	return it2->second;
	
}

int Project_SalvageEU5::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
	if ("q_Acute" == varName) {
		argModelParam._qolData.q_acute = argVal;
	}
	else if ("q_F0" == varName) {
		argModelParam._qolData.q_F0 = argVal;
	}
	else if ("q_F1" == varName) {
		argModelParam._qolData.q_F1 = argVal;
	}
	else if ("q_F2" == varName) {
		argModelParam._qolData.q_F2 = argVal;
	}
	else if ("q_F3" == varName) {
		argModelParam._qolData.q_F3 = argVal;
	}
	else if ("q_CoCirr" == varName) {
		argModelParam._qolData.q_CoCirr = argVal;
	}
	else if ("q_DeCirr" == varName) {
		argModelParam._qolData.q_DeCirr = argVal;
	}
	else if ("q_HCC" == varName) {
		argModelParam._qolData.q_HCC = argVal;
	}
	else if ("q_LivTr" == varName) {
		argModelParam._qolData.q_LivTr = argVal;
	}
	else if ("q_PostLivT" == varName) {
		argModelParam._qolData.q_PostLivTr = argVal;
	}
	else if ("q_SVR" == varName) {
		argModelParam._qolData.q_SVR = argVal;
	}
	else if ("q_Anemia" == varName) {
		argModelParam._qolData.q_Dec_Anemia = argVal;
	}
	else if ("q_TX_oSOC" == varName) {
		argModelParam._qolData.q_TX_oSOC = argVal;
	}
	else if ("q_Tx_DAA" == varName) {
		argModelParam._qolData.q_TX_DAA = argVal;
	}
	else if ("c_Acute" == varName) {
		argModelParam._costData.c_acute = argVal;
	}
	else if ("c_F0" == varName) {
		argModelParam._costData.c_F0 = argVal;
	}
	else if ("c_F1" == varName) {
		argModelParam._costData.c_F1 = argVal;
	}
	else if ("c_F2" == varName) {
		argModelParam._costData.c_F2 = argVal;
	}
	else if ("c_F3" == varName) {
		argModelParam._costData.c_F3 = argVal;
	}
	else if ("c_CoCirr" == varName) {
		argModelParam._costData.c_CoCirr = argVal;
	}
	else if ("c_DeCirr" == varName) {
		argModelParam._costData.c_DeCirr = argVal;
	}
	else if ("c_DeCirr1yrPlus" == varName) {
		argModelParam._costData.c_DeCirr1yrPlus = argVal;
	}
	else if ("c_HCC" == varName) {
		argModelParam._costData.c_HCC = argVal;
	}
	else if ("c_LivTr" == varName) {
		argModelParam._costData.c_LivTr = argVal;
	}
	else if ("c_PostLivTr" == varName) {
		argModelParam._costData.c_PostLivTr = argVal;
	}
	else if ("pF0_F1_SA" == varName) {
		argModelParam._transData.pr_F0_F1 = argVal;
	}
	else if ("pF1_F2_SA" == varName) {
		argModelParam._transData.pr_F1_F2 = argVal;
	}
	else if ("pF2_F3_SA" == varName) {
		argModelParam._transData.pr_F2_F3 = argVal;
	}
	else if ("pF3_F4_SA" == varName) {
		argModelParam._transData.pr_F3_CoCirr = argVal;
	}
	else if ("pF4_DC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_DeCirr = argVal;
	}
	else if ("pF4_HCC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_HCC = argVal;
	}
	else if ("pDC_HCC_SA" == varName) {
		argModelParam._transData.pr_DeCirr_HCC = argVal;
	}
	else if ("pDC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_DeCirr_LivTr = argVal;
	}
	else if ("pMort_dc_cyc_1_SA" == varName) {
		argModelParam._transData.pr_DeCirr_DeathLiv = argVal;
	}
	else if ("pMort_dc_cyc_2_SA" == varName) {
		argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = argVal;
	}
	else if ("pHCC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_HCC_LivTr = argVal;
	}
	else if ("pMort_hcc_cyc_SA" == varName) {
		argModelParam._transData.pr_HCC_DeathLiv = argVal;
	}
	else if ("pMort_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr_DeathLiv = argVal;
	}
	else if ("pMort_Post_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = argVal;
	}
	else if ("pr_SVR_CoCirr_DeCirr" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_DeCirr = argVal;
	}
	else if ("pr_SVR_CoCirr_HCC" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_HCC = argVal;
	}
	else if ("pSVR_Delta_oSOC" == varName) {
		argModelParam._transData.pr_SVR_Delta_oSOC = argVal;
	}
	else if ("pSVR_Delta_DAA" == varName) {
		argModelParam._transData.pr_SVR_Delta_DAA = argVal;
	}
	else if ("c_testing_genotype" == varName) {
		argModelParam._costData.c_testing_genotype = argVal;
	}
	else if ("c_testing_RNA" == varName) {
		argModelParam._costData.c_testing_RNA = argVal;
	}
	else if ("pAcuteSelfClearance" == varName) {
		argModelParam._transData._pr_acute_spont_clearance = argVal;
	}
	else if ("pr_SVR" == varName) {
		//argModelParam._transData._pr_svr_for_SA = argVal;
		argModelParam._transData.pr_SVR_Delta_DAA = 1.0 -   argVal / argModelParam._transData._pr_svr_for_SA;
	}
	else if ("pr_SVR_acute" == varName) {
		argModelParam._transData._pr_svr_for_SA_acute = argVal;
	}
	else if ("pr_SVR_F0" == varName) {
		argModelParam._transData._pr_svr_for_SA_F0 = argVal;		
	}
	else {
		ExitWithMsg("[Error] Project_Acute::ChangeValue(string varName, double argVal): Unknown parameter " + varName);
	}
	return 0;
}





int Project_SalvageEU5::BaseCaseRun()
{
	BurdenModelSim myModel;
	myModel.Run(_modelParam);

	return 0;
}


int Project_SalvageEU5::Initialize_Model_Parameter(string argStrCountry)
{
	ReadTable_Param_by_Country("project_disburdn_eu5_input_paramByCountry.txt");
	_modelParam._disBurdnData._num_initial_population = GetValue_byCountry("numInitialPpl", argStrCountry);

	//ReadTable_BackgroundMortality("Input_mortality_male_female_burdn.in");
	//ReadTable_BackgroundMortality("Input_mortality_male_female.in");
	ReadTable_BackgroundMortality_5yAgeGroup_EU5("Input_mortality_male_female_eu5.in", argStrCountry);
	//ReadTable_BackgroundMortality("Input_mortality_male_female_2011_DE.in");

	//ReadTable_TreatmentCapacity("project_disburdn_txfailure_input_tx_cap.txt");
	ReadTable_TreatmentCapacity("project_disburdn_eu5_input_tx_cap.txt", argStrCountry);
	//ReadTable_TreatmentCapacity("project_disburdn_eu5_input_tx_cap_noTxBefore2004.txt", argStrCountry);


	//ReadTable_NS5AMarketShare("project_disburdn_txfailure_input_ns5a_marketshare.txt");
	ReadTable_NS5AMarketShare("project_disburdn_eu5_input_ns5a_marketshare_"+argStrCountry+".txt");
	ReadTable_ExtendedUsePRandPI("project_disburdn_eu5_input_extended_use_pr_pi.txt", argStrCountry);
	
	//ReadTable_Incidence_Num("project_disburdn_txfailure_input_incidence.txt");
	ReadTable_Incidence_Num("project_disburdn_eu5_input_incidence.txt", argStrCountry);

	//ReadTable_Incidence_AgeDistr("project_disburdn_txfailure_input_inci_age_distr.txt");
	ReadTable_Incidence_AgeDistr("project_disburdn_eu5_input_inci_age_distr_"+argStrCountry+".txt");

	// [QC comment] override the probability of becoming aware from usual care:
	ReadTable_Prob_Aware("project_disburdn_txfailure_input_probAware.txt");
	for (map<typeStateDisBurdnModel, map<bool, map<double, double>>>::iterator it = _modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.begin();
		it != _modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.end(); it++) {
		for (map<bool, map<double, double>> ::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			for (map<double, double>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
				it3->second = GetValue_byCountry("rateScreening", argStrCountry);
			}
		}
	}
	

	if (RECOVER_SALVAGE_POSTER_VERSION) {
		Initialize_SVR_Table("project_disburdn_txfailure_input_svr_verPoster.txt", _modelParam._disBurdnData._table_SVR);
	}
	else {
		Initialize_SVR_Table("project_disburdn_txfailure_input_svr.txt", _modelParam._disBurdnData._table_SVR);
	}


	// ---- update transplant related parameters ------------
	ReadTable_TransplantCapacity("project_disburdn_eu5_input_transplant_cap.txt", argStrCountry);
	_modelParam._transData.pr_DeCirr_LivTr = 1.0;
	_modelParam._transData.pr_HCC_LivTr = 1.0;
	_modelParam._disBurdnData._flag_use_constant_transplant_capacity = false;


	// ---- update disBurdnDataType variables for EU5 analysis ------------
	_modelParam._disBurdnData._flag_extended_use_PR_PI = true;

	_modelParam._disBurdnData._prob_male = GetValue_byCountry("distrMale", argStrCountry); // 0.6422;
	_modelParam._disBurdnData._prob_receivable = 1; // 0.915;
	//_prob_TrR = 1.0; // unchanged
	
	_modelParam._disBurdnData._prob_chronic = 0.78; // unchanged
	//_modelParam._disBurdnData._prob_chronic = GetValue_byCountry("viremicRate", argStrCountry);

	_modelParam._disBurdnData._prob_tx_coverage_by_medicarePartD = 1.0; // 0.9;
	_modelParam._disBurdnData._prob_naive_initialization = 1.0; // unchanged // 0.65, or 0.17???
	_modelParam._disBurdnData._prob_chronic_contra = 0.346; // unchanged
	_modelParam._disBurdnData._prob_mod_given_contra = 0.231 / 0.346; // unchanged
	_modelParam._disBurdnData._prob_screen_accept = 1.0; // 0.819;		//.91 * .90, Rein 2012. (73.71-90.09)
	//									//percentage of people who accept to be screened for HepC
	//									//Rein 2012 assumed that 28%of the chronic cases are aware of their infection, so screening would not benefit them.
	//_prob_birthCohortScr_insured = 0.9; // unchanged
	//_prob_birthCohortScr_uninsured = 0.1; // unchanged
	//_prob_delay_tx_until_wave1_f0f2 = 0.75; // unchanged
	//_prob_delay_tx_until_wave1_f3 = 0.25; // unchanged


	// Assumption: use the same awareness rate for [new incidences] in US	
	double p_aware_in_new_incidence = 0.07; // 7% CDC data // diagnosed rate = 57% [Razavi2014]
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[40] = p_aware_in_new_incidence;// 0.0581;// 0.42 * 0.1383; // lower_bounds(): < 40
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[50] = p_aware_in_new_incidence;// 0.1950;// 0.71*0.2747;	// 40-49
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[60] = p_aware_in_new_incidence;// 0.1629;// 0.61*0.2670;	// 50-59
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[200] = p_aware_in_new_incidence;// 0.0577;// 0.27*0.2136; // 60-100
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[40] = p_aware_in_new_incidence;// 0.1763;// 0.53*0.3326;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[50] = p_aware_in_new_incidence;// 0.6045;// 0.915*0.6606;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[60] = p_aware_in_new_incidence;// 0.5459;// 0.85*0.6423;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[200] = p_aware_in_new_incidence;// 0.3494;// 0.68*0.5138;

	// Read the (calibrated) initial awareness rate from the input file
	//_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured = _modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured;
	//_modelParam._disBurdnData._table_aware_prob_initial_population_insured = _modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured;
	_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured[40] = GetValue_byCountry("awarenessRateInitialPop_40", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured[50] = GetValue_byCountry("awarenessRateInitialPop_50", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured[60] = GetValue_byCountry("awarenessRateInitialPop_60", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured[200] = GetValue_byCountry("awarenessRateInitialPop_200", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_insured[40] = GetValue_byCountry("awarenessRateInitialPop_40", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_insured[50] = GetValue_byCountry("awarenessRateInitialPop_50", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_insured[60] = GetValue_byCountry("awarenessRateInitialPop_60", argStrCountry);
	_modelParam._disBurdnData._table_aware_prob_initial_population_insured[200] = GetValue_byCountry("awarenessRateInitialPop_200", argStrCountry);


	// [QC comment, 03/07/2017]
	// Read country specific genotype 
	// <genotype-int, probability>
	_modelParam._disBurdnData._table_distr_genotype[1] = GetValue_byCountry("distrG1", argStrCountry); // 0.73;
	_modelParam._disBurdnData._table_distr_genotype[2] = GetValue_byCountry("distrG2", argStrCountry); //0.14;
	_modelParam._disBurdnData._table_distr_genotype[3] = GetValue_byCountry("distrG3", argStrCountry); //0.08;
	_modelParam._disBurdnData._table_distr_genotype[456] = GetValue_byCountry("distrG456", argStrCountry); //0.05;



	// <fibrosis, probability>	
	_modelParam._disBurdnData._table_distr_fib[state_CoCirr] = GetValue_byCountry("distrCirr", argStrCountry); //0; // 0.0554;
	_modelParam._disBurdnData._table_distr_fib[state_DeCirr] = 0; //0.00206;
	_modelParam._disBurdnData._table_distr_fib[state_DeCirr1yrPlus] = 0; //0.00928;
	_modelParam._disBurdnData._table_distr_fib[state_HCC] = GetValue_byCountry("distrHCC", argStrCountry); //0; //0;
	_modelParam._disBurdnData._table_distr_fib[state_F0] = GetValue_byCountry("distrF0", argStrCountry); //1.0; // 0.44494;
	_modelParam._disBurdnData._table_distr_fib[state_F1] = GetValue_byCountry("distrF1", argStrCountry); // 0; //0.26474;
	_modelParam._disBurdnData._table_distr_fib[state_F2] = GetValue_byCountry("distrF2", argStrCountry); //0; //0.1357;
	_modelParam._disBurdnData._table_distr_fib[state_F3] = GetValue_byCountry("distrF3", argStrCountry); //0; //1 - 0.91212;



	// ------------------------------------------------
	// [QC comment] 
	// There will be no "insurance type" in EU5 analysis
	// "Pretend" all patients having "private" before age 65, and medicare after age 65
	// but no difference will be made by the insurance status.
	// ------------------------------------------------
	//_table_distr_insurance[insr_uninsured] = 0.267;
	//_table_distr_insurance[insr_private] = 0.498; // 0.7649 - 0.2666;
	//_table_distr_insurance[insr_medicaid] = 0.143; // 0.9078 - 0.7649;
	//_table_distr_insurance[insr_military] = 0.092; // 1.0 - 0.9078;
	_modelParam._disBurdnData._table_distr_insurance[insr_private] = 1.0; 
	_modelParam._disBurdnData._table_distr_insurance[insr_uninsured] = 0.0;
	_modelParam._disBurdnData._table_distr_insurance[insr_medicaid] = 0.0;
	_modelParam._disBurdnData._table_distr_insurance[insr_military] = 0.0;

	//// [QC comment] No ACA changes in EU5 analysis
	// Assume change is 0
	_modelParam._disBurdnData.ACATotalratio[2100] = 0.0;
	_modelParam._disBurdnData.ACAPrivateratio[2100] = 0.0;
	_modelParam._disBurdnData.ACAMedicaidratio[2100] = 0.0;
	_modelParam._disBurdnData.ACAPrivatechange[2100] = 0.0;
	_modelParam._disBurdnData.ACAMedicaidchange[2100] = 0.0;
	
	double p_uninsured = _modelParam._disBurdnData._table_distr_insurance[insr_uninsured];
	double p_private = _modelParam._disBurdnData._table_distr_insurance[insr_private];
	double p_medicaid = _modelParam._disBurdnData._table_distr_insurance[insr_medicaid];
	double p_mil = _modelParam._disBurdnData._table_distr_insurance[insr_military];
	for (int y = ACA_START_YR; y <= END_YR; y++) {
		map<typeInsurance, double> insr_distr;
		insr_distr[insr_uninsured] = p_uninsured * (1 - _modelParam._disBurdnData.ACAPrivateratio.lower_bound(y)->second - _modelParam._disBurdnData.ACAMedicaidratio.lower_bound(y)->second);
		insr_distr[insr_private] = p_private + p_uninsured *  _modelParam._disBurdnData.ACAPrivateratio.lower_bound(y)->second;
		insr_distr[insr_medicaid] = p_medicaid + p_uninsured * _modelParam._disBurdnData.ACAMedicaidratio.lower_bound(y)->second;
		insr_distr[insr_military] = p_mil;

		double s = insr_distr[insr_uninsured] + insr_distr[insr_private] + insr_distr[insr_medicaid] + insr_distr[insr_military];
		if (abs(s - 1.0) > EPSILON) {
			ExitWithMsg("ERROR @ disBurdnDataType(): adjusting insurance distribution given ACA: distributio does not sum to 1 (sum = " + basicToStr(s) + ")");
		}
		_modelParam._disBurdnData._table_distr_insurance_ACA_by_year.insert(pair<int, map<typeInsurance, double>>(y, insr_distr));
	}


	// ------------ [QC comment] Following remain unchanged ----------------------------------
	//_table_respStatus_distr_naive[trStatus_contraInd_mod] = _prob_mod_given_contra * _prob_chronic_contra;
	//_table_respStatus_distr_naive[trStatus_contraInd_nonmod] = (1 - _prob_mod_given_contra) * _prob_chronic_contra;
	//_table_respStatus_distr_naive[trStatus_naive] = 1 - _prob_chronic_contra;

	//_table_respStatus_distr_expr_G1[trStatus_relap] = 0.53;
	//_table_respStatus_distr_expr_G1[trStatus_partial] = 0.19;
	//_table_respStatus_distr_expr_G1[trStatus_null] = 1 - 0.53 - 0.19;

	//_table_respStatus_distr_expr_G234[trStatus_relap] = 0.47;
	//_table_respStatus_distr_expr_G234[trStatus_partial] = 0.16;
	//_table_respStatus_distr_expr_G234[trStatus_null] = 1 - 0.47 - 0.16;


	// -------------------------------------------------
	// [QC comment] 
	// age categories and distribution are read from table directly
	// See Project_SalvageEU5::ReadTable_Incidence_AgeDistr()
	// -------------------------------------------------
	//_vec_age_category_forInitialization.push_back(18);
	//_vec_age_category_forInitialization.push_back(24);
	//_vec_age_category_forInitialization.push_back(34);
	//_vec_age_category_forInitialization.push_back(44);
	//_vec_age_category_forInitialization.push_back(54);
	//_vec_age_category_forInitialization.push_back(64);

	//_vec_age_category_forNewIncidence.push_back(18);
	//_vec_age_category_forNewIncidence.push_back(20);
	//_vec_age_category_forNewIncidence.push_back(30);
	//_vec_age_category_forNewIncidence.push_back(40);
	//_vec_age_category_forNewIncidence.push_back(50);
	//_vec_age_category_forNewIncidence.push_back(60);
	//_vec_age_category_forNewIncidence.push_back(70);



	// [QC comment] Below will not be used. On the safe side, set all of them = 0
	_modelParam._disBurdnData._table_universal_scr_rate[-1] = 0;//0;
	_modelParam._disBurdnData._table_universal_scr_rate[0] = 0;// 0.142857;
	_modelParam._disBurdnData._table_universal_scr_rate[1] = 0;//0.16667;
	_modelParam._disBurdnData._table_universal_scr_rate[2] = 0;//0.2;
	_modelParam._disBurdnData._table_universal_scr_rate[3] = 0;//0.25;
	_modelParam._disBurdnData._table_universal_scr_rate[4] = 0;//0.333;
	_modelParam._disBurdnData._table_universal_scr_rate[5] = 0;//0.5;
	_modelParam._disBurdnData._table_universal_scr_rate[6] = 0;//1;
	_modelParam._disBurdnData._table_universal_scr_rate[7] = 0;//0;
	_modelParam._disBurdnData._table_universal_scr_rate[100] = 0;//0;

	return 0;
}


int Project_SalvageEU5::ReadTable_Param_by_Country(string argStr)
{
	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	// read table header, register countries
	vector<string> listCountry;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word;
	while (iss >> word) {
		listCountry.push_back(word);
	}


	string str;
	while (inFile >> str) {	
		if (str == "//") {
			getline(inFile, str);
			continue;
		}
		map<string, double> row;
		for (int k = 0; k < (int)listCountry.size(); k++) {
			inFile >> row[listCountry[k]];
		}
		_param_by_country[str] = row;
	}
	inFile.close();
	return 0;


}
int Project_SalvageEU5::ReadTable_NS5AMarketShare(string argStr)
{
	// Read the NS5A/nonNS5A market share input data
	ifstream inFileNS5AMarketShare;
	inFileNS5AMarketShare.open(argStr);
	if (inFileNS5AMarketShare.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	while(!inFileNS5AMarketShare.eof()){
		int year;
		map<int,double> share;
		inFileNS5AMarketShare>>year;
		double s1,s2,s3,s4;
		inFileNS5AMarketShare>>s1>>s2>>s3>>s4;
		share[1]=s1;
		share[2]=s2;
		share[3]=s3;
		share[456]=s4;
		_modelParam._disBurdnData._table_ns5a_marketshare[year]=share;
	}
	//cout<<"testing: _tfa_market_share = " << _modelParam._disBurdnData._table_ns5a_marketshare[2015][1]<<"\t"<<_modelParam._disBurdnData._table_ns5a_marketshare[2016][3]<<endl;
	inFileNS5AMarketShare.close();
	return 0;
}

int Project_SalvageEU5::ReadTable_ExtendedUsePRandPI(string argStr, string argCountry)
{
	ifstream inFile;
	inFile.open(argStr.c_str());		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	// read table header, register countries
	vector<string> listCountry;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word, line;

	map<string, map<int, map<int, double>>> table_pr_pi, table_pr; // [country][year][genotype]
	
	map<int, map<int, double>> map1; 
	
	while (iss >> word) {
		listCountry.push_back(word);
		table_pr_pi[word] = map1;
		table_pr[word] = map1;
	}

	vector<int> listGt;
	listGt.push_back(1);
	listGt.push_back(2);
	listGt.push_back(3);
	listGt.push_back(456);

	vector<double> ageDistr;
	while (inFile >> word) {
		if (word == "//") {
			getline(inFile, line);
		}
		else {
			int gt, year;
			year = atoi(word.c_str());			
			
			for (int k = 0; k < (int)listCountry.size(); k++) {
				for(int j=0; j<(int)listGt.size(); j++){
				double v;
				inFile >> v;
				table_pr_pi[listCountry[k]][year][listGt[j]] = v;
				inFile >> v;
				table_pr[listCountry[k]][year][listGt[j]] = v;
				}
			}
		}
	}
	
	_modelParam._disBurdnData._table_pr_extended_use_pr_pi = table_pr_pi[argCountry];	
	_modelParam._disBurdnData._table_pr_extended_use_pr = table_pr[argCountry];
	inFile.close();

	return 0;

}

int Project_SalvageEU5::ReadTable_Incidence_Num(string argStr, string argCountry)
{
	//ReadTable_2Col(argStr, _modelParam._disBurdnData._table_incidence);	
	map<string, map<int, long>> tableIncidence;

	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	// read table header, register countries
	vector<string> listCountry;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word;
	map<int, long> emptyTable;
	while (iss >> word) {
		listCountry.push_back(word);
		tableIncidence[word] = emptyTable;
	}

	int year, val;
	while (inFile >> year) {		
		for (int k = 0; k < (int)listCountry.size(); k++) {
			inFile >> val;
			tableIncidence.find(listCountry[k])->second.insert(pair<int, long>(year, val));
		}		
	}
	inFile.close();
	
	_modelParam._disBurdnData._table_incidence = tableIncidence[argCountry];
	return 0;
}


int Project_SalvageEU5::ReadTable_TreatmentCapacity(string argStr, string argCountry)
{
	//ReadTable_2Col(argStr, _modelParam._disBurdnData._table_tx_capacity);
	map<string, map<int, int>> tableCap;

	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}
								// read table header, register countries
	vector<string> listCountry;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word;
	map<int, int> emptyTable;
	while (iss >> word) {
		listCountry.push_back(word);
		tableCap[word] = emptyTable;
	}

	int year, val;
	while (inFile >> year) {
		for (int k = 0; k < (int)listCountry.size(); k++) {
			inFile >> val;
			tableCap.find(listCountry[k])->second.insert(pair<int, int>(year, val));
		}
	}
	inFile.close();

	_modelParam._disBurdnData._table_tx_capacity = tableCap[argCountry];
	return 0;
}


int Project_SalvageEU5::ReadTable_TransplantCapacity(string argStr, string argCountry)
{
	//ReadTable_2Col(argStr, _modelParam._disBurdnData._table_tx_capacity);
	map<string, map<int, int>> tableCap;

	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}
	// read table header, register countries
	vector<string> listCountry;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word;
	map<int, int> emptyTable;
	while (iss >> word) {
		listCountry.push_back(word);
		tableCap[word] = emptyTable;
	}

	int year, val;
	while (inFile >> year) {
		for (int k = 0; k < (int)listCountry.size(); k++) {
			inFile >> val;
			tableCap.find(listCountry[k])->second.insert(pair<int, int>(year, val));
		}
	}
	inFile.close();

	_modelParam._disBurdnData._table_transplant_capacity = tableCap[argCountry];
	return 0;
}

int Project_SalvageEU5::ReadTable_Incidence_AgeDistr(string argStr)
{


	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.clear();
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(18);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(20);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(30);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(40);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(50);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(60);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(70);

	//ReadTable("project_disburdn_txfailure_input_inci_age_distr.txt", _modelParam._disBurdnData._table_inci_age_distr);

	//// reset the default age categories
	//_modelParam._disBurdnData._vec_age_category_forInitialization.clear();
	//_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(18);
	//_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(24);
	//_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(34);
	//_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(64);
	//
	//vector<double> ad; 	
	//ad.push_back(0.23);
	//ad.push_back(0.07);
	//ad.push_back(0);
	//ad.push_back(0.7);
	//_modelParam._disBurdnData._table_inci_age_distr[0] = ad;


	// reset the default age categories
	_modelParam._disBurdnData._vec_age_category_forInitialization.clear();
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.clear();

	//ReadTable(argStr, _modelParam._disBurdnData._table_inci_age_distr);
	ifstream inFile;
	inFile.open(argStr.c_str());		
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	string word, line;
	double v;
	vector<double> ageDistr;
	while (inFile >> word) {
		if (word == "//") {
			getline(inFile, line);
		}
		else if (word == "INITIAL") {
			ageDistr.clear();
			while(inFile >> word) {
				if (word == "END") {
					_modelParam._disBurdnData._table_inci_age_distr[0] = ageDistr;
					break;
				}
				_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(atoi(word.c_str()));								
				inFile >> v;
				ageDistr.push_back(v);
			}
		}
		else if (word == "NEW_INCIDENCE") {
			ageDistr.clear();
			while (inFile >> word) {
				if (word == "END") {
					_modelParam._disBurdnData._table_inci_age_distr[2100] = ageDistr;
					break;
				}
				_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(atoi(word.c_str()));
				inFile >> v;
				ageDistr.push_back(v);
			}
		}
		else {
			ExitWithMsg("ERROR @ Readtable_Incidence_AgeDistr(): Error when reading the input file: " + argStr);
		}
	}

	inFile.close();

	return 0;
}

int Project_SalvageEU5::ReadTable_Prob_Aware(string argStr)
{
	int vfib;
	bool vInsured;
	double vAge;
	double vProb;
	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	string str;
	while (inFile >> str) {
		if ("//" == str) { getline(inFile, str); continue; }
		int tmp = atoi(str.c_str());
		for (int f = 0; f < 5; f++) {
			map<bool, map<double, double>> mp2;
			for (int b = 1; b >=0; b--) {
				map<double, double> mp;
				for (int idxAgeCat = 0; idxAgeCat < 4; idxAgeCat++) {
					vfib = tmp;
					inFile >> vInsured;
					inFile >> vAge;
					inFile >> vProb;

					assert(b == (int)vInsured && f == vfib);
					mp[vAge] = vProb;
					inFile >> tmp;
				}
				mp2[(bool)b] = mp;

			}
			_modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age[(typeStateDisBurdnModel)f] = mp2;
		}
		
	}
	inFile.close();

	return 0;
}





int Project_SalvageEU5::ReadTable_BackgroundMortality(string argStr)
{
	ifstream inFile2;
	inFile2.open(argStr);		//US life tables of 2007, published in 2011
	int k;
	while (inFile2 >> k) {
		//cout << "k= " << k << endl;
		inFile2 >> _modelParam._disBurdnData._bgMort_male[k]
			>> _modelParam._disBurdnData._bgMort_female[k]
			>> _modelParam._disBurdnData._qol_male[k]
			>> _modelParam._disBurdnData._qol_female[k];
	}
	inFile2.close();
	return 0;
}
int Project_SalvageEU5::ReadTable_BackgroundMortality_5yAgeGroup_EU5(string argStr, string argCountry)
{

	// ======================================
	// NOTE: life table is for 5-year age groups
	// need to use .lowerbound() function, 
	// instead of .find() function to retrive the background mortality
	// ======================================
	map<string, map<int, double>> table_mFemale, table_mMale, table_qFemale, table_qMale;
		
	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	// read table header, register countries
	vector<string> listCountry;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word;
	map<int, double> emptyTable;
	while (iss >> word) {
		listCountry.push_back(word);
		table_mFemale[word] = emptyTable;
		table_mMale[word] = emptyTable;
		table_qFemale[word] = emptyTable;
		table_qMale[word] = emptyTable;
	}

	string str;
	double val;
	while (inFile >> str) {
		if (str == "//") {
			getline(inFile, str);
			continue;
		}

		for (int k = 0; k < (int)listCountry.size(); k++) {
			inFile >> val;
			table_mFemale.find(listCountry[k])->second.insert(pair<int, double>(atoi(str.c_str()), val));
			inFile >> val;
			table_mMale.find(listCountry[k])->second.insert(pair<int, double>(atoi(str.c_str()), val));
			inFile >> val;
			table_qFemale.find(listCountry[k])->second.insert(pair<int, double>(atoi(str.c_str()), val));
			inFile >> val;
			table_qMale.find(listCountry[k])->second.insert(pair<int, double>(atoi(str.c_str()), val));
		}
	}
	inFile.close();


	// fill-in n the table
	_modelParam._disBurdnData._bgMort_male.clear();
	_modelParam._disBurdnData._bgMort_female.clear();
	_modelParam._disBurdnData._qol_male.clear();
	_modelParam._disBurdnData._qol_female.clear();

	for (int age = 0; age < 120; age++) {
		double p_male = table_mMale[argCountry].lower_bound(age)->second;
		double p_female = table_mFemale[argCountry].lower_bound(age)->second;

		if (age == 0) {
			// do  nothing
		}
		else if (age >= 1 && age <= 4) {
			p_male = 1.0 - pow(1.0 - p_male, 1.0 / 4.0);
			p_female = 1.0 - pow(1.0 - p_female, 1.0 / 4.0);
		}
		else {
			p_male = 1.0 - pow(1.0 - p_male, 1.0 / 5.0);
			p_female = 1.0 - pow(1.0 - p_female, 1.0 / 5.0);
		}

		_modelParam._disBurdnData._bgMort_male[age] = p_male;
		_modelParam._disBurdnData._bgMort_female[age] = p_female;
		_modelParam._disBurdnData._qol_male[age] = table_qMale[argCountry].lower_bound(age)->second;
		_modelParam._disBurdnData._qol_female[age] = table_qFemale[argCountry].lower_bound(age)->second;
	}


	//_modelParam._disBurdnData._bgMort_male = table_mMale[argCountry];
	//_modelParam._disBurdnData._bgMort_female = table_mFemale[argCountry];
	//_modelParam._disBurdnData._qol_male = table_qMale[argCountry];
	//_modelParam._disBurdnData._qol_female = table_qFemale[argCountry];
	return 0;


}



int Project_SalvageEU5::Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap)
{
	// read the table from file
	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	string strg;
	while (inFile >> strg) {
		typeTxResponseStatus resp;
		if ("treatment_naive" == strg) {
			resp = trStatus_naive;
		}
		else if ("contraind_mod" == strg) {
			resp = trStatus_contraInd_mod;
		}
		else if ("contraind_nonmod" == strg) {
			resp = trStatus_contraInd_nonmod;
		}
		else if ("relapse" == strg) {
			resp = trStatus_relap;
		}
		else if ("partial" == strg) {
			resp = trStatus_partial;
		}
		else if ("null" == strg) {
			resp = trStatus_null;
		}
		else if ("failed_firstgenpi" == strg) {
			resp = trStatus_failed_PI1;
		}
		else if ("failed_daa1_nonns5a" == strg) {
			resp = trStatus_failed_DAA1_nonNS5A;
		}
		else if ("failed_daa2_nonns5a" == strg) {
			resp = trStatus_failed_DAA2_nonNS5A;
		}
		else if ("failed_daa2_ns5a" == strg) {
			resp = trStatus_failed_DAA2_NS5A;
		}
		else if ("failed_daa3_nonns5a" == strg) {
			resp = trStatus_failed_DAA3_nonNS5A;
		}
		else if ("failed_daa3_ns5a" == strg) {
			resp = trStatus_failed_DAA3_NS5A;
		}
		else {
			ExitWithMsg("Error @ reading SVR table: unknown treatment response type = " + strg);
		}

		int gt; inFile >> gt;
		string strFib; inFile >> strFib;

		//cout << strg << "\t" << gt << "\t" << strFib << endl;
		string val;
		// read values for all treatment categories 
		// typedef map<int, map<typeTxCategory,map<typeTxResponseStatus,map<typeStateDisBurdnModel, double>>>> map_svr_table;

		if ("0-2" == strFib) {
			vector<double> vecval;
			for (int tx = 0; tx <= (int)txCat_DAA3_NS5A; tx++) { string str; inFile >> str; vecval.push_back("NA" == str ? -100 : atof(str.c_str())); }
			for (int k = 0; k <= (int)state_F2; k++) {

				argSVRMap[gt][txCat_PEGRBV][resp][(typeStateDisBurdnModel)k] = vecval[0];
				argSVRMap[gt][txCat_PI1][resp][(typeStateDisBurdnModel)k] = vecval[1];
				argSVRMap[gt][txCat_DAA1_nonNS5A][resp][(typeStateDisBurdnModel)k] = vecval[2];
				argSVRMap[gt][txCat_DAA2_nonNS5A][resp][(typeStateDisBurdnModel)k] = vecval[3];
				argSVRMap[gt][txCat_DAA2_NS5A][resp][(typeStateDisBurdnModel)k] = vecval[4];
				argSVRMap[gt][txCat_DAA3_nonNS5A][resp][(typeStateDisBurdnModel)k] = vecval[5];
				argSVRMap[gt][txCat_DAA3_NS5A][resp][(typeStateDisBurdnModel)k] = vecval[6];
			}
		}
		else {
			inFile >> val;  argSVRMap[gt][txCat_PEGRBV][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());
			inFile >> val;  argSVRMap[gt][txCat_PI1][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());
			inFile >> val;  argSVRMap[gt][txCat_DAA1_nonNS5A][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());
			inFile >> val;  argSVRMap[gt][txCat_DAA2_nonNS5A][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());
			inFile >> val;  argSVRMap[gt][txCat_DAA2_NS5A][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());
			inFile >> val;  argSVRMap[gt][txCat_DAA3_nonNS5A][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());
			inFile >> val;  argSVRMap[gt][txCat_DAA3_NS5A][resp][(typeStateDisBurdnModel)atoi(strFib.c_str())] = "NA" == val ? -100 : atof(val.c_str());

		}
	}
	inFile.close();
	


	// ******************************   INTERNAL CHECK POINTS *******************************
	// check-1: no PEGRBV and PI for contra_mod-G1, and all contra_nonmod patients


	// check-2: DAA1_nonNS5A: contra-mod = contra-nomod, across all genotypes

	// check-3: DAA1_nonNS5A: G2, G3, G456: identical SVR for relapse/partial/null responses 


	// check-4: DAA1_nonNS5A: treatment naive: identical SVR for G1, G2, G456

	// check-5: DAA1_nonNS5A, DAA2_nonNS5A, DAA3_nonNS5A are identical

	// check-6: no data for PI-failed & G2-6

	// check-7: no FirstGenPI for G2-6

	// check-8: [assumption] DAA2_NS5A: identical across naive/contraind_mod/contraind_nonmod/relapse/partial/null


	// check-9: DAA2_NS5A = DAA3_NS5A

	return 0;
}
