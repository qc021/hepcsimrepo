#include "project_disburdn_statex.h"
using namespace std;

Project_DisBurd_StateX::Project_DisBurd_StateX()
{
	_drugCostReduction = 0;



	_modelParam._transData.pr_F3_HCC = ENABLE_HCC_FROM_F3 ? 0.008 : 0;
	_modelParam._transData.pr_F3SVR_HCC = ENABLE_HCC_FROM_F3SVR ? 0.001926 : 0;
	_modelParam._transData.pr_regression = ENABLE_REGRESSION ? 0.4615 : 0;



	////// Overwrite DALY data ////////////////////////////////////////
	_modelParam._dalyData._discountRate = 0.0;
	_modelParam._dalyData._beta = 0.04;
	_modelParam._dalyData._C = 0.1658;
	_modelParam._dalyData._K = 0;

	//_dalyData.ReadLifeExpectancyData("./Input_life_expectancy_India.txt");
	_modelParam._dalyData.ReadLifeExpectancyData("./Input_life_expectancy_WHO.txt");
	// Egypt study: Estes et al 205 APT
	_modelParam._dalyData._dw_f0 = 0;
	_modelParam._dalyData._dw_f1 = 0;
	_modelParam._dalyData._dw_f2 = 0;
	_modelParam._dalyData._dw_f3 = 0;
	_modelParam._dalyData._dw_CoCirr = 0;
	_modelParam._dalyData._dw_DeCirr = 0.194;
	_modelParam._dalyData._dw_DeCirr1yrPlus = 0.194;
	_modelParam._dalyData._dw_HCC = 0.508;
	_modelParam._dalyData._dw_LivTr = 0;
	_modelParam._dalyData._dw_LivTr1yrPlus = 0;
	_modelParam._dalyData._dw_SVR = 0;


	// ------------ additional StateX parameters ------------------
	ReadTable_Param();
	_paramPK_tx_cap_statusquo = GetValue("txCapStatusQuo");
	_paramPK_tx_cap_future = GetValue("txCapFrom2017");
	_paramPK_universal_screeningRate = GetValue("screeningRate");
	_paramPK_universal_screeningCap = GetValue("screeningCap");


}


int Project_DisBurd_StateX::Initialize_StateX_Parameter(type_statex_incidence argInciMode)
{
	_modelParam._disBurdnData._num_initial_population = GetValue("InitialPplNum"); // initial population at year 1995
	_paramPK_start_year = GetValue("InitialPplYear");
	_paramPK_inci_const = GetValue("incidenceConst");
	// ****** Incidence of new HCV cases ************************************************
	_paramPK_inci_mode = argInciMode;
	if (_paramPK_inci_mode == STATEX_INCI_PROPORTIONAL) {
		_modelParam._disBurdnData._flag_inci_proportional_to_prevalence = true;
		for (int y = START_YR; y <= END_YR; y++) { _modelParam._disBurdnData._table_incidence_proportion[y] = _paramPK_inci_proportion; }

	}
	else if (_paramPK_inci_mode == STATEX_INCI_CONSTANT) {
		for (int y = START_YR; y <= END_YR; y++) { _modelParam._disBurdnData._table_incidence_proportion[y] = _paramPK_inci_proportion; }
		_modelParam._disBurdnData._year_start_constant_new_incidence = 2010;
		for (int y = _modelParam._disBurdnData._year_start_constant_new_incidence; y <= END_YR; y++) {
			_modelParam._disBurdnData._table_incidence[y] = _paramPK_inci_const;
		}
	}
	else if (_paramPK_inci_mode == STATEX_INCI_INCREASING) {
		for (int y = START_YR; y <= END_YR; y++) { _modelParam._disBurdnData._table_incidence_proportion[y] = _paramPK_inci_proportion; }
		_modelParam._disBurdnData._year_start_constant_new_incidence = 2010;
		for (int y = _modelParam._disBurdnData._year_start_constant_new_incidence; y <= END_YR; y++) {
			_modelParam._disBurdnData._table_incidence[y] = (int)(_paramPK_inci_const * (1.0 + 0.02 * (y - _modelParam._disBurdnData._year_start_constant_new_incidence)));
		}

	}
	else if (_paramPK_inci_mode == STATEX_INCI_TABLE) {
		ReadTable_2Col("project_disburdn_statex_input_incidence.txt", _modelParam._disBurdnData._table_incidence);

	}
	else {
		ExitWithMsg("Error @ Initialize_StateX_parameter(): Unknown incidence mode parameter = " + basicToStr(_paramPK_inci_mode));
	}



	// ****** background mortality ************************************************
	ReadTable_BackgroundMortality("Input_mortality_male_female_burdn.in");




	// ****** NS5A market share: simplified with 0 ****************************************************
	ReadTable_NS5AMarketShare("project_disburdn_txfailure_input_ns5a_marketshare.txt");
	//ReadTable_NS5AMarketShare("project_disburdn_eu5_input_ns5a_marketshare_"+argStrCountry+".txt");
	//ReadTable_ExtendedUsePRandPI("project_disburdn_eu5_input_extended_use_pr_pi.txt", argStrCountry);


	// ****** Treatment timeline ***************************************************************
	// PEG-RBV: 1998 - 2014
	// DAA:     2015 - end
	// See chanages in the section of "#elif defined(__PAKISTAN_SETTING__)" in disburdn_config.h


	// ****** Age distribution of new HCV incidence ************************************************
	ReadTable_Incidence_AgeDistr("project_disburdn_txfailure_input_inci_age_distr_byCDCacuteReport.txt");

	//vector<int> ageCategoryLowerBound_newInci;
	////ageCategoryLowerBound_newInci.push_back(10);
	//ageCategoryLowerBound_newInci.push_back(0);
	//ageCategoryLowerBound_newInci.push_back(5);
	//ageCategoryLowerBound_newInci.push_back(20);
	//ageCategoryLowerBound_newInci.push_back(30);
	//ageCategoryLowerBound_newInci.push_back(40);
	//ageCategoryLowerBound_newInci.push_back(50);
	//ageCategoryLowerBound_newInci.push_back(60);
	////ageCategoryLowerBound_newInci.push_back(70);
	//_modelParam._disBurdnData._vec_age_category_forNewIncidence = ageCategoryLowerBound_newInci;
	//vector<double> ageDistr;
	////ageDistr.push_back(0.064);
	////ageDistr.push_back(0.391);
	////ageDistr.push_back(0.315);
	////ageDistr.push_back(0.170);
	////ageDistr.push_back(0.052);
	////ageDistr.push_back(0.008);
	////ageDistr.push_back(0.000);
	//ageDistr.push_back(0.065);
	//	ageDistr.push_back(0.210);
	//	ageDistr.push_back(0.309);
	//	ageDistr.push_back(0.184);
	//	ageDistr.push_back(0.123);
	//	ageDistr.push_back(0.062);
	//	ageDistr.push_back(0.047);




	vector<int> ageCategoryLowerBound_initial;
	//ageCategoryLowerBound_initial.push_back(10); //0.4728
	ageCategoryLowerBound_initial.push_back(20); //0.4728
	ageCategoryLowerBound_initial.push_back(30); //0.3467
	ageCategoryLowerBound_initial.push_back(40); //0.2644
	ageCategoryLowerBound_initial.push_back(50); //0.1063
	ageCategoryLowerBound_initial.push_back(60); //0.0216	
	_modelParam._disBurdnData._vec_age_category_forInitialization = ageCategoryLowerBound_initial;
	vector<double> ageDistr2;
	// -----------------------------
	// ---- calibrated results -----
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation20_29"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation30_39"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation40_49"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation50_59"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation60"));



	_modelParam._disBurdnData._table_inci_age_distr[0] = ageDistr2;
	//_modelParam._disBurdnData._table_inci_age_distr[1] = ageDistr;
	//_modelParam._disBurdnData._table_inci_age_distr[2100] = ageDistr;



	// ****** probability of becoming aware by usual care ************************************************


	// [QC comment] override the probability of becoming aware from usual care:
	ReadTable_Prob_Aware("project_disburdn_txfailure_input_probAware.txt");
	for (map<typeStateDisBurdnModel, map<bool, map<double, double>>>::iterator it = _modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.begin();
		it != _modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.end(); it++) {
		for (map<bool, map<double, double>> ::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			for (map<double, double>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
				if (it->first == state_F0) {
					it3->second = GetValue("probBecomingAwareF0");
				}
				else if (it->first == state_F1) {
					it3->second = GetValue("probBecomingAwareF1");
				}
				else if (it->first == state_F2) {
					it3->second = GetValue("probBecomingAwareF2");
				}
				else if (it->first == state_F3) {
					it3->second = GetValue("probBecomingAwareF3");
				}
				else if (it->first == state_CoCirr) {
					it3->second = GetValue("probBecomingAwareF4");
				}
				else {
					ExitWithMsg("Error @ Initialize_StateX_parameter(): Incorrect fibrosis state to assign prob(aware) for usual care ");
				}

			}
		}
	}


	Initialize_SVR_Table("project_disburdn_txfailure_input_svr.txt", _modelParam._disBurdnData._table_SVR);


	_modelParam._disBurdnData._flag_extended_use_PR_PI = false;

	_modelParam._disBurdnData._prob_male = GetValue("distrMale");

	_modelParam._disBurdnData._prob_tx_coverage_by_medicarePartD = 1.0; // 0.9;
	//_modelParam._disBurdnData._prob_screen_accept = 1.0; // 0.819;		//.91 * .90, Rein 2012. (73.71-90.09)
	_modelParam._disBurdnData._prob_screen_accept = 0.819;		//.91 * .90, Rein 2012. (73.71-90.09)



	// ****** awareness in new HCV incidences ************************************************
	//double p_aware_in_new_incidence = 0;	
	
	_modelParam._disBurdnData._flag_override_awareness_for_initial_population = true;
	_modelParam._disBurdnData._awareness_for_initial_population_if_overrided = 0.5;
	
	double p_aware_in_new_incidence = GetValue("probAwareInNewInci");
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[40] = p_aware_in_new_incidence;// 0.0581;// 0.42 * 0.1383; // lower_bounds(): < 40
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[50] = p_aware_in_new_incidence;// 0.1950;// 0.71*0.2747;	// 40-49
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[60] = p_aware_in_new_incidence;// 0.1629;// 0.61*0.2670;	// 50-59
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[200] = p_aware_in_new_incidence;// 0.0577;// 0.27*0.2136; // 60-100
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[40] = p_aware_in_new_incidence;// 0.1763;// 0.53*0.3326;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[50] = p_aware_in_new_incidence;// 0.6045;// 0.915*0.6606;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[60] = p_aware_in_new_incidence;// 0.5459;// 0.85*0.6423;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[200] = p_aware_in_new_incidence;// 0.3494;// 0.68*0.5138;

	// assign the value for bookkeeping purposes, not used
	_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured = _modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured;
	_modelParam._disBurdnData._table_aware_prob_initial_population_insured = _modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured;


	// ****** fibrosis distribution in new HCV incidences ************************************************
	//// <fibrosis, probability>	
	_modelParam._disBurdnData._table_distr_fib[state_CoCirr] = GetValue("distrCirr"); //0; // 0.0554;
	_modelParam._disBurdnData._table_distr_fib[state_DeCirr] = GetValue("distrDC"); //0.00206;
	_modelParam._disBurdnData._table_distr_fib[state_DeCirr1yrPlus] = 0; //0.00928;
	_modelParam._disBurdnData._table_distr_fib[state_HCC] = GetValue("distrHCC"); //0; //0;
	_modelParam._disBurdnData._table_distr_fib[state_F0] = GetValue("distrF0"); //1.0; // 0.44494;
	_modelParam._disBurdnData._table_distr_fib[state_F1] = GetValue("distrF1"); // 0; //0.26474;
	_modelParam._disBurdnData._table_distr_fib[state_F2] = GetValue("distrF2"); //0; //0.1357;
	_modelParam._disBurdnData._table_distr_fib[state_F3] = GetValue("distrF3"); //0; //1 - 0.91212;




	// ------------------------------------------------
	// [QC comment] 
	// There will be no "insurance type" in EU5 analysis
	// "Pretend" all patients having "private" before age 65, and medicare after age 65
	// but no difference will be made by the insurance status.
	// ------------------------------------------------
	_modelParam._disBurdnData._table_distr_insurance[insr_uninsured] = 0; //0.267;
	_modelParam._disBurdnData._table_distr_insurance[insr_private] = 0; //0.498; // 0.7649 - 0.2666;
	_modelParam._disBurdnData._table_distr_insurance[insr_medicaid] = 0; //0.143; // 0.9078 - 0.7649;
	_modelParam._disBurdnData._table_distr_insurance[insr_military] = 0; //0.092; // 1.0 - 0.9078;
	_modelParam._disBurdnData._table_distr_insurance[insr_private] = 1.0;


	//// [QC comment] No ACA changes in pakistan analysis
	// Assume change is 0
	_modelParam._disBurdnData.ACATotalratio[2100] = 0.0;
	_modelParam._disBurdnData.ACAPrivateratio[2100] = 0.0;
	_modelParam._disBurdnData.ACAMedicaidratio[2100] = 0.0;
	_modelParam._disBurdnData.ACAPrivatechange[2100] = 0.0;
	_modelParam._disBurdnData.ACAMedicaidchange[2100] = 0.0;

	for (int y = ACA_START_YR; y <= END_YR; y++) {
		_modelParam._disBurdnData._table_distr_insurance_ACA_by_year.insert(pair<int, map<typeInsurance, double>>(y, _modelParam._disBurdnData._table_distr_insurance));
	}


	// ****** Treatment capacity *********************************************************************************** 
	ReadTable_2Col("project_disburdn_txfailure_input_tx_cost.txt", _modelParam._costData._table_tx_cost);
	for (int y = 1995; y <= 2009; y++) { _modelParam._costData._table_tx_cost[y] = 0; }
	for (int y = 2010; y <= END_YR; y++) { _modelParam._costData._table_tx_cost[y] = GetValue("daaTxCost"); }
	_modelParam._costData.scalingCoef = 1000000;

	ReadTable_2Col("project_disburdn_statex_input_tx_cap.txt", _modelParam._disBurdnData._table_tx_capacity);
	for (int y = 1995; y <= 2009; y++) { _modelParam._disBurdnData._table_tx_capacity[y] = 0; }
	for (int y = 2010; y <= 2016; y++) { _modelParam._disBurdnData._table_tx_capacity[y] = _paramPK_tx_cap_statusquo; }
	for (int y = 2017; y <= END_YR; y++) { _modelParam._disBurdnData._table_tx_capacity[y] = _paramPK_tx_cap_future; }


	// ****** Screening rate *********************************************************************************** 
	// [QC comment] Below will not be used. On the safe side, set all of them = 0
	_modelParam._disBurdnData._table_universal_scr_rate.clear();
	_modelParam._disBurdnData._table_universal_scr_rate[-1] = 0;//0;
	_modelParam._disBurdnData._table_universal_scr_rate[0] = _paramPK_universal_screeningRate;
	_modelParam._disBurdnData._table_universal_scr_rate[100] = _paramPK_universal_screeningRate;//0;

	_modelParam._disBurdnData._table_universal_scr_cap.clear();
	_modelParam._disBurdnData._table_universal_scr_cap[-1] = 0;//0;
	//_modelParam._disBurdnData._table_universal_scr_cap[0] = _paramPK_universal_screeningCap;
	_modelParam._disBurdnData._table_universal_scr_cap[(int)GetValue("screeningStartYr") - START_YR_SCREEN] = 0;//0;
	_modelParam._disBurdnData._table_universal_scr_cap[(int)GetValue("screeningStartYr") - START_YR_SCREEN + 1] = _paramPK_universal_screeningCap;//0;

	_modelParam._disBurdnData._table_universal_scr_cap[100] = _paramPK_universal_screeningCap;//0;



	_modelParam._disBurdnData._record_ageDistr_year = 2010;
	_modelParam._disBurdnData._record_ageDistr_outputFolder = "output_project_statex";

	return 0;
}


int Project_DisBurd_StateX::ReadTable_Param()
{
	ifstream inFile;
	string argStr = "project_disburdn_statex_input_params.txt";
	inFile.open(argStr.c_str());		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	// read table header, register countries


	string str;
	while (inFile >> str) {
		if ("//" == str) {
			getline(inFile, str);
			continue;
		}
		inFile >> _paramPK_map[str];
	}
	inFile.close();
	return 0;


}

double Project_DisBurd_StateX::GetValue(string argVarName)
{
	map<string, double > ::const_iterator it = _paramPK_map.find(argVarName);
	if (it == _paramPK_map.end()) {
		ExitWithMsg("Error @ Project_DisBurd_StateX::GetValue(string argVarName): Can't find parameter: " + argVarName);
	}
	return it->second;

}


int Project_DisBurd_StateX::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
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
		argModelParam._transData.pr_SVR_Delta_DAA = 1.0 - argVal / argModelParam._transData._pr_svr_for_SA;
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





int Project_DisBurd_StateX::BaseCaseRun()
{
	BurdenModelSim myModel;
	myModel.RedfineStartYear(_paramPK_start_year);
	myModel.Run(_modelParam);

	//getchar();
	return 0;
}


int Project_DisBurd_StateX::ROI_Analysis(int argIdx)
{
	string fileStr = "output_project_statex/roi_" + basicToStr(argIdx) + ".txt";
	ofstream out_target;
	out_target.open(fileStr.c_str());
	if (out_target.fail()) {
		ExitWithMsg("Error @ RunBaseCase(): output_project_statex/roi_" + basicToStr(argIdx) + ".txt");
	}


	int cycle_base = 2017 - START_YR;


	// creat a list of treatment capacity
	vector<int> vecTxCap;
	vecTxCap.push_back(0); // no treatment
	vecTxCap.push_back(GetValue("txCapStatusQuo")); // statusquo
	//for (int l = 0; l < 10; l++) {
	//	vecTxCap.push_back(1000+ 1000 * (l));
	//}
	vecTxCap.push_back(1000);
	vecTxCap.push_back(6000);

	for (int m = 0; m < 2; m++) {
		for (int c = 0; c <20; c++) {
			double daaCost = 3000 + c * 1000;
			

			vector<double> vecFixedResults;
			for (int l = 0; l < (int)vecTxCap.size(); l++) {
				_paramPK_universal_screeningRate = 0;
				// _paramPK_universal_screeningCap = l<2 ? 0 : 100 * (m + 1);
				_paramPK_universal_screeningCap = l<2 ? 0 :  (500 + m * 1500);

				_paramPK_tx_cap_future = vecTxCap[l];
				_paramPK_map["daaTxCost"] = daaCost;
				
				SetOutFile("output_project_statex/roi/roi_"
					+ basicToStr(_paramPK_universal_screeningRate) + "_"
					+ basicToStr(_paramPK_universal_screeningCap) + "_"
					+ basicToStr(_paramPK_tx_cap_future) + ".txt");


				Initialize_StateX_Parameter();
				SetScreeningOption(typeScreeningScenario::screen_universal_capacity);

				SIM_SEED = 2;
				BurdenModelSim myModel;
				myModel.RedfineStartYear(_paramPK_start_year);
				myModel.Run(_modelParam);

				const disBurdnCounterType simData = myModel.GetSimCounters();



				vector<double> vecCostScr = simData._counterCost_byCategory.find("c_screening")->second;
				vector<double> vecCostTx = simData._counterCost_byCategory.find("c_treatment")->second;
				vector<double> vecCostDis = simData._counterCost_byCategory.find("c_natr_hist")->second;

				if (l >= 2) {
					out_target << fixed << setprecision(4) << _paramPK_universal_screeningRate << "\t"
						<< _paramPK_universal_screeningCap << "\t"
						<< daaCost << "\t"
						<< _paramPK_tx_cap_future << "\t";
					//for (int idx = 0; idx < (int)vecFixedResults.size(); idx++) {
					for (int idx = 0; idx < (int)vecFixedResults.size()/2; idx++) {
						out_target << vecFixedResults[idx] << "\t";
					}
				}

					int horz[] = { 5,10 };
					for (int h = 0; h < sizeof(horz) / sizeof(int); h++) {
						double total_num_tx = 0;
						double total_num_scr = 0;
						double total_cost = 0;
						double total_cost_excl_scr = 0;
						double tx_cost = 0;
						for (int k = cycle_base; k < cycle_base + horz[h]; k++) {
							total_num_tx += simData._counterTxFailure_numTx[k];
							total_num_scr += simData._counter_screening[k];
							total_cost_excl_scr += (vecCostTx[k] + vecCostDis[k]);
							total_cost += (vecCostTx[k] + vecCostDis[k] + vecCostScr[k]);
							tx_cost += vecCostTx[k];
						}
						if (l <= 1) {
							vecFixedResults.push_back(total_num_tx);
							vecFixedResults.push_back(total_cost_excl_scr);
							vecFixedResults.push_back(total_cost);							
							vecFixedResults.push_back(tx_cost);

						}
						else {
							out_target << total_num_tx << "\t"
								<< total_cost_excl_scr << "\t"
								<< total_cost <<"\t"
								<< tx_cost << "\t";
						}
					}// end of horizons
					if(l>=2) out_target << endl;

				}

			// end of various treatment capacity
		}// end of treatment cost
	}// end of screening capacities

	out_target.close();
	return 0;
}




int Project_DisBurd_StateX::CompareCost()
{
	string fileStr = "output_project_statex/compare_cost.txt";
	ofstream out_target;
	out_target.open(fileStr.c_str());
	if (out_target.fail()) {
		ExitWithMsg("Error @ RunBaseCase(): output_project_statex/compare_cost.txt");
	}


	int cycle_base = 2017 - START_YR;

	int listScrCap[] = { 500, 2000 };
	int listTxCap[] = { 1000, 6000 };
	double listTxCost[] = {0, 5000, 10000, 15000 };


	vector<vector<double>> resultsVec2d;

	for (int idxScr = 0; idxScr < sizeof(listScrCap) / sizeof(int); idxScr++) {
		for (int idxTx = 0; idxTx < sizeof(listTxCap) / sizeof(int); idxTx++) {
			vector<double> resultsPrev1, resultsPrev2;

			for (int idxCost = 0; idxCost < sizeof(listTxCost) / sizeof(double); idxCost++) {

				_paramPK_universal_screeningRate = 0;
				_paramPK_universal_screeningCap = idxCost == 0 ? 0 : listScrCap[idxScr];
				_paramPK_tx_cap_future = idxCost == 0 ? 0 : listTxCap[idxTx];
				_paramPK_map["daaTxCost"] = listTxCost[idxCost];

				SetOutFile("output_project_statex/compare_cost/cmpcost_"
					+ basicToStr(_paramPK_universal_screeningRate) + "_"
					+ basicToStr(_paramPK_universal_screeningCap) + "_"
					+ basicToStr(_paramPK_tx_cap_future) + "_"
					+ basicToStr(GetValue("daaTxCost")) + ".txt");


				Initialize_StateX_Parameter();
				SetScreeningOption(typeScreeningScenario::screen_universal_capacity);

				SIM_SEED = 2;
				BurdenModelSim myModel;
				myModel.RedfineStartYear(_paramPK_start_year);
				myModel.Run(_modelParam);

				const disBurdnCounterType simData = myModel.GetSimCounters();



				vector<double> vecCostScr = simData._counterCost_byCategory.find("c_screening")->second;
				vector<double> vecCostTx = simData._counterCost_byCategory.find("c_treatment")->second;
				vector<double> vecCostDis = simData._counterCost_byCategory.find("c_natr_hist")->second;

				vector<double> oneRow;
				for (int k = 0; k < END_YR - START_YR; k++) {
					oneRow.push_back(vecCostDis[k] + vecCostTx[k] + vecCostScr[k]);

					if (idxCost == 0) {
						resultsPrev1.push_back(simData._counterHCV[k]);
					}
					else if (idxCost == 1) {
						resultsPrev2.push_back(simData._counterHCV[k]);
					}
				}
				resultsVec2d.push_back(oneRow);

				
			}// cost
			resultsVec2d.push_back(resultsPrev1);
			resultsVec2d.push_back(resultsPrev2);
		}// tx
	}// scr
			

	int nrows = resultsVec2d.size();
	int ncols = resultsVec2d[0].size();
	for (int c = 0; c < ncols; c++) {
		out_target << fixed<<setprecision(4) << START_YR + c <<"\t"; // year
		for (int r = 0; r < nrows; r++) {
			out_target << resultsVec2d[r][c] << "\t";
		}
		out_target << endl;
	}
	out_target.close();
	return 0;
}


int Project_DisBurd_StateX::CalibInitialSizeAndInciRate(int argIdx)
{
	string filestr = "output_project_pakistan/calib/output_calibrate_n0_inciRate_" + basicToStr(argIdx) + ".txt";
	ofstream out_target;
	out_target.open(filestr.c_str());
	if (out_target.fail()) {
		ExitWithMsg("Error @ CalibInitialSizeAndInciRate(): " + filestr);
	}
	out_target.close();


	//double listInciRate[] = {0.02, 0.025,0.03, 0.035, 0.04, 0.045, 0.05, 0.06};
	//int listN0[] = {40, 45, 50, 55, 60, 65, 70,75, 80,85,90,95,100,105,110};

	double listInciRate[] = { 0.025, 0.026,0.027,0.028,0.029,0.030,0.031,0.032,0.033,0.034,0.035 };
	int listN0[] = { 70,75, 80,85,90,95,100,105,110,115,120 };


	//for (int m = 0; m < sizeof(listInciRate) / sizeof(double); m++) {
	int m = argIdx - 1; // idx = 1,2,... 11
	for (int l = 0; l < sizeof(listN0) / sizeof(int); l++) {
		_paramPK_map["incidenceProportion"] = listInciRate[m];
		_paramPK_map["numInitialPpl"] = listN0[l] * 100000;


		SetOutFile("output_project_pakistan/calib/out_calib_" + basicToStr(listInciRate[m]) + "_" + basicToStr(listN0[l]) + ".txt");

		Initialize_StateX_Parameter();
		BurdenModelSim myModel;
		myModel.Run(_modelParam);
		const disBurdnCounterType simData = myModel.GetSimCounters();

		/////// write output file /////////
		ofstream out_target;
		out_target.open(filestr.c_str(), ios::app);
		if (out_target.fail()) {
			ExitWithMsg("Error @ CalibInitialSizeAndInciRate(): " + filestr);
		}

		out_target << fixed << setprecision(4) << listInciRate[m] << "\t" << listN0[l] << "\t";
		out_target << simData._counterHCV[2008 - START_YR] << "\t"
			<< simData._counter_new_hcv_incidence[2014 - START_YR] << "\t";
		out_target << endl;
		out_target.close();
	}
	//}

	return 0;
}




int Project_DisBurd_StateX::ReadTable_NS5AMarketShare(string argStr)
{
	// Read the NS5A/nonNS5A market share input data
	ifstream inFileNS5AMarketShare;
	inFileNS5AMarketShare.open(argStr);
	while (!inFileNS5AMarketShare.eof()) {
		int year;
		map<int, double> share;
		inFileNS5AMarketShare >> year;
		double s1, s2, s3, s4;
		inFileNS5AMarketShare >> s1 >> s2 >> s3 >> s4;
		share[1] = s1;
		share[2] = s2;
		share[3] = s3;
		share[456] = s4;
		_modelParam._disBurdnData._table_ns5a_marketshare[year] = share;
	}
	//cout<<"testing: _tfa_market_share = " << _modelParam._disBurdnData._table_ns5a_marketshare[2015][1]<<"\t"<<_modelParam._disBurdnData._table_ns5a_marketshare[2016][3]<<endl;
	inFileNS5AMarketShare.close();
	return 0;
}

int Project_DisBurd_StateX::ReadTable_BackgroundMortality_5yAgeGroup(string argStr, string argCountry)
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


int Project_DisBurd_StateX::ReadTable_Prob_Aware(string argStr)
{
	int vfib;
	bool vInsured;
	double vAge;
	double vProb;
	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	string str;
	while (inFile >> str) {
		if ("//" == str) { getline(inFile, str); continue; }
		int tmp = atoi(str.c_str());
		for (int f = 0; f < 5; f++) {
			map<bool, map<double, double>> mp2;
			for (int b = 1; b >= 0; b--) {
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





int Project_DisBurd_StateX::ReadTable_BackgroundMortality(string argStr)
{
	ifstream inFile2;
	inFile2.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile2.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

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



int Project_DisBurd_StateX::Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap)
{
	// read the table from file
	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

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

int Project_DisBurd_StateX::ReadTable_Incidence_AgeDistr(string argStr)
{
	_modelParam._disBurdnData._vec_age_category_forInitialization.clear();
	_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(18);
	_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(24);
	_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(34);
	_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(44);
	_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(54);
	_modelParam._disBurdnData._vec_age_category_forInitialization.push_back(64);

	_modelParam._disBurdnData._vec_age_category_forNewIncidence.clear();
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(18);
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(20);
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(30);
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(40);
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(50);
	_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(60);
	if (argStr == "project_disburdn_txfailure_input_inci_age_distr.txt") {
		_modelParam._disBurdnData._vec_age_category_forNewIncidence.push_back(70);
	}
	else if (argStr == "project_disburdn_txfailure_input_inci_age_distr_byCDCacuteReport.txt") {
		// do not include age category 70+
	}
	else {
		ExitWithMsg("Project_TxFailure::ReadTable_Incidence_AgeDistr(string argStr): Recheck the age categories for the age incidence data file: " + argStr);
	}

	ReadTable(argStr, _modelParam._disBurdnData._table_inci_age_distr);
	return 0;
}
