#include "project_disburdn_pakistan.h"
using namespace std;

Project_DisBurd_Pakistan::Project_DisBurd_Pakistan()
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


	// ------------ additional pakistan parameters ------------------

	ReadTable_Param();

	_paramPK_tx_cap_future = GetValue("txCapFrom2017");
	_paramPK_universal_screeningRate = GetValue("screeningRate");
	_paramPK_universal_screeningCap = GetValue("screeningCap");
}


int Project_DisBurd_Pakistan::Initialize_Pakistan_Parameter(type_pakistan_incidence argInciMode)
{
	//_paramPK_inci_mode = INCI_CONSTANT;
	//_paramPK_inci_mode = INCI_INCREASING;
	//_paramPK_inci_mode = INCI_PROPORTIONAL;
	_paramPK_inci_mode = argInciMode;

	_paramPK_inci_const = GetValue("incidenceConst");
	_paramPK_inci_proportion = GetValue("incidenceProportion");

	_modelParam._disBurdnData._num_initial_population = GetValue("numInitialPpl"); // initial population at year 1995

	ReadTable_BackgroundMortality_5yAgeGroup("Input_mortality_male_female_pakistan.in", "PK");
	//ReadTable_BackgroundMortality("Input_mortality_male_female_pakistan_wrong.in");	



	// ****** NS5A market share: simplified with 0 ****************************************************
	//ReadTable_NS5AMarketShare("project_disburdn_txfailure_input_ns5a_marketshare.txt");
	//ReadTable_NS5AMarketShare("project_disburdn_eu5_input_ns5a_marketshare_"+argStrCountry+".txt");
	//ReadTable_ExtendedUsePRandPI("project_disburdn_eu5_input_extended_use_pr_pi.txt", argStrCountry);
	for(int y=START_YR; y<=END_YR; y++){
		_modelParam._disBurdnData._table_ns5a_marketshare[y][1] = 0 ;
		_modelParam._disBurdnData._table_ns5a_marketshare[y][2] = 0 ;
		_modelParam._disBurdnData._table_ns5a_marketshare[y][3] = 0 ;
		_modelParam._disBurdnData._table_ns5a_marketshare[y][456] = 0 ;
	} // assuming all non-NS5A


	// ****** Treatment timeline ***************************************************************
	// PEG-RBV: 1998 - 2014
	// DAA:     2015 - end
	// See chanages in the section of "#elif defined(__PAKISTAN_SETTING__)" in disburdn_config.h

	// ****** Incidence of new HCV cases ************************************************
	//ReadTable_Incidence_Num("project_disburdn_txfailure_input_incidence.txt");
	//ReadTable_Incidence_Num("project_disburdn_eu5_input_incidence.txt", argStrCountry);
	if (_paramPK_inci_mode == INCI_PROPORTIONAL){
		_modelParam._disBurdnData._flag_inci_proportional_to_prevalence = true;
		for(int y=START_YR; y<=END_YR; y++){_modelParam._disBurdnData._table_incidence_proportion[y]=_paramPK_inci_proportion;}

	}else if(_paramPK_inci_mode == INCI_CONSTANT){
		for(int y=START_YR; y<=END_YR; y++){_modelParam._disBurdnData._table_incidence_proportion[y]=_paramPK_inci_proportion;}
		_modelParam._disBurdnData._year_start_constant_new_incidence = 2015;
		for(int y=_modelParam._disBurdnData._year_start_constant_new_incidence; y<=END_YR; y++){
			_modelParam._disBurdnData._table_incidence[y]=_paramPK_inci_const;
		}
	}else if (_paramPK_inci_mode == INCI_INCREASING){
		for(int y=START_YR; y<=END_YR; y++){_modelParam._disBurdnData._table_incidence_proportion[y]=_paramPK_inci_proportion;}
		_modelParam._disBurdnData._year_start_constant_new_incidence = 2015;
		for(int y=_modelParam._disBurdnData._year_start_constant_new_incidence; y<=END_YR; y++){
			_modelParam._disBurdnData._table_incidence[y]=(int)(_paramPK_inci_const * (1.0+0.02 * (y-_modelParam._disBurdnData._year_start_constant_new_incidence)));
		}

	}else{
		ExitWithMsg("Error @ Initialize_Pakistan_parameter(): Unknown incidence mode parameter = "+basicToStr(_paramPK_inci_mode));
	}

	// ****** Age distribution of new HCV incidence ************************************************
	//ReadTable_Incidence_AgeDistr("project_disburdn_txfailure_input_inci_age_distr.txt");
	//ReadTable_Incidence_AgeDistr("project_disburdn_eu5_input_inci_age_distr_"+argStrCountry+".txt");
	
	vector<int> ageCategoryLowerBound_newInci;
	//ageCategoryLowerBound_newInci.push_back(10);
	ageCategoryLowerBound_newInci.push_back(0);
	ageCategoryLowerBound_newInci.push_back(5);
	ageCategoryLowerBound_newInci.push_back(20);
	ageCategoryLowerBound_newInci.push_back(30);
	ageCategoryLowerBound_newInci.push_back(40);
	ageCategoryLowerBound_newInci.push_back(50);
	ageCategoryLowerBound_newInci.push_back(60);
	//ageCategoryLowerBound_newInci.push_back(70);
	_modelParam._disBurdnData._vec_age_category_forNewIncidence = ageCategoryLowerBound_newInci;
	vector<double> ageDistr;
	//ageDistr.push_back(0.064);
	//ageDistr.push_back(0.391);
	//ageDistr.push_back(0.315);
	//ageDistr.push_back(0.170);
	//ageDistr.push_back(0.052);
	//ageDistr.push_back(0.008);
	//ageDistr.push_back(0.000);
	ageDistr.push_back(0.065);
		ageDistr.push_back(0.210);
		ageDistr.push_back(0.309);
		ageDistr.push_back(0.184);
		ageDistr.push_back(0.123);
		ageDistr.push_back(0.062);
		ageDistr.push_back(0.047);




	vector<int> ageCategoryLowerBound_initial;
	//ageCategoryLowerBound_initial.push_back(10); //0.4728
	ageCategoryLowerBound_initial.push_back(0); //0.4728
	ageCategoryLowerBound_initial.push_back(5); //0.4728
	ageCategoryLowerBound_initial.push_back(20); //0.2132
	ageCategoryLowerBound_initial.push_back(30); //0.3467
	ageCategoryLowerBound_initial.push_back(40); //0.2644
	ageCategoryLowerBound_initial.push_back(50); //0.1063
	ageCategoryLowerBound_initial.push_back(60); //0.0216	
	_modelParam._disBurdnData._vec_age_category_forInitialization = ageCategoryLowerBound_initial;
	vector<double> ageDistr2;
	// -----------------------------
	// ---- calibrated results -----
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation0_5"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation5_20"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation20_30"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation30_40"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation40_50"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation50_60"));
	ageDistr2.push_back(GetValue("ageDistrInitialPopulation60"));
	


	_modelParam._disBurdnData._table_inci_age_distr[0] = ageDistr2;
	_modelParam._disBurdnData._table_inci_age_distr[1] = ageDistr;
	_modelParam._disBurdnData._table_inci_age_distr[2100] = ageDistr;




	// [QC comment] override the probability of becoming aware from usual care:
	ReadTable_Prob_Aware("project_disburdn_txfailure_input_probAware.txt");
	for (map<typeStateDisBurdnModel, map<bool, map<double, double>>>::iterator it = _modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.begin();
		it != _modelParam._disBurdnData._table_prob_becoming_aware_by_fib_insur_age.end(); it++) {
			for (map<bool, map<double, double>> ::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
				for (map<double, double>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
					if(it->first == state_F0){
						it3->second = GetValue("probBecomingAwareF0");
					}else if (it->first == state_F1){
						it3->second = GetValue("probBecomingAwareF1");
					}else if (it->first == state_F2){
						it3->second = GetValue("probBecomingAwareF2");
					}else if (it->first == state_F3){
						it3->second = GetValue("probBecomingAwareF3");
					}else if (it->first == state_CoCirr){
						it3->second = GetValue("probBecomingAwareF4");
					}else{
						ExitWithMsg("Error @ Initialize_Pakistan_parameter(): Incorrect fibrosis state to assign prob(aware) for usual care ");
					}
					
				}
			}
	}


	if (RECOVER_SALVAGE_POSTER_VERSION) {
		Initialize_SVR_Table("project_disburdn_txfailure_input_svr_verPoster.txt", _modelParam._disBurdnData._table_SVR);
	}
	else {
		Initialize_SVR_Table("project_disburdn_txfailure_input_svr.txt", _modelParam._disBurdnData._table_SVR);
	}


	// ---- update disBurdnDataType variables for EU5 analysis ------------
	_modelParam._disBurdnData._flag_extended_use_PR_PI = false;

	_modelParam._disBurdnData._prob_male = 0.5248; // Qureshi 2010; see parameters_forpakistan.xlsx
	_modelParam._disBurdnData._prob_receivable = 0.915;
	//_prob_TrR = 1.0; // unchanged
	_modelParam._disBurdnData._prob_chronic = 0.741; // see parameters forpakistan.xlsx
	_modelParam._disBurdnData._prob_tx_coverage_by_medicarePartD = 1.0; // 0.9;
	_modelParam._disBurdnData._prob_naive_initialization = 1.0; // unchanged // 0.65, or 0.17???
	_modelParam._disBurdnData._prob_chronic_contra = 0.346; // unchanged
	_modelParam._disBurdnData._prob_mod_given_contra = 0.231 / 0.346; // unchanged
	//_modelParam._disBurdnData._prob_screen_accept = 1.0; // 0.819;		//.91 * .90, Rein 2012. (73.71-90.09)
	_modelParam._disBurdnData._prob_screen_accept = 0.819;		//.91 * .90, Rein 2012. (73.71-90.09)
	//									//percentage of people who accept to be screened for HepC
	//									//Rein 2012 assumed that 28%of the chronic cases are aware of their infection, so screening would not benefit them.
	//_prob_birthCohortScr_insured = 0.9; // unchanged
	//_prob_birthCohortScr_uninsured = 0.1; // unchanged
	//_prob_delay_tx_until_wave1_f0f2 = 0.75; // unchanged
	//_prob_delay_tx_until_wave1_f3 = 0.25; // unchanged


	//double p_aware_in_new_incidence = 0;
	double p_aware_in_new_incidence = GetValue("probAwareInNewInci"); // assumption for Pakistan // in the US: diagnosed rate = 57% [Razavi2014]	
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[40] = p_aware_in_new_incidence;// 0.0581;// 0.42 * 0.1383; // lower_bounds(): < 40
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[50] = p_aware_in_new_incidence;// 0.1950;// 0.71*0.2747;	// 40-49
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[60] = p_aware_in_new_incidence;// 0.1629;// 0.61*0.2670;	// 50-59
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[200] = p_aware_in_new_incidence;// 0.0577;// 0.27*0.2136; // 60-100
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[40] = p_aware_in_new_incidence;// 0.1763;// 0.53*0.3326;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[50] = p_aware_in_new_incidence;// 0.6045;// 0.915*0.6606;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[60] = p_aware_in_new_incidence;// 0.5459;// 0.85*0.6423;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[200] = p_aware_in_new_incidence;// 0.3494;// 0.68*0.5138;

	_modelParam._disBurdnData._table_aware_prob_initial_population_uninsured = _modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured;
	_modelParam._disBurdnData._table_aware_prob_initial_population_insured = _modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured;


	// [QC comment, 03/07/2017]
	// Read country specific genotype 
	// <genotype-int, probability>
	_modelParam._disBurdnData._table_distr_genotype[1] = GetValue("distrG1"); // see parameters_forpakistan.xlsx
	_modelParam._disBurdnData._table_distr_genotype[2] = GetValue("distrG2"); // see parameters_forpakistan.xlsx
	_modelParam._disBurdnData._table_distr_genotype[3] = GetValue("distrG3"); // see parameters_forpakistan.xlsx
	_modelParam._disBurdnData._table_distr_genotype[456] = GetValue("distrG456"); // see parameters_forpakistan.xlsx



	// <fibrosis, probability>	
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


	// ****** Treatment capacity *********************************************************************************** 
	// 1998-2010: 10000
	// 2011-2016: 85000
	// 2017+    : variable
	//ReadTable_TreatmentCapacity("project_disburdn_txfailure_input_tx_cap.txt");
	//ReadTable_TreatmentCapacity("project_disburdn_eu5_input_tx_cap.txt", argStrCountry);

	// --- v0 -----
	//for(int y=1995; y<=1997; y++){_modelParam._disBurdnData._table_tx_capacity[y] = 0;}
	//for(int y=1998; y<=2015; y++){_modelParam._disBurdnData._table_tx_capacity[y] = 10000;}	
	//for(int y=2016; y<=END_YR; y++){_modelParam._disBurdnData._table_tx_capacity[y] = _paramPK_tx_cap_2016onwards;}
	// --- v1: email correspondence with Naveed ----
	for(int y=1995; y<=2003; y++){_modelParam._disBurdnData._table_tx_capacity[y] = 0;}
	_modelParam._disBurdnData._table_tx_capacity[2004] = 11809;
	_modelParam._disBurdnData._table_tx_capacity[2005] = 19828;
	_modelParam._disBurdnData._table_tx_capacity[2006] = 30675;
	_modelParam._disBurdnData._table_tx_capacity[2007] = 50857;
	_modelParam._disBurdnData._table_tx_capacity[2008] = 70499;
	for(int y=2009; y<=2014; y++){_modelParam._disBurdnData._table_tx_capacity[y] = 85000;}
	_modelParam._disBurdnData._table_tx_capacity[2015] = 65385;
	_modelParam._disBurdnData._table_tx_capacity[2016] = 160650;
	_modelParam._disBurdnData._table_tx_capacity[2017] = 160650;
	for(int y=2018; y<=END_YR; y++){_modelParam._disBurdnData._table_tx_capacity[y] = _paramPK_tx_cap_future;}
	// ****** Screening rate *********************************************************************************** 

	// [QC comment] Below will not be used. On the safe side, set all of them = 0
	_modelParam._disBurdnData._table_universal_scr_rate.clear();
	_modelParam._disBurdnData._table_universal_scr_rate[-1] = 0;//0;
	_modelParam._disBurdnData._table_universal_scr_rate[0] = _paramPK_universal_screeningRate;
	_modelParam._disBurdnData._table_universal_scr_rate[100] = _paramPK_universal_screeningRate;//0;

	_modelParam._disBurdnData._table_universal_scr_cap.clear();
	_modelParam._disBurdnData._table_universal_scr_cap[-1] = 0;//0;
	_modelParam._disBurdnData._table_universal_scr_cap[0] = _paramPK_universal_screeningCap;
	_modelParam._disBurdnData._table_universal_scr_cap[100] = _paramPK_universal_screeningCap;//0;



	_modelParam._transData.pr_DeCirr_LivTr = 0.0;
	_modelParam._transData.pr_HCC_LivTr = 0.0;

	return 0;
}


int Project_DisBurd_Pakistan::ReadTable_Param()
{
	ifstream inFile;
	string argStr	= "project_disburdn_pakistan_input_params.txt";
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
		inFile>> _paramPK_map[str];
	}
	inFile.close();
	return 0;


}

double Project_DisBurd_Pakistan::GetValue(string argVarName)
{
	map<string, double > ::const_iterator it = _paramPK_map.find(argVarName);
	if (it == _paramPK_map.end()) {
		ExitWithMsg("Error @ Project_DisBurd_Pakistan::GetValue(string argVarName): Can't find parameter: "+argVarName);
	}
	return it->second;

}


int Project_DisBurd_Pakistan::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
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





int Project_DisBurd_Pakistan::BaseCaseRun()
{
	BurdenModelSim myModel;
	myModel.Run(_modelParam);
	
	//getchar();
	return 0;
}


int Project_DisBurd_Pakistan::TestWHOTarget(int argIdx)
{
	string fileStr = "output_project_pakistan/target/WHO_target_" + basicToStr(argIdx) + ".txt";
	ofstream out_target;
	out_target.open(fileStr.c_str());
	if (out_target.fail()) {
		ExitWithMsg("Error @ RunBaseCase(): output_project_pakistan/target/WHO_target_" + basicToStr(argIdx) + ".txt");
	}
	out_target.close();

	int cycle_base = 2015 - START_YR;
	int cycle_target = 2030 - START_YR;

	//double listScrRate[] = { 0.05,0.1,0.15,0.2,0.3,0.4 };
	//int listTxCap[] = {100000, 300000,500000, 700000, 900000};
		
	//double listScrRate[] = {0.1, 0.15, 0.2, 0.25,0.3};
	//int listTxCap[] = {300000, 400000, 500000, 600000, 700000,900000};



	//for (int m = 0; m < sizeof(listScrRate) / sizeof(double); m++) {
	//	for (int l = 0; l < sizeof(listTxCap) / sizeof(int); l++) {
	//		_paramPK_universal_screeningRate = listScrRate[m];
	//		_paramPK_tx_cap_2016onwards = listTxCap[l];

	//for(int m=0; m<10; m++){
	int m = argIdx - 1; // argIdx = 1, 2, ..., 10
		for(int l=0; l<10; l++){
			_paramPK_universal_screeningRate = 0;
			_paramPK_universal_screeningCap = 400000 + 50000 * (m + 1);
			_paramPK_tx_cap_future = 200000+ 100000 * (l + 1);

			SetOutFile("output_project_pakistan/target/out_who_target_"
				+basicToStr(_paramPK_universal_screeningRate)+"_"
				+basicToStr(_paramPK_universal_screeningCap/1000)+"_"
				+basicToStr(_paramPK_tx_cap_future/1000)+".txt");


			Initialize_Pakistan_Parameter();
			SetScreeningOption(typeScreeningScenario::screen_universal_capacity);

			SIM_SEED = 2;

			BurdenModelSim myModel;
			myModel.Run(_modelParam);

			const disBurdnCounterType simData = myModel.GetSimCounters();			
			// WHO target 1: # new infection
			double target1_reduction_new_infection = 1.0 - ((double)simData._counter_new_hcv_incidence[cycle_target]) / simData._counter_new_hcv_incidence[cycle_base];

			// WHO target 2: # liver-related deaths
			double target2_reduction_lrd = 1.0 - ((double)simData._counter_liverRelatedDeath[cycle_target]) / simData._counter_liverRelatedDeath[cycle_base];

			// WHO target 3: % treated among eligible people
			double target3_treated_pct_among_eligible = ((double)simData._counter_treated_alive[cycle_target]) / (simData._counter_txElig_all[cycle_target] + simData._counter_treated_alive[cycle_target]);

			// WHO target 4: % awareness/diagnose rate
			//double target4_awareness_rate = ((double)simData._counter_aware[cycle_target]) / simData._counterHCV[cycle_target];
			double target4_awareness_rate = ((double)simData._counter_aware[cycle_target] + simData._counterHCV_complement[cycle_target]) / ((double)simData._counterHCV[cycle_target] + simData._counterHCV_complement[cycle_target]);


			// actual number of screening and treatment
			int actualNumScr = 0;
			int actualNumTx = 0;
			for (int k = cycle_base; k <= cycle_target; k++) {
				actualNumScr += simData._counter_screening[k];
				actualNumTx += simData._counterTxFailure_numTx[k];
			}

			/////// write output file /////////
			ofstream out_target;
			out_target.open(fileStr.c_str(), ios::app);
			if (out_target.fail()) {
				ExitWithMsg("Error @ RunBaseCase(): "+fileStr);
			}

			out_target <<fixed <<setprecision(4) << _paramPK_universal_screeningRate << "\t" << _paramPK_universal_screeningCap << "\t" << _paramPK_tx_cap_future << "\t";
			out_target << actualNumScr << "\t" << actualNumTx << "\t";
			out_target << target1_reduction_new_infection << "\t"
				<< target2_reduction_lrd << "\t"
				<< target3_treated_pct_among_eligible << "\t"
				<< target4_awareness_rate << "\t";
			out_target << endl;
			out_target.close();
		}
	//}

	return 0;
}


int Project_DisBurd_Pakistan::CalibInitialSizeAndInciRate(int argIdx)
{
	string filestr = "output_project_pakistan/calib/output_calibrate_n0_inciRate_"+basicToStr(argIdx)+".txt";
	ofstream out_target;
	out_target.open(filestr.c_str());
	if (out_target.fail()) {
		ExitWithMsg("Error @ CalibInitialSizeAndInciRate(): " + filestr);
	}
	out_target.close();


	//double listInciRate[] = {0.02, 0.025,0.03, 0.035, 0.04, 0.045, 0.05, 0.06};
	//int listN0[] = {40, 45, 50, 55, 60, 65, 70,75, 80,85,90,95,100,105,110};
	
	double listInciRate[] = {0.030,0.031,0.032,0.033,0.034,0.035};
	int listN0[] = {900,905, 910, 915, 920,925,930,935,940,945,950,955,960,965,970,975,980,985,990,1000}; // in 10,000

	
	//for (int m = 0; m < sizeof(listInciRate) / sizeof(double); m++) {
	int m = argIdx -1 ; // idx = 1,2,... 11
		for (int l = 0; l < sizeof(listN0) / sizeof(int); l++) {
			_paramPK_map["incidenceProportion"] = listInciRate[m];
			_paramPK_map["numInitialPpl"] = listN0[l] * 10000;
			

			SetOutFile("output_project_pakistan/calib/out_calib_" + basicToStr(listInciRate[m]) + "_" + basicToStr(listN0[l]) + ".txt");

			Initialize_Pakistan_Parameter();
			BurdenModelSim myModel;
			myModel.Run(_modelParam);
			const disBurdnCounterType simData = myModel.GetSimCounters();

			/////// write output file /////////
			ofstream out_target;
			out_target.open(filestr.c_str(), ios::app);
			if (out_target.fail()) {
				ExitWithMsg("Error @ CalibInitialSizeAndInciRate(): "+filestr);
			}

			out_target << fixed << setprecision(4) << listInciRate[m] << "\t" << listN0[l] << "\t";
			out_target << simData._counterHCV[2008-START_YR] << "\t" 
				<< simData._counter_new_hcv_incidence[2014-START_YR] << "\t";
			out_target << endl;
			out_target.close();
		}
	//}

	return 0;
}




int Project_DisBurd_Pakistan::ReadTable_BackgroundMortality_5yAgeGroup(string argStr, string argCountry)
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


int Project_DisBurd_Pakistan::ReadTable_Prob_Aware(string argStr)
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





int Project_DisBurd_Pakistan::ReadTable_BackgroundMortality(string argStr)
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



int Project_DisBurd_Pakistan::Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap)
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
