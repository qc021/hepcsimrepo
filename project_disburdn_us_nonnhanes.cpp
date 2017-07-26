#include "project_disburdn_us_nonnhanes.h"
using namespace std;

Project_DisBurd_UScomprh::Project_DisBurd_UScomprh()
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



}


int Project_DisBurd_UScomprh::Initialize_UScomprh_Parameter()
{
	ReadTable_Param();


	// <genotype-int, probability>, use updated estimates from He2016 (Nainin 2006)
	_modelParam._disBurdnData._table_distr_genotype[1] = 0.796;
	_modelParam._disBurdnData._table_distr_genotype[2] = 0.13;
	_modelParam._disBurdnData._table_distr_genotype[3] = 0.063;
	_modelParam._disBurdnData._table_distr_genotype[456] = 0.011;


	double p_aware_in_new_incidence = 0.072; // 7% CDC data
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[40] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[50] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[60] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_uninsured[200] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[40] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[50] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[60] = p_aware_in_new_incidence;
	_modelParam._disBurdnData._table_aware_prob_among_new_incidence_insured[200] = p_aware_in_new_incidence;

	_modelParam._disBurdnData._flag_multipleCohort = true;

	// ****** immigrant parameters ********************************************************************
	_modelParam._disBurdnData._flag_include_immigrants_lpr = true;
	ReadTable_immigrants_hcv_cases("project_disburdn_uscomprh_input_immigrants.txt");
	_modelParam._nonNHANESPplData._immigrants_age_category.clear();
	
	ReadTable_immigrant_age_distr("project_disburdn_uscomprh_input_immigrants_age_distr.txt",_modelParam._nonNHANESPplData._immigrants_age_category, _modelParam._nonNHANESPplData._immigrant_distr_age);
	_modelParam._nonNHANESPplData._immigrants_male_pct = 0.458;
	_modelParam._nonNHANESPplData._immigrants_year_start = 2010;
	// ****** NS5A market share: simplified with 0 ****************************************************
	ReadTable_BackgroundMortality("Input_mortality_male_female_burdn.in");
	ReadTable_BackgroundMortality_5yAgeGroup_IndianResv("Input_mortality_male_female_american_indian.in");

	ReadTable_TreatmentCapacity("project_disburdn_txfailure_input_tx_cap.txt");
	ReadTable_NS5AMarketShare("project_disburdn_txfailure_input_ns5a_marketshare.txt");
	
	//ReadTable_Incidence_Num("project_disburdn_txfailure_input_incidence_updated2015.txt");
	ReadTable_Incidence_Num("project_disburdn_txfailure_input_incidence_updated2015_linearGrowth.txt");
		
	ReadTable_Incidence_AgeDistr("project_disburdn_txfailure_input_inci_age_distr_byCDCacuteReport.txt");

	ReadTable_Prob_Aware("project_disburdn_txfailure_input_probAware.txt");


	Initialize_SVR_Table("project_disburdn_txfailure_input_svr.txt", _modelParam._disBurdnData._table_SVR);
	_modelParam._disBurdnData._prob_tx_coverage_by_medicarePartD = 1.0; // 0.9;


	double p_uninsured = _modelParam._disBurdnData._table_distr_insurance[insr_uninsured];
	double p_private = _modelParam._disBurdnData._table_distr_insurance[insr_private];
	double p_medicaid = _modelParam._disBurdnData._table_distr_insurance[insr_medicaid];
	double p_mil = _modelParam._disBurdnData._table_distr_insurance[insr_military];
	for (int y = ACA_START_YR; y <= END_YR; y++) {
		map<typeInsurance, double> insr_distr;
		// TO-DO: adjust the insurance distribution

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




	return 0;
}


int Project_DisBurd_UScomprh::ReadTable_Param()
{

	string strFile = "project_disburdn_uscomprh_input_by_cohort.txt";
	ifstream inFile;
	inFile.open(strFile.c_str());		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + strFile);
	}

	// read table header, register countries
	vector<typeCohortDisBurdnModel> listCohort;
	string t;
	getline(inFile, t);
	istringstream iss(t);
	string word;
	while (iss >> word) {		
		listCohort.push_back(_modelParam._nonNHANESPplData.GetCohortType(word));
	}

	_modelParam._nonNHANESPplData._mapNonNHANESPplData.clear();
	_modelParam._nonNHANESPplData._distr_age_initial_category_by_cohort.clear();
	_modelParam._nonNHANESPplData._distr_age_initial_by_cohort.clear();
	_modelParam._nonNHANESPplData._distr_age_new_inci_category_by_cohort.clear();
	_modelParam._nonNHANESPplData._distr_age_new_inci_by_cohort.clear();

	string str, str2;
	while (inFile >> str) {	
		if ("//" == str) {
			getline(inFile, str);
			continue;
		}
		else if ("pplAdded_distrAge" == str) {
			inFile >> str;
			getline(inFile, str2);
			// first line is age category
			vector<int> vecCat = ReadRowAsIntVector(inFile);
			// second line is age distribution
			vector<double> vecDist = ReadRowAsDoubleVector(inFile);
			if (vecCat.size() != vecDist.size()) {
				ExitWithMsg("Error @ ReadParam(): Unmatched vectors for age distribution (cateogory and distribution)");
			}
			_modelParam._nonNHANESPplData._distr_age_initial_category_by_cohort[_modelParam._nonNHANESPplData.GetCohortType(str)] = vecCat;
			_modelParam._nonNHANESPplData._distr_age_initial_by_cohort[_modelParam._nonNHANESPplData.GetCohortType(str)] = vecDist;

		}
		else if ("new_inci_distrAge" == str) {
			inFile >> str;
			getline(inFile, str2);
			// first line is age category
			vector<int> vecCat = ReadRowAsIntVector(inFile);
			// second line is age distribution
			vector<double> vecDist = ReadRowAsDoubleVector(inFile);
			if (vecCat.size() != vecDist.size()) {
				ExitWithMsg("Error @ ReadParam(): Unmatched vectors for age distribution (cateogory and distribution)");
			}
			_modelParam._nonNHANESPplData._distr_age_new_inci_category_by_cohort[_modelParam._nonNHANESPplData.GetCohortType(str)] = vecCat;
			_modelParam._nonNHANESPplData._distr_age_new_inci_by_cohort[_modelParam._nonNHANESPplData.GetCohortType(str)] = vecDist;
		}
		else if ("annual_TxCap" == str) {
			inFile >> str;
			map<typeCohortDisBurdnModel, int> row;
			for (int k = 0; k < (int)listCohort.size(); k++) {
				inFile >> row[listCohort[k]];
			}
			_modelParam._nonNHANESPplData._txCap_by_year_cohort[atoi(str.c_str())] = row;
		}
		else if ("annual_scr_rate" == str) {
			inFile >> str;
			map<typeCohortDisBurdnModel, double> row;
			for (int k = 0; k < (int)listCohort.size(); k++) {
				inFile >> row[listCohort[k]];
			}
			_modelParam._nonNHANESPplData._scrRate_by_year_cohort[atoi(str.c_str())] = row;

		}
		else {
			map<typeCohortDisBurdnModel, double> row;
			for (int k = 0; k < (int)listCohort.size(); k++) {
				inFile >> row[listCohort[k]];
			}
			_modelParam._nonNHANESPplData._mapNonNHANESPplData[str] = row;

		}
	}
	inFile.close();


	for (map<string, typeCohortDisBurdnModel>::const_iterator it = _modelParam._nonNHANESPplData._mapNonNHANESCohortName.begin();
		it != _modelParam._nonNHANESPplData._mapNonNHANESCohortName.end(); it++) {
		map<int, double> distr_genotype;
		distr_genotype[1] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrGenotype_G1", it->second);
		distr_genotype[2] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrGenotype_G2", it->second);
		distr_genotype[3] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrGenotype_G3", it->second);
		distr_genotype[456] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrGenotype_G456", it->second);
		_modelParam._nonNHANESPplData._distr_genotype_by_cohort[it->second] = distr_genotype;

		map<typeInsurance, double> distr_insr;
		distr_insr[insr_uninsured] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("PplAdded_distrInsurance_uninsured", it->second);
		distr_insr[insr_private] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("PplAdded_distrInsurance_private", it->second);
		distr_insr[insr_medicaid] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("PplAdded_distrInsurance_medicaid", it->second);
		distr_insr[insr_military] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("PplAdded_distrInsurance_military", it->second);
		_modelParam._nonNHANESPplData._distr_insur_by_cohort[it->second] = distr_insr;

		map<typeStateDisBurdnModel, double> distr_state;
		distr_state[state_F0] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_F0", it->second);
		distr_state[state_F1] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_F1", it->second);
		distr_state[state_F2] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_F2", it->second);
		distr_state[state_F3] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_F3", it->second);
		distr_state[state_CoCirr] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_F4", it->second);
		distr_state[state_DeCirr] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_DC", it->second);
		distr_state[state_DeCirr1yrPlus] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_DC+", it->second);
		distr_state[state_LivTr] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_LivTr", it->second);
		distr_state[state_HCC] = _modelParam._nonNHANESPplData.GetNonNHANESPplData("pplAdded_distrHCVState_HCC", it->second);
		_modelParam._nonNHANESPplData._distr_fib_by_cohort[it->second] = distr_state;
	}

	return 0;
}

double Project_DisBurd_UScomprh::GetValue_byCohort(string argVarName, const  typeCohortDisBurdnModel & argCohort)
{
	map<string, map<typeCohortDisBurdnModel, double >>::const_iterator it = _modelParam._nonNHANESPplData._mapNonNHANESPplData.find(argVarName);
	if (it == _modelParam._nonNHANESPplData._mapNonNHANESPplData.end()) {
		ExitWithMsg("Error @ Project_DisBurd_nonNHANES::GetValue(string argVarName, typeCohortDisBurdnModel argCohort): Can't find parameter: "+argVarName);
	}

	map<typeCohortDisBurdnModel, double>::const_iterator it2 = it->second.find(argCohort);
	if (it2 == it->second.end()) {
		ExitWithMsg("Error @ Project_DisBurd_nonNHANES::GetValue(string argVarName, typeCohortDisBurdnModel argCohort): Can't find parameter: "
			+argVarName + " in cohort "+basicToStr((int)argCohort));
	}


	return it2->second;

}


int Project_DisBurd_UScomprh::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
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





int Project_DisBurd_UScomprh::BaseCaseRun()
{
	BurdenModelSim myModel;
	myModel.Run(_modelParam);

	//getchar();
	return 0;
}



int Project_DisBurd_UScomprh::ReadTable_BackgroundMortality_5yAgeGroup_IndianResv(string argStr)
{

	// ======================================
	// NOTE: life table is for 5-year age groups
	// need to use .lowerbound() function, 
	// instead of .find() function to retrive the background mortality
	// ======================================
	map<int, double> table_mFemale, table_mMale, table_qFemale, table_qMale;

	ifstream inFile;
	inFile.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	string str;
	double val;
	while (inFile >> str) {
		if (str == "//") {
			getline(inFile, str);
			continue;
		}
		inFile >> val;
		table_mFemale.insert(pair<int, double>(atoi(str.c_str()), val));
		inFile >> val;
		table_mMale.insert(pair<int, double>(atoi(str.c_str()), val));
		inFile >> val;
		table_qFemale.insert(pair<int, double>(atoi(str.c_str()), val));
		inFile >> val;
		table_qMale.insert(pair<int, double>(atoi(str.c_str()), val));
	}
	inFile.close();


	// fill-in n the table
	_modelParam._disBurdnData._bgMort_AI_male.clear();
	_modelParam._disBurdnData._bgMort_AI_female.clear();
	_modelParam._disBurdnData._qol_AI_male.clear();
	_modelParam._disBurdnData._qol_AI_female.clear();

	for (int age = 0; age < 120; age++) {
		double p_male = table_mMale.lower_bound(age)->second;
		double p_female = table_mFemale.lower_bound(age)->second;

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

		_modelParam._disBurdnData._bgMort_AI_male[age] = p_male;
		_modelParam._disBurdnData._bgMort_AI_female[age] = p_female;
		_modelParam._disBurdnData._qol_AI_male[age] = table_qMale.lower_bound(age)->second;
		_modelParam._disBurdnData._qol_AI_female[age] = table_qFemale.lower_bound(age)->second;
	}


	//_modelParam._disBurdnData._bgMort_male = table_mMale[argCountry];
	//_modelParam._disBurdnData._bgMort_female = table_mFemale[argCountry];
	//_modelParam._disBurdnData._qol_male = table_qMale[argCountry];
	//_modelParam._disBurdnData._qol_female = table_qFemale[argCountry];
	return 0;


}


int Project_DisBurd_UScomprh::ReadTable_Prob_Aware(string argStr)
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





int Project_DisBurd_UScomprh::ReadTable_BackgroundMortality(string argStr)
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



int Project_DisBurd_UScomprh::Initialize_SVR_Table(string argStr, map_svr_table & argSVRMap)
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


int Project_DisBurd_UScomprh::ReadTable_immigrants_hcv_cases(string argStr)
{
	ifstream inf;
	inf.open(argStr.c_str());
	if (inf.fail())      ExitWithMsg("ReadTable_2Col(): Fail to read " + argStr);
	typeStateDisBurdnModel listState[] = { state_F0,state_F1,state_F2,state_F3,state_CoCirr };
	int listGenotype[] = { 1,2,3,456 };
	string word;
	int num;
	while (inf >> word) {
		if (word == "//") {
			getline(inf, word);
			continue;
		}
		int yr = atoi(word.c_str());
		// skip the next three columns
		inf >> word;
		inf >> word;
		inf >> word;

		map<typeStateDisBurdnModel, map<int, int>> map_fib_gt;
		for (int k = 0; k<sizeof(listState)/sizeof(typeStateDisBurdnModel); k++) {
			map<int, int> map_g;
			for (int j = 0; j < sizeof(listGenotype) / sizeof(int); j++) {				
				double v;
				inf >> v;
				map_g.insert(pair<int, int>(listGenotype[j], (int)v));
			}
			map_fib_gt[listState[k]] = map_g;
		}
		_modelParam._nonNHANESPplData._immigrants_hcv_num_by_year_fib_genotype[yr] = map_fib_gt;
	}
	inf.close();
	return 0;
}

int Project_DisBurd_UScomprh::ReadTable_immigrant_age_distr(string argStr, vector<int> & argAgeCat, vector<double>& argAgeDistr)
{
	ifstream inf;
	inf.open(argStr.c_str());
	if (inf.fail())      ExitWithMsg("ReadTable_immigrant_age_distr(): Fail to read " + argStr);
	
	argAgeCat.clear();
	argAgeDistr.clear();

	string word;
	double distr;
	while (inf >> word) {
		if (word == "//") {
			getline(inf, word);
			continue;
		}
		int ageCat = atoi(word.c_str());
		
		inf >> distr;
		argAgeCat.push_back(ageCat);
		argAgeDistr.push_back(distr);
	}
	inf.close();
	return 0; 
	
}

int Project_DisBurd_UScomprh::ReadTable_TreatmentCapacity(string argStr)
{
	ReadTable_2Col(argStr, _modelParam._disBurdnData._table_tx_capacity);
	return 0;

}
int Project_DisBurd_UScomprh::ReadTable_Incidence_Num(string argStr)
{
	ReadTable_2Col(argStr, _modelParam._disBurdnData._table_incidence);
	return 0;
}


int Project_DisBurd_UScomprh::ReadTable_Incidence_AgeDistr(string argStr)
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
int Project_DisBurd_UScomprh::ReadTable_NS5AMarketShare(string argStr)
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
