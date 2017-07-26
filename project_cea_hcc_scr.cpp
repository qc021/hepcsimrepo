#include "project_cea_hcc_scr.h"
using namespace std;





int Project_CEA_HCC_Screening::ReadTable_InputParam()
{
	ifstream inFile;
	string strFile = "project_cea_hcc_screening_input_params.txt";
	inFile.open(strFile.c_str());		//US life tables of 2007, published in 2011
	if (inFile.fail()) {
		ExitWithMsg("Fail to open file " + strFile);
	}

	string str;
	
	while (inFile >> str) {
		if ("//" == str) {
			getline(inFile, str);
			continue;
		}
		if ("distr_tx" == str) {
			string condTx;
			inFile >> condTx;
			getline(inFile, str);
			istringstream iss(str);
			int ky;
			iss >> ky;
			vector<double> vecDistr;
			double val;
			while (iss >> val) {
				vecDistr.push_back(val);
			}
			_param_map_txDistr[condTx] = vecDistr;
		}
		else if ("distr_initial_fib" == str) {
			inFile >> _param_map_distrInitialFib[scr_fib_F3];
			inFile >> _param_map_distrInitialFib[scr_fib_F3_SVR];
			inFile >> _param_map_distrInitialFib[scr_fib_F4];
			inFile >> _param_map_distrInitialFib[scr_fib_F4_SVR];
			inFile >> _param_map_distrInitialFib[scr_fib_DC];

		}
		else {
			inFile >> _param_map[str];
		}
		
	}
	inFile.close();
	return 0;
}

int Project_CEA_HCC_Screening::ReadTable_BackgroundMortality(string argStr)
{
	ifstream inFile2;
	inFile2.open(argStr);		//US life tables of 2007, published in 2011
	if (inFile2.fail()) {
		ExitWithMsg("Fail to open file " + argStr);
	}

	int k;
	while (inFile2 >> k) {
		//cout << "k= " << k << endl;
		inFile2 >> _modelParam._hccScreeningData._bgMort_male[k]
			>> _modelParam._hccScreeningData._bgMort_female[k]
			>> _modelParam._hccScreeningData._qol_male[k]
			>> _modelParam._hccScreeningData._qol_female[k];
	}
	inFile2.close();
	return 0;


}

double Project_CEA_HCC_Screening::GetValue(string argVarName)
{
	map<string, double > ::const_iterator it = _param_map.find(argVarName);
	if (it == _param_map.end()) {
		ExitWithMsg("Error @ Project_CEA_HCC_Screening::GetValue(string argVarName): Can't find parameter: " + argVarName);
	}
	return it->second;
}

vector<double> Project_CEA_HCC_Screening::EvaluateScreeningPolicy(int argScrInterval)
{
	modelParamType theParam = _modelParam;
	cout << "==== Evaluating screening policy - every " << argScrInterval<< "\tmonths ("<<_note<<")"<<endl;
	HCCScreenSim * mySim;
	mySim = new HCCScreenSim;
	mySim->SetRandomSeed(SIM_SEED);
	mySim->SetScreeningInterval(argScrInterval);	
	mySim->Run(theParam);

	vector<double> r;
	r.push_back(mySim->GetAvgLY());
	r.push_back(mySim->GetAvgQALY());
	r.push_back(mySim->GetAvgCost());
	r.push_back(mySim->GetCounter().countDeCirr);
	r.push_back(mySim->GetCounter().countHCC);
	r.push_back(mySim->GetCounter().countLivTr);
	r.push_back(mySim->GetCounter().countDeathLiv);
	r.push_back(mySim->GetTestCost());
	

	FreeMem(mySim);

	_outFileSummary.open(_outFile, ios::app);
	_outFileSummary <<fixed<<setprecision(4) << _note << "\t" << argScrInterval << "\t";
	for (int k = 0; k < (int)r.size(); k++) {
		_outFileSummary << r[k] << "\t";
	}
	_outFileSummary << endl;
	_outFileSummary.close();
	return r;
}

int Project_CEA_HCC_Screening::SetInitialPopulation(type_scr_fib_state argFibState)
{
	for (map<type_scr_fib_state, double>::iterator it = _modelParam._hccScreeningData._distr_initial_fib.begin();
		it != _modelParam._hccScreeningData._distr_initial_fib.end();
		it++) {
		if (it->first == argFibState) {
			it->second = 1.0;
		}
		else {
			it->second = 0.0;
		}
	}
	return 0;
}




int Project_CEA_HCC_Screening::Initialize()
{
	cout << "==== Initialize parameters " << endl;
	ReadTable_InputParam();
	ReadTable_BackgroundMortality("Input_mortality_male_female_2011.in");
	_modelParam._hccScreeningData._distr_male = GetValue("distr_male");
	_modelParam._hccScreeningData._initial_age = GetValue("initial_age");
	
	_modelParam._hccScreeningData._distr_genotype[1] = GetValue("distr_genotype_1");
	_modelParam._hccScreeningData._distr_genotype[2] = GetValue("distr_genotype_2");
	_modelParam._hccScreeningData._distr_genotype[3] = GetValue("distr_genotype_3");
	_modelParam._hccScreeningData._distr_genotype[456] = GetValue("distr_genotype_456");


	// ---------- HCC treatment parameters: qol, cost, and survival outcomes -------------------
	_modelParam._qolData.q_tx_transpl_1y = GetValue("qol_transpl_1y");
	_modelParam._qolData.q_tx_transpl_1yPlus = GetValue("qol_transpl_1y+");
	_modelParam._qolData.qChange_tx_res_1y = GetValue("qol_change_resection_1y");
	_modelParam._qolData.qChange_tx_res_1yPlus = GetValue("qol_change_resection_1y+");
	_modelParam._qolData.qChange_tx_abl_1y = GetValue("qol_change_ablation_1y");
	_modelParam._qolData.qChange_tx_abl_1yPlus = GetValue("qol_change_ablation_1y+");
	_modelParam._qolData.qChange_tx_pal_1y = GetValue("qol_change_pal_1y");
	_modelParam._qolData.qChange_tx_pal_1yPlus = GetValue("qol_change_pal_1y+");
	
	_modelParam._costData.hccscr_c_tx_transpl_1y = GetValue("cost_transpl_1y");
	_modelParam._costData.hccscr_c_tx_transpl_1yPlus = GetValue("cost_transpl_1y+");
	_modelParam._costData.hccscr_c_tx_res = GetValue("cost_resection");
	_modelParam._costData.hccscr_c_tx_abl = GetValue("cost_ablation");
	_modelParam._costData.hccscr_c_tx_pal = GetValue("cost_palliative");


	// 1-year and 5-year OS input data
	_modelParam._transData._mort_tx_transpl_1y = 1.0-GetValue("OS_transpl_1y");
	_modelParam._transData._mort_tx_transpl_1yPlus = funcConvertProbForCycle(1.0 - GetValue("OS_transpl_5y")/GetValue("OS_transpl_1y"), 4.0);
	//  5-year OS input data
	_modelParam._transData._mort_tx_res_small = funcConvertProbForCycle(1.0 - GetValue("OS_resection_small"), 5.0);
	_modelParam._transData._mort_tx_res_med = funcConvertProbForCycle(1.0 - GetValue("OS_resection_med"), 5.0);
	_modelParam._transData._mort_tx_res_large = funcConvertProbForCycle(1.0 - GetValue("OS_resection_large"), 5.0);
	_modelParam._transData._mort_tx_abl_cc_small = funcConvertProbForCycle(1.0 - GetValue("OS_ablation_cc_small"), 5.0);
	_modelParam._transData._mort_tx_abl_cc_med = funcConvertProbForCycle(1.0 - GetValue("OS_ablation_cc_med"), 5.0);
	_modelParam._transData._mort_tx_abl_cc_large = funcConvertProbForCycle(1.0 - GetValue("OS_ablation_cc_large"), 5.0);
	_modelParam._transData._mort_tx_abl_dc_small = funcConvertProbForCycle(1.0 - GetValue("OS_ablation_dc_small"), 5.0);
	_modelParam._transData._mort_tx_abl_dc_med = funcConvertProbForCycle(1.0 - GetValue("OS_ablation_dc_med"), 5.0);
	_modelParam._transData._mort_tx_abl_dc_large = funcConvertProbForCycle(1.0 - GetValue("OS_ablation_dc_large"), 5.0);
	// 2-year OS input data
	_modelParam._transData._mort_tx_pal = funcConvertProbForCycle(1.0 - GetValue("OS_pal"), 2.0);


	_modelParam._hccScreeningData._distr_hcc_tx = _param_map_txDistr;
	_modelParam._hccScreeningData._list_hcc_tx.push_back(scr_abs_res);
	_modelParam._hccScreeningData._list_hcc_tx.push_back(scr_abs_trspl);
	_modelParam._hccScreeningData._list_hcc_tx.push_back(scr_abs_abl);
	_modelParam._hccScreeningData._list_hcc_tx.push_back(scr_abs_pal);

	// -------- HCC screening testing --------------
	_modelParam._costData.hccscr_c_test_surveillance = GetValue("cost_surveillance");
	_modelParam._costData.hccscr_c_test_diagnostic = GetValue("cost_diagnostic");
	_modelParam._transData._sens_surveillance = GetValue("sens_surveillance");
	_modelParam._transData._spec_surveillance = GetValue("spec_surveillance");
	_modelParam._transData._sens_diagnostic = GetValue("sens_diagnostic");

	// ------ use some old estimates -----
	_modelParam._costData.c_F3 = 1394;
	_modelParam._costData.c_CoCirr = 1626;
	_modelParam._costData.c_DeCirr = 18064;

	_modelParam._transData._pr_transpl_hcv = funcConvertProbForCycle(GetValue("prob_transpl_hcv"),3.0);
	_modelParam._transData._pr_transpl_hcc = funcConvertProbForCycle(GetValue("prob_transpl_hcc"), 3.0);

	_modelParam._hccScreeningData._distr_initial_fib = _param_map_distrInitialFib;

	return 0;
}

int Project_CEA_HCC_Screening::CEA_BaseResults()
{

	//EvaluateScreeningPolicy(scr_fib_F3, 10000); // no screening
	//EvaluateScreeningPolicy(scr_fib_F3, 20000); // no screening
	//EvaluateScreeningPolicy(scr_fib_F3, 3);
	//EvaluateScreeningPolicy(scr_fib_F3, 6);
	//EvaluateScreeningPolicy(scr_fib_F3, 12);
	//EvaluateScreeningPolicy(scr_fib_F3, 24);
	//EvaluateScreeningPolicy(scr_fib_F3, 36);


	_note = "mixed";
	SetInitialPopulationMixed();
	EvaluateScreeningPolicy(10000);// no screening
	EvaluateScreeningPolicy(36);
	EvaluateScreeningPolicy(24);
	EvaluateScreeningPolicy(12);
	EvaluateScreeningPolicy(6);
	EvaluateScreeningPolicy(3);

	_note = "F4_cohort";
	SetInitialPopulation(scr_fib_F4);
	EvaluateScreeningPolicy(10000);// no screening
	EvaluateScreeningPolicy(36);
	EvaluateScreeningPolicy(24);
	EvaluateScreeningPolicy(12);
	EvaluateScreeningPolicy(6);
	EvaluateScreeningPolicy(3);


	_note = "F3_cohort";
	SetInitialPopulation(scr_fib_F3);
	EvaluateScreeningPolicy(10000);// no screening
	EvaluateScreeningPolicy(36);
	EvaluateScreeningPolicy(24);
	EvaluateScreeningPolicy(12);
	EvaluateScreeningPolicy(6);
	EvaluateScreeningPolicy(3); 

	return 0;
}



int Project_CEA_HCC_Screening::CEA_OneWay(int argCmpIdx, string argFileOnewayParam, int argIdxForBatch)
{

	//int nBatch = 20;


	//ONEWAYSA_OPTION_PROJ_ACUTE = 1;
	//_drugCostReduction_bak = _drugCostReduction;

	//ofstream outf_low, outf_high;
	//outf_low.open("output_project_acute/SA_LOW_" + _listCmp[argCmpIdx] + "_"+ basicToStr(argIdxForBatch) + "_0505.txt");
	//outf_high.open("output_project_acute/SA_HIGH_" + _listCmp[argCmpIdx] + "_" + basicToStr(argIdxForBatch) + "_0505.txt");
	//outf_low << fixed << showpoint;
	//outf_high << fixed << showpoint;

	//vector<string> sa_varName;
	//vector<SParamTriplet> sa_varVal;
	//ReadOneWaySAParamValueRange(argFileOnewayParam, sa_varName, sa_varVal);

	//int sa_counter = 1;
	//// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	///********************************** Define trial arms **************************************************************/
	//// === scenarios ===
	//int nScenarios = 2;
	//stateType listState[] = { s_Acute,s_Acute };
	//char listGender[] = { 'M','F' };
	//double listAge[2];
	//for (int k = 0; k < nScenarios; k++) { listAge[k] = _age_initial; }


	//for (int k = 0; k < nScenarios; k++) {
	//	// define cohort
	//	// create patient profile:
	//	// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
	//	baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', _listGenotype[argCmpIdx], 'N', 'A');

	//	// ===== baseline ======

	//	if (argIdxForBatch == 0) { // run base case only for the first batch

	//		_modelParam = _baselineModelParam;// recover the baseline values

	//		vector<double> r1_b = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
	//		vector<double> r2_b = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

	//		assert(r1_b.size() == r2_b.size());
	//		outf_low << "base\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
	//			<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
	//			<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
	//			<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
	//		outf_high << "base\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
	//			<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
	//			<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
	//			<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
	//		for (int i = 2; i < r1_b.size(); i++) {
	//			outf_low << r1_b[i] << "\t" << r2_b[i] << "\t";
	//			outf_high << r1_b[i] << "\t" << r2_b[i] << "\t";
	//		}
	//		outf_low << endl;
	//		outf_high << endl;

	//		sa_counter++;
	//	}

	//	int batchSize = (int)sa_varName.size() / nBatch + 1;
	//	for (int j = 0; j < sa_varName.size(); j++) {

	//		if (j < argIdxForBatch * batchSize || j >= (argIdxForBatch + 1) * batchSize) {
	//			continue;
	//		}

	//		// ===== low value ======
	//		cout << endl << " - change value of paramter: " << sa_varName[j] << " = " << sa_varVal[j]._lb << "\t";
	//		_modelParam = _baselineModelParam;// recover the baseline values
	//		ChangeValue(sa_varName[j], sa_varVal[j]._lb, _modelParam);
	//		vector<double> r1_low = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
	//		vector<double> r2_low = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

	//		outf_low << sa_varName[j] << "\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
	//			<< r1_low[0] << "\t" << r2_low[0] << "\t" << r2_low[0] - r1_low[0] << "\t"
	//			<< r1_low[1] << "\t" << r2_low[1] << "\t" << r2_low[1] - r1_low[1] << "\t"
	//			<< (r2_low[1] - r1_low[1]) / (r2_low[0] - r1_low[0]) << "\t";

	//		for (int i = 2; i < r1_low.size(); i++) {
	//			outf_low << r1_low[i] << "\t" << r2_low[i] << "\t";
	//		}
	//		outf_low << endl;

	//		// ===== high value ======
	//		cout << endl << " - change value of paramter: " << sa_varName[j] << " = " << sa_varVal[j]._ub << "\t";
	//		_modelParam = _baselineModelParam;// recover the baseline values
	//		ChangeValue(sa_varName[j], sa_varVal[j]._ub, _modelParam);
	//		vector<double> r1_high = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
	//		vector<double> r2_high = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

	//		outf_high << sa_varName[j] << "\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
	//			<< r1_high[0] << "\t" << r2_high[0] << "\t" << r2_high[0] - r1_high[0] << "\t"
	//			<< r1_high[1] << "\t" << r2_high[1] << "\t" << r2_high[1] - r1_high[1] << "\t"
	//			<< (r2_high[1] - r1_high[1]) / (r2_high[0] - r1_high[0]) << "\t";
	//		for (int i = 2; i < r1_high.size(); i++) {
	//			outf_high << r1_high[i] << "\t" << r2_high[i] << "\t";
	//		}
	//		outf_high << endl;

	//		sa_counter++;
	//	}


	//	_drugCostReduction = _drugCostReduction_bak;

	//}


	//_modelParam = _baselineModelParam;//restore the nominal values
	//ONEWAYSA_OPTION_PROJ_ACUTE = 0;


	//outf_low.close(); outf_high.close();
	return 0;
}


int Project_CEA_HCC_Screening ::CEA_PSA(int argIdx, int argNBatches) {
	//ofstream outf_psa;
	//outf_psa.open("output_project_acute/PSA_" + basicToStr(argIdx) + ".txt");
	//outf_psa << fixed << showpoint;

	//psaDistrType psaDistr;
	//psaDistr.ReadDistrForPSA("project_acute_input_PSA_distr.txt");

	//HepCSim mySim_psaSampler;		//used for sampling
	//mySim_psaSampler.PSA_initializeSampler(psaDistr);

	//int nIterOfThisBatch = NUM_RUNS_PSA / argNBatches;
	//for (int n = 0; n < nIterOfThisBatch; n++) {
	//	if ((n + 1) % (nIterOfThisBatch / 10) == 0) {
	//		cout << (n + 1) << " PSA iterations finished...(" << (n + 1)*1.0 / (1.0 * nIterOfThisBatch) * 100 << "%)" << endl;
	//	}
	//	// sample
	//	mySim_psaSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
	//	// solve


	//	//if (argIdx == 0 && n == 0) {
	//	//	// print table head
	//	//	outf_psa << "iter\t"
	//	//		<< "QALY1\tCost1\t"
	//	//		<< "QALY2\tCost2\t"
	//	//		<< endl;
	//	//}

	//	vector<double> vecResult;
	//	outf_psa << argIdx * nIterOfThisBatch + n << "\t";
	//	vecResult = GetAggregatedResults();
	//	for (int k = 0; k < vecResult.size(); k++) {
	//		outf_psa << vecResult[k] << "\t";
	//	}
	//	outf_psa << endl;



	//}
	//outf_psa.close();

	return 0;


}


int Project_CEA_HCC_Screening::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
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
	else if ("c_drug_discount" == varName) {
		_drugCostReduction = argVal;
	}
	else {
		ExitWithMsg("[Error] Project_Acute::ChangeValue(string varName, double argVal): Unknown parameter " + varName);
	}
	return 0;
}

Project_CEA_HCC_Screening::Project_CEA_HCC_Screening()
{
	_note = "";
	_drugCostReduction = 0;
	_param_map_distrInitialFib[scr_fib_F3] = 0;
	_param_map_distrInitialFib[scr_fib_F3_SVR] = 0;
	_param_map_distrInitialFib[scr_fib_F4] = 0;
	_param_map_distrInitialFib[scr_fib_F4_SVR] = 0;
	_param_map_distrInitialFib[scr_fib_DC] = 0;
}





