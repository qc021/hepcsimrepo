#include "project_acute.h"
using namespace std;







int Project_Acute::CEA_BaseResults()
{
	// See exact correspondence of treatment regimen in trials.cpp: 
	// txProfileType GetTx(const modelParamType & argParam, const patientType & argPat)

	double pClearance_bak = _modelParam._transData._pr_acute_spont_clearance;



	_note = "base";
	Compare("F0_G1", "Acute_G1", 1);

	_note = "all_self_clear";
	_modelParam._transData._pr_acute_spont_clearance = 1;
	Compare("F0_G1", "Acute_G1", 1);

	_note = "none_self_clear_TxF0";
	_modelParam._transData._pr_acute_spont_clearance = 0;
	Compare("F0_G1", "Acute_G1", 1);

	_note = "none_self_clear_NoTxF0";
	OVERRIDE_NOTREATMENT = true;
	Compare("F0_G1", "Acute_G1", 1);


	OVERRIDE_NOTREATMENT = false;
	_modelParam._transData._pr_acute_spont_clearance = pClearance_bak;

	_note = "healthy_population";
	Compare("HealthyPopulation", "HealthyPopulation", 1);

	OVERRIDE_NOTREATMENT = false;
	_modelParam._transData._pr_acute_spont_clearance = pClearance_bak;


	return 0;
}


vector<double> Project_Acute::Compare(string strArm1, string strArm2, int argGenotype)
{
	vector<double> outVec;

	//_outFileSummary.open("project_acute_output_compare.txt",ios::app);
	//_outFileSummary<<endl<<"===== "<<currentDateTime()<<" ===== "<<strArm1<<" vs "<<strArm2<<" ==="<<endl
	//	<< fixed << showpoint;
	// print the table header:

	_outFileSummary.open(_outFile, ios::app);
	_outFileSummary << fixed << showpoint;


	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/

	char listGender[] = { 'M','F' };



	for (int k = 0; k < sizeof(listGender) / sizeof(char); k++) {
		double wtScenario = _mfDistr[k];


		cout << "Comparing\t[" << strArm1 << "] - [" << strArm2 << "]\tAge=" << _age_initial
			<< "\tGender=" << listGender[k] << "\tGenotype=" << argGenotype << "\t";
		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]

		baseCohortType testCohort(s_Acute, _age_initial, listGender[k], 'W', argGenotype, 'N', 'A');

		//baseCohortType testCohort(s_F0, _age_initial, listGender[k], 'W', argGenotype, 'N', 'A');
		//strArm1 = "NoTx"; strArm2 = "TN_G1_LDP_SOF12";


		vector<double> r1 = EvaluateOneArm(strArm1, testCohort);
		vector<double> r2 = EvaluateOneArm(strArm2, testCohort);

		cout << fixed
			<< setprecision(3) << r1[0] << " (" << setprecision(0) << r1[1] << ")" << "\t"
			<< setprecision(3) << r2[0] << " (" << setprecision(0) << r2[1] << ")" << "\t";
		cout << endl;

		assert(r1.size() == r2.size());
		//_outFileSummary<<TIME_HORIZON<<"\t"<<_age_initial<<"\t"<<r1[0]<<"\t"<<r2[0]<<"\t"<<r2[0]-r1[0]<<"\t"
		//	<<r1[1]<<"\t"<<r2[1]<<"\t"<<r2[1]-r1[1]<<"\t"
		//	<<(r2[1]-r1[1])/(r2[0]-r1[0])<<"\t";

		//for(int k=2;k<r1.size(); k++){
		//	_outFileSummary<<r1[k]<<"\t"<<r2[k]<<"\t";
		//}
		//_outFileSummary<<endl;




		// initialize outVec
		if (k == 0) {
			for (int l = 0; l < 2 * r1.size(); l++) {
				outVec.push_back(0.0);
			}
		}
		for (int l = 0; l < r1.size(); l++) {
			outVec[l] = outVec[l] + wtScenario * r1[l];
			outVec[r1.size() + l] = outVec[r1.size() + l] + wtScenario * r2[l];
		}

	}

	_outFileSummary << NUM_PATIENTS << "\t" << SIM_SEED << "\t" << _note << "\t" << TIME_HORIZON << "\t" << _age_initial << "\t" << strArm1 << "\t";
	for (int l = 0; l < outVec.size() / 2; l++) { _outFileSummary << outVec[l] << "\t"; }
	_outFileSummary << endl;
	_outFileSummary << NUM_PATIENTS << "\t" << SIM_SEED << "\t" << _note << "\t" << TIME_HORIZON << "\t" << _age_initial << "\t" << strArm2 << "\t";
	for (int l = 0; l < outVec.size() / 2; l++) { _outFileSummary << outVec[l + outVec.size() / 2] << "\t"; }
	_outFileSummary << endl;


	//_outFileSummary << endl;
	_outFileSummary.close();
	return outVec;
}

vector<double> Project_Acute::EvaluateOneArm(string strArm, const baseCohortType & testCohort)
{
	modelParamType arm = _modelParam;
	arm._armName = strArm;
	arm._cohortData = testCohort;
	arm.ReduceDrugCost(_drugCostReduction);
	arm.SetAcuteTxDuration(_nWks_acute_tx);
	arm.SetF0TxDuration(_nWks_f0_tx);

	HepCSim * mySim;
	mySim = new HepCSim;
	mySim->SetRandomSeed(SIM_SEED);

	if (strArm == "HealthyPopulation") {// added 4/25/2017 for healthy population only.
		mySim->Run_HealthyPopulation(arm);
	}
	else {
		mySim->Run(arm);
	}
	// output one line summary of this arm
	vector<double> r;
	r.push_back(mySim->GetAvgQALY());
	r.push_back(mySim->GetAvgCost());
	r.push_back(mySim->GetCounter().countDeCirr);
	r.push_back(mySim->GetCounter().countHCC);
	r.push_back(mySim->GetCounter().countLivTr);
	r.push_back(mySim->GetCounter().countDeathLiv);
	r.push_back(mySim->GetTxCost());
	r.push_back(mySim->GetAvgLY());

	FreeMem(mySim);
	return r;

}

vector<double> Project_Acute::EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort)
{

	if ("TN_G1_LDP_SOF" == txArm) {
		vector<double> r;
		vector<double> r_8week = EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
		vector<double> r_12week = EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
		assert(r_8week.size() == r_12week.size());
		for (int l = 0; l < r_8week.size(); l++) {
			r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1 - _modelParam._transData._p8Week_LDV)*r_12week[l]);
		}
		return r;
	}
	else {
		return EvaluateOneArm(txArm, testCohort);
	}

}

int Project_Acute::CEA_OneWay(int argCmpIdx, string argFileOnewayParam, int argIdxForBatch)
{

	int nBatch = 20;


	ONEWAYSA_OPTION_PROJ_ACUTE = 1;
	_drugCostReduction_bak = _drugCostReduction;

	ofstream outf_low, outf_high;
	outf_low.open("output_project_acute/SA_LOW_" + _listCmp[argCmpIdx] + "_"+ basicToStr(argIdxForBatch) + "_0505.txt");
	outf_high.open("output_project_acute/SA_HIGH_" + _listCmp[argCmpIdx] + "_" + basicToStr(argIdxForBatch) + "_0505.txt");
	outf_low << fixed << showpoint;
	outf_high << fixed << showpoint;

	vector<string> sa_varName;
	vector<SParamTriplet> sa_varVal;
	ReadOneWaySAParamValueRange(argFileOnewayParam, sa_varName, sa_varVal);

	int sa_counter = 1;
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios = 2;
	stateType listState[] = { s_Acute,s_Acute };
	char listGender[] = { 'M','F' };
	double listAge[2];
	for (int k = 0; k < nScenarios; k++) { listAge[k] = _age_initial; }


	for (int k = 0; k < nScenarios; k++) {
		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', _listGenotype[argCmpIdx], 'N', 'A');

		// ===== baseline ======

		if (argIdxForBatch == 0) { // run base case only for the first batch

			_modelParam = _baselineModelParam;// recover the baseline values

			vector<double> r1_b = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_b = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

			assert(r1_b.size() == r2_b.size());
			outf_low << "base\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
				<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
				<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
				<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
			outf_high << "base\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
				<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
				<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
				<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
			for (int i = 2; i < r1_b.size(); i++) {
				outf_low << r1_b[i] << "\t" << r2_b[i] << "\t";
				outf_high << r1_b[i] << "\t" << r2_b[i] << "\t";
			}
			outf_low << endl;
			outf_high << endl;

			sa_counter++;
		}

		int batchSize = (int)sa_varName.size() / nBatch + 1;
		for (int j = 0; j < sa_varName.size(); j++) {

			if (j < argIdxForBatch * batchSize || j >= (argIdxForBatch + 1) * batchSize) {
				continue;
			}

			// ===== low value ======
			cout << endl << " - change value of paramter: " << sa_varName[j] << " = " << sa_varVal[j]._lb << "\t";
			_modelParam = _baselineModelParam;// recover the baseline values
			ChangeValue(sa_varName[j], sa_varVal[j]._lb, _modelParam);
			vector<double> r1_low = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_low = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

			outf_low << sa_varName[j] << "\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
				<< r1_low[0] << "\t" << r2_low[0] << "\t" << r2_low[0] - r1_low[0] << "\t"
				<< r1_low[1] << "\t" << r2_low[1] << "\t" << r2_low[1] - r1_low[1] << "\t"
				<< (r2_low[1] - r1_low[1]) / (r2_low[0] - r1_low[0]) << "\t";

			for (int i = 2; i < r1_low.size(); i++) {
				outf_low << r1_low[i] << "\t" << r2_low[i] << "\t";
			}
			outf_low << endl;

			// ===== high value ======
			cout << endl << " - change value of paramter: " << sa_varName[j] << " = " << sa_varVal[j]._ub << "\t";
			_modelParam = _baselineModelParam;// recover the baseline values
			ChangeValue(sa_varName[j], sa_varVal[j]._ub, _modelParam);
			vector<double> r1_high = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_high = EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

			outf_high << sa_varName[j] << "\t" << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
				<< r1_high[0] << "\t" << r2_high[0] << "\t" << r2_high[0] - r1_high[0] << "\t"
				<< r1_high[1] << "\t" << r2_high[1] << "\t" << r2_high[1] - r1_high[1] << "\t"
				<< (r2_high[1] - r1_high[1]) / (r2_high[0] - r1_high[0]) << "\t";
			for (int i = 2; i < r1_high.size(); i++) {
				outf_high << r1_high[i] << "\t" << r2_high[i] << "\t";
			}
			outf_high << endl;

			sa_counter++;
		}


		_drugCostReduction = _drugCostReduction_bak;

	}


	_modelParam = _baselineModelParam;//restore the nominal values
	ONEWAYSA_OPTION_PROJ_ACUTE = 0;


	outf_low.close(); outf_high.close();
	return 0;
}


int Project_Acute::CEA_PSA(int argIdx, int argNBatches) {
	ofstream outf_psa;
	outf_psa.open("output_project_acute/PSA_" + basicToStr(argIdx) + ".txt");
	outf_psa << fixed << showpoint;

	psaDistrType psaDistr;
	psaDistr.ReadDistrForPSA("project_acute_input_PSA_distr.txt");

	HepCSim mySim_psaSampler;		//used for sampling
	mySim_psaSampler.PSA_initializeSampler(psaDistr);

	int nIterOfThisBatch = NUM_RUNS_PSA / argNBatches;
	for (int n = 0; n < nIterOfThisBatch; n++) {
		if ((n + 1) % (nIterOfThisBatch / 10) == 0) {
			cout << (n + 1) << " PSA iterations finished...(" << (n + 1)*1.0 / (1.0 * nIterOfThisBatch) * 100 << "%)" << endl;
		}
		// sample
		mySim_psaSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
		// solve


		//if (argIdx == 0 && n == 0) {
		//	// print table head
		//	outf_psa << "iter\t"
		//		<< "QALY1\tCost1\t"
		//		<< "QALY2\tCost2\t"
		//		<< endl;
		//}

		vector<double> vecResult;
		outf_psa << argIdx * nIterOfThisBatch + n << "\t";
		vecResult = GetAggregatedResults();
		for (int k = 0; k < vecResult.size(); k++) {
			outf_psa << vecResult[k] << "\t";
		}
		outf_psa << endl;



	}
	outf_psa.close();

	return 0;


}



int Project_Acute::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
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

Project_Acute::Project_Acute()
{
	_note = "";
	_nWks_acute_tx = 6;
	_nWks_f0_tx = 8;
	_age_initial = 26.0; // base case average age.
	_drugCostReduction = 0;
	// read all arms for comparisons
	ifstream inf;
	inf.open("project_acute_input_comparators.txt");
	string word;
	while (inf >> word) {
		_listCmp.push_back(word);
		inf >> word; _listArm1.push_back(word);
		inf >> word; _listArm2.push_back(word);
		int g;
		inf >> g;	_listGenotype.push_back(g);
	}
	inf.close();

}






int Project_Acute::CEA_VaryAge(int argIdx)
{
	double age_bak = _age_initial;
	ofstream outf;
	outf.open("output_project_acute/out_change_age_0504.txt");
	outf << fixed << showpoint;

	double listInitialAge[] = { 20,25,30,35,40,45,50,60,70 };


	for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {

		SetInitialAge(listInitialAge[k]);
		vector<double> rs = GetAggregatedResults();

		outf << listInitialAge[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

	}
	outf.close();
	_age_initial = age_bak;
	return 0;
}


int Project_Acute::CEA_VaryTimeHorizon(int argIdx)
{
	double horz_bak = TIME_HORIZON;
	ofstream outf;
	outf.open("output_project_acute/out_change_horz_0504.txt");
	outf << fixed << showpoint;

	double listHorz[] = { 1,5,10,20,40, 150 };
	for (int k = 0; k < sizeof(listHorz) / sizeof(double); k++) {
		TIME_HORIZON = listHorz[k];
		vector<double> rs = GetAggregatedResults();
		outf << _age_initial << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;
	}

	outf.close();
	TIME_HORIZON = horz_bak;
	return 0;
}

vector<double> Project_Acute::GetAggregatedResults()
{
	return Compare("F0_G1", "Acute_G1", 1);
	// return Compare("F0_G1", "HealthyPopulation", 1);
}


int Project_Acute::CEA_VaryAcuteTxDuration_probSelfClearance(int argIdx)
{
	//SIM_SEED = 1;
	int dur_bak = _nWks_acute_tx;
	double pClr_bak = _modelParam._transData._pr_acute_spont_clearance;

	ofstream outf;
	outf.open("output_project_acute/out_change_acuteTxDur_prClr" + basicToStr(argIdx) + "_0417.txt");
	outf << fixed << showpoint;

	double listDur[] = { 4,5,6,7,8,9,10,11,12 };
	//int k = 2;
	//int k = argIdx;
	for (int k = 0; k < sizeof(listDur) / sizeof(double); k++) {



		//double listPrClr[] = { 0, 1, 0.5};
		double listPrClr[] = { 0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 };
		for (int m = 0; m < sizeof(listPrClr) / sizeof(double); m++) {
			double p = listPrClr[m];

			//for (double p = 0.1; p < 0.5; ) {

			_nWks_acute_tx = listDur[k];
			_modelParam._transData._pr_acute_spont_clearance = p;


			vector<double> rs = Compare("F0_G1", "Acute_G1_ChangableDuration", 1);

			outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << _nWks_acute_tx << "\t" << p << "\t";
			for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
			outf << endl;

		}
	}

	outf.close();

	_nWks_acute_tx = dur_bak;
	_modelParam._transData._pr_acute_spont_clearance = pClr_bak;

	return 0;
}



int Project_Acute::CEA_VaryAcuteTxDuration_SVR(int argIdx)
{


	//SIM_SEED = 1;
	int dur_bak = _nWks_acute_tx;

	ofstream outf;
	outf.open("output_project_acute/out_change_acuteTxDur_SVRacute" + basicToStr(argIdx) + "_0417.txt");
	outf << fixed << showpoint;

	double listDur[] = { 4,5,6,7,8,9,10,11,12 };
	//int k = 2;
	//int k = argIdx;
	for (int k = 0; k < sizeof(listDur) / sizeof(double); k++) {



		//double listPrClr[] = { 0, 1, 0.5};
		double listSVRacute[] = { 0.85,0.88,0.9,0.92,0.94,0.96,0.98 };
		for (int m = 0; m < sizeof(listSVRacute) / sizeof(double); m++) {
			double p = listSVRacute[m];

			ONEWAYSA_OPTION_PROJ_ACUTE = 1; // need this to change SVR rates for acute HCV treatment temporarily

			_nWks_acute_tx = listDur[k];
			_modelParam._transData._pr_svr_for_SA_acute = p;

			vector<double> rs = Compare("F0_G1", "Acute_G1_ChangableDuration", 1);

			outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << _nWks_acute_tx << "\t" << p << "\t";
			for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
			outf << endl;

		}
	}

	outf.close();

	_nWks_acute_tx = dur_bak;
	ONEWAYSA_OPTION_PROJ_ACUTE == 0;

	return 0;
}



int Project_Acute::CEA_VaryTxDuration_Acute_ChronicF0(int argIdx)
{


	//SIM_SEED = 1;
	int dur_bak_acute = _nWks_acute_tx;
	int dur_bak_f0 = _nWks_f0_tx;
	ofstream outf;
	outf.open("output_project_acute/out_change_TxDur_Acute_ChronicF0_" + basicToStr(argIdx) + "_0504.txt");
	outf << fixed << showpoint;

	double listDurAcute[] = { 6, 4 };
	double listDurChronic[] = { 8, 12, 24 };

	for (int k = 0; k < sizeof(listDurAcute) / sizeof(double); k++) {

		for (int m = 0; m < sizeof(listDurChronic) / sizeof(double); m++) {

			_nWks_acute_tx = listDurAcute[k];
			_nWks_f0_tx = listDurChronic[m];

			vector<double> rs = Compare("F0_G1_ChangableDuration", "Acute_G1_ChangableDuration", 1);

			outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << _nWks_f0_tx << "\t" << _nWks_acute_tx << "\t";
			for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
			outf << endl;

		}
	}

	outf.close();

	_nWks_acute_tx = dur_bak_acute;
	_nWks_f0_tx = dur_bak_f0;


	return 0;
}



int Project_Acute::CEA_Vary_TxCost_SVRacute(int argIdx)
{

	ofstream outf;
	outf.open("output_project_acute/out_change_TxCost_SVRacute" + basicToStr(argIdx) + "_0417.txt");
	outf << fixed << showpoint;

	double listCostDiscount[] = { -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5 };
	double listSVRacute[] = { 0.85,0.88,0.9,0.92,0.94,0.96,0.98, 0.99, 1.0 };
	for (int k = 0; k < sizeof(listCostDiscount) / sizeof(double); k++) {
		for (int m = 0; m < sizeof(listSVRacute) / sizeof(double); m++) {
			double p = listSVRacute[m];

			ONEWAYSA_OPTION_PROJ_ACUTE = 1; // need this to change SVR rates for acute HCV treatment temporarily

			SetDrugCostReduction(listCostDiscount[k]);
			_modelParam._transData._pr_svr_for_SA_acute = p;

			vector<double> rs = Compare("F0_G1", "Acute_G1", 1);

			outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << listCostDiscount[k] << "\t" << p << "\t";
			for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
			outf << endl;

		}
	}

	outf.close();

	ONEWAYSA_OPTION_PROJ_ACUTE == 0;
	SetDrugCostReduction(0);

	return 0;
}



int Project_Acute::CEA_Vary_TxCost_ClearanceRate(int argIdx)
{
	double pClr_bak = _modelParam._transData._pr_acute_spont_clearance;

	ofstream outf;
	outf.open("output_project_acute/out_change_TxCost_ClrRate" + basicToStr(argIdx) + "_0417.txt");
	outf << fixed << showpoint;

	double listCostDiscount[] = { -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5 };
	double listPrClr[] = { 0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 };
	for (int k = 0; k < sizeof(listCostDiscount) / sizeof(double); k++) {

		for (int m = 0; m < sizeof(listPrClr) / sizeof(double); m++) {
			double p = listPrClr[m];

			SetDrugCostReduction(listCostDiscount[k]);
			_modelParam._transData._pr_acute_spont_clearance = p;


			vector<double> rs = Compare("F0_G1", "Acute_G1", 1);

			outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << listCostDiscount[k] << "\t" << p << "\t";
			for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
			outf << endl;

		}
	}

	outf.close();


	SetDrugCostReduction(0);
	_modelParam._transData._pr_acute_spont_clearance = pClr_bak;

	return 0;
}



int Project_Acute::CEA_Vary_TxCost_SVRacute_ClrRate(int argIdx)
{
	double pClr_bak = _modelParam._transData._pr_acute_spont_clearance;
	ofstream outf;
	outf.open("output_project_acute/out_change_TxCost_SVRacute_ClrRate_" + basicToStr(argIdx) + "_0504.txt");
	outf << fixed << showpoint;

	double listCostDiscount[] = { 0,-0.2,-0.4 };
	double listSVRacute[] = { 0.85,0.88,0.9,0.92,0.94,0.96,0.98, 0.99, 1.0 };
	double listPrClr[] = { 0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 };

	for (int k = 0; k < sizeof(listCostDiscount) / sizeof(double); k++) {
		for (int m = 0; m < sizeof(listSVRacute) / sizeof(double); m++) {
			int l = argIdx; // 0, ..., 10
			//for (int l = 0; l < sizeof(listPrClr) / sizeof(double); l++) {
				SetDrugCostReduction(listCostDiscount[k]);

				ONEWAYSA_OPTION_PROJ_ACUTE = 1; // need this to change SVR rates for acute HCV treatment temporarily
				_modelParam._transData._pr_svr_for_SA_acute = listSVRacute[m];

				_modelParam._transData._pr_acute_spont_clearance = listPrClr[l];


				vector<double> rs = Compare("F0_G1", "Acute_G1", 1);

				outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << listCostDiscount[k] << "\t" << listSVRacute[m] << "\t" << listPrClr[l] << "\t";
				for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
				outf << endl;
			//}
		}
	}

	outf.close();

	ONEWAYSA_OPTION_PROJ_ACUTE == 0;
	SetDrugCostReduction(0);
	_modelParam._transData._pr_acute_spont_clearance = pClr_bak;

	return 0;
}

int Project_Acute::CEA_Vary_AccessProbability(int argIdx)
{

	double prob_bak = _modelParam._transData._pr_chronic_prob_follow_up;
	//SIM_SEED = 1;

	ofstream outf;
	outf.open("output_project_acute/out_change_accessProb_chronicArm_0504.txt");
	outf << fixed << showpoint;

	double listAccessProb[] = {1.0, 0.95,0.9,0.85,0.8,0.7, 0.6,0.5,0.4,0.3, 0.2,0.1, 0.01,0 };
	//double listAccessProb[] = { 1.0, 0.8};
	for (int k = 0; k < sizeof(listAccessProb) / sizeof(double); k++) {
		_modelParam._transData._pr_acute_prob_complete_tx = listAccessProb[k];
		_modelParam._transData._pr_chronic_prob_follow_up = listAccessProb[k];

		vector<double> rs = Compare("F0_G1", "Acute_G1", 1);
		//vector<double> rs = Compare("Acute_G1","F0_G1", 1);
		outf << SIM_SEED << "\t" << _age_initial << "\t" << TIME_HORIZON << "\t" << listAccessProb[k] << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;
	}

	outf.close();

	_modelParam._transData._pr_chronic_prob_follow_up = prob_bak;
	return 0;
}
