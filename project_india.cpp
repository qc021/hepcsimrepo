#include "project_india.h"
using namespace std;




int Project_India::CEA_BaseResults()
{
	////TN TOL
	////Compare("TN_TOL_G1_BOC","TN_TOL_G1_SOF_PEG_RBV12",1);
	////Compare("TN_TOL_G1_TEL","TN_TOL_G1_SOF_PEG_RBV12",1);

	////Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF8",1);
	////Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF12",1);
	//Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF",1);
	////Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF8",1);
	////Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF12",1);
	//Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF",1);
	//Compare("TN_TOL_G2_PEG_RBV24","TN_TOL_G2_SOF_RBV12",2);
	//Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_RBV24",3);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_PEG_RBV12",4);

	//// TN NOT-TOL
	////Compare("NoTx","TN_NOT_G1_SOF_SMV12",1);

	////Compare("NoTx","TN_G1_LDP_SOF8",1);
	////Compare("NoTx","TN_G1_LDP_SOF12",1);
	//Compare("NoTx","TN_G1_LDP_SOF",1);
	//Compare("NoTx","TN_NOT_G2_SOF_RBV12",2);
	//Compare("NoTx","TN_NOT_G3_SOF_RBV24",3);
	//Compare("NoTx","TN_NOT_G4_SOF_RBV24",4);

	//// TE
	////Compare("TE_G1_BOC","TE_G1_SOF_SMV12",1);
	////Compare("TE_G1_TEL","TE_G1_SOF_SMV12",1);

	//Compare("TE_G1_BOC","TE_G1_LDP_SOF12_24",1);
	//Compare("TE_G1_TEL","TE_G1_LDP_SOF12_24",1);
	//Compare("TE_G2_PEG_RBV24","TE_G2_SOF_RBV12",2);
	//Compare("TE_G3_PEG_RBV24","TE_G3_SOF_RBV24",3);
	//Compare("TE_G4_PEG_RBV48","TE_G4_SOF_PEG_RBV12",4);


	// ========================================================================
	// All arms used in the calculator
	// Updated on Aug 27, 2015
	// genotype 11 = 1a
	// genotype 12 = 1b
	// - 1a and 1b share the same treatment effectiveness of TN_G1_LDP_SOF and TE_G1_LDV_SOF12_24
	// ========================================================================

	Compare("NoTx", "G1_SOF_LDV12", 1);
	//Compare("G1_SOF_LDV12", "NoTx", 1);
	//Compare("NoTx","G1_SOF_DCV12_24",1);

	//Compare("NoTx","TN_G1_LDP_SOF12",11);	
	//Compare("NoTx","TN_TOL_G1a_DCV_SOF12_24",11);
	//Compare("NoTx","TN_TOL_G1a_PrOD_RBV12_24",11);

	//Compare("NoTx","TN_G1_LDP_SOF12",12);
	//Compare("NoTx","TN_TOL_G1b_DCV_SOF12_24",12);
	//Compare("NoTx","TN_TOL_G1b_PrOD12",12);
	//Compare("NoTx","TN_TOL_G1b_SOF_SMV12",12);

	Compare("NoTx", "TN_TOL_G2_SOF_RBV12_16", 2);


	Compare("NoTx", "G3_SOF_DCV12_24", 3);
	//Compare("NoTx","TN_TOL_G3_DCV_SOF12_24",3);
	//Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_RBV24",3);
	//Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_PEG_RBV12",3);

	Compare("NoTx", "G4_SOF_LDV12", 4);
	//Compare("NoTx","TN_TOL_G4_LDV_SOF12",4);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_PrOD_RBV12",4);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_RBV24",4);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_PEG_RBV12",4);


	////_strTxExperience="TN_INTOL";
	//Compare("NoTx","TN_G1_LDP_SOF",11);
	//Compare("NoTx","TN_TOL_G1a_DCV_SOF12_24",11);
	//Compare("NoTx","TN_TOL_G1a_PrOD_RBV12_24",11);

	//Compare("NoTx","TN_G1_LDP_SOF",12);
	//Compare("NoTx","TN_TOL_G1b_DCV_SOF12_24",12);
	//Compare("NoTx","TN_TOL_G1b_PrOD12",12);
	//Compare("NoTx","TN_TOL_G1b_SOF_SMV12",12);

	//Compare("NoTx","TN_TOL_G2_SOF_RBV12_16",2);

	//Compare("NoTx","TN_TOL_G3_SOF_RBV24",3);
	//Compare("NoTx","TN_TOL_G3_DCV_SOF12_24",3);
	//Compare("NoTx","TN_TOL_G3_SOF_PEG_RBV12",3);


	//Compare("NoTx","TN_TOL_G4_LDV_SOF12",4);
	//Compare("NoTx","TN_TOL_G4_PrOD_RBV12",4);
	//Compare("NoTx","TN_TOL_G4_SOF_RBV24",4);
	//Compare("NoTx","TN_TOL_G4_SOF_PEG_RBV12",4);

	return 0;
}

vector<double> Project_India::GetAggregatedResults() {
	vector<double> outVec;
	vector<double> v1, v2, v3, v4;
	v1 = Compare("NoTx", "G1_SOF_LDV12", 1);
	v2 = Compare("NoTx", "TN_TOL_G2_SOF_RBV12_16", 2);
	v3 = Compare("NoTx", "G3_SOF_DCV12_24", 3);
	v4 = Compare("NoTx", "G4_SOF_LDV12", 4);
	for (int k = 0; k < v1.size(); k++) {
		outVec.push_back(_modelParam._pplData._distr_genotype["G1"] * v1[k]
			+ _modelParam._pplData._distr_genotype["G2"] * v2[k]
			+ _modelParam._pplData._distr_genotype["G3"] * v3[k]
			+ _modelParam._pplData._distr_genotype["G4"] * v4[k]);
	}	
	return outVec;
}

vector<double> Project_India::Compare(string strArm1, string strArm2, int argGenotype)
{

	vector<double> outVec;
	_outFileSummary.open(_outFile, ios::app);
	//_outFileSummary.open("project_india_output_compare.txt",ios::app);
	//_outFileSummary<<endl<<"===== "<<currentDateTime()<<" ===== "<<strArm1<<" vs "<<strArm2<<" ==="<<endl
	//	<< fixed << showpoint;

	_outFileSummary << fixed << showpoint;

	// print the table header:



	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios = 10;
	stateType listState[] = { s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr };
	char listGender[] = { 'M','F','M','F','M','F','M','F','M','F' };

	//double listAge[] = {50,50,56,56,58,58,58,58,59,59};
	//double listAge[] = { 25.9,25.9,		37.9,37.9,		42.8,42.8,		37.7,37.7,		48.4,48.4 };
	//double listAge[]={35,35,35,35,35,35,35,35,35,35};
	double listAge[10];
	for (int k = 0; k < sizeof(listAge) / sizeof(double); k++) {
		listAge[k] = _age_initial;
	}

	//int nScenarios=30;
	//stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr,
	//	s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr,
	//	s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	//char listGender[]={'M','F','M','F','M','F','M','F','M','F',
	//	'M','F','M','F','M','F','M','F','M','F',
	//	'M','F','M','F','M','F','M','F','M','F'};
	//double listAge[]={35,35,35,35,35,35,35,35,35,35,
	//	20,20,20,20,20,20,20,20,20,20,		
	//	60,60,60,60,60,60,60,60,60,60};

	_modelParam.ApplyIndiaParameters();
	for (int k = 0; k < nScenarios; k++) {
		double wtScenario = _modelParam._pplData._distr_fib[listState[k]] * _modelParam._pplData._distr_gender[listGender[k]];


		cout << "Comparing\t" << strArm1 << "-" << strArm2 << ":\tFibrosis=" << listState[k] << "\tAge=" << listAge[k]
			<< "\tGender=" << listGender[k] << "\tGenotype=" << argGenotype << "\t" << endl;

		int genotype_noSubtype = (argGenotype == 11 || argGenotype == 12) ? 1 : argGenotype;

		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		//baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', argGenotype, 'N', 'A');	
		baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', genotype_noSubtype, 'N', 'A');

		//vector<double> r1=EvaluateOneArm(strArm1, testCohort);
		//vector<double> r2=EvaluateOneArm(strArm2, testCohort);

		// changed, 7/9/2016
		//vector<double> r1=EvaluateOneArm_CombineLDV8Or12WkOnline(strArm1, testCohort);
		//vector<double> r2=EvaluateOneArm_CombineLDV8Or12WkOnline(strArm2, testCohort);
		vector<double> r1 = EvaluateOneArm_Averaged(strArm1, testCohort);
		vector<double> r2 = EvaluateOneArm_Averaged(strArm2, testCohort);

		assert(r1.size() == r2.size());
		_outFileSummary << _age_initial << "\t" << TIME_HORIZON << "\t"
			<< r1[0] << "\t" << r2[0] << "\t" << r2[0] - r1[0] << "\t"
			<< r1[1] << "\t" << r2[1] << "\t" << r2[1] - r1[1] << "\t"
			<< (r2[1] - r1[1]) / (r2[0] - r1[0]) << "\t";

		for (int l = 2; l < r1.size(); l++) {
			_outFileSummary << r1[l] << "\t" << r2[l] << "\t";
		}

		_outFileSummary << _modelParam._costData.c_SOF << "\t";
		_outFileSummary << endl;


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

	_outFileSummary << endl;
	_outFileSummary.close();



	return outVec;
}


vector<double> Project_India::EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort)
{
	vector<double> r;
	if ("TN_TOL_G1_BOC_TEL" == txArm) {
		vector<double> r_boc = EvaluateOneArm("TN_TOL_G1_BOC", testCohort);
		vector<double> r_tel = EvaluateOneArm("TN_TOL_G1_TEL", testCohort);
		r = LinearCombTwoVec(r_boc, r_tel, _modelParam._transData._pBOC);

	}
	else if ("TE_G1_BOC_TEL" == txArm) {
		vector<double> r_boc = EvaluateOneArm("TE_G1_BOC", testCohort);
		vector<double> r_tel = EvaluateOneArm("TE_G1_TEL", testCohort);
		r = LinearCombTwoVec(r_boc, r_tel, _modelParam._transData._pBOC);

	}
	else if ("TN_G1_LDP_SOF" == txArm) {
		vector<double> r_8week = EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
		vector<double> r_12week = EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
		r = LinearCombTwoVec(r_8week, r_12week, _modelParam._transData._p8Week_LDV);
		//assert(r_8week.size()==r_12week.size());
		//for(int l=0; l<r_8week.size(); l++){
		//	r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1-_modelParam._transData._p8Week_LDV)*r_12week[l]);
		//}
	}
	else {
		r = EvaluateOneArm(txArm, testCohort);
	}
	return r;
}



vector<double> Project_India::EvaluateOneArm(string strArm, const baseCohortType & testCohort)
{
	modelParamType arm = _modelParam;
	arm._armName = strArm;
	arm._cohortData = testCohort;
	arm.ReduceDrugCost(_drugCostReduction);
	//arm.ApplyIndiaParameters();

	HepCSim * mySim;
	mySim = new HepCSim;
	mySim->SetRandomSeed(SIM_SEED);
	mySim->Run(arm);
	// output one line summary of this arm
	vector<double> r;
	r.push_back(mySim->GetAvgQALY());
	r.push_back(mySim->GetAvgCost());
	r.push_back(mySim->GetCounter().countDeCirr);
	r.push_back(mySim->GetCounter().countHCC);
	r.push_back(mySim->GetCounter().countLivTr);
	r.push_back(mySim->GetCounter().countDeathLiv);
	r.push_back(mySim->GetTxCost());
	r.push_back(mySim->GetDALY_YLL());
	r.push_back(mySim->GetDALY_YLD());
	r.push_back(mySim->GetDALY_Total());
	r.push_back(mySim->GetTestCost());

	FreeMem(mySim);
	return r;

}

//vector<double> Project_India::EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort )
//{
//
//	if("TN_G1_LDP_SOF" == txArm){
//		vector<double> r;
//		vector<double> r_8week=EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
//		vector<double> r_12week=EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
//		assert(r_8week.size()==r_12week.size());
//		for(int l=0; l<r_8week.size(); l++){
//			r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1-_modelParam._transData._p8Week_LDV)*r_12week[l]);
//		}
//		return r;
//	}else{
//		return EvaluateOneArm(txArm, testCohort);
//	}
//
//}

int Project_India::CEA_OneWay(int argCmpIdx, string argFileOnewayParam)
{



	ofstream outf_low, outf_high;
	outf_low.open("output_project_india/SA_LOW_" + _listCmp[argCmpIdx] + ".txt");
	outf_high.open("output_project_india/SA_HIGH_" + _listCmp[argCmpIdx] + ".txt");
	outf_low << fixed << showpoint;
	outf_high << fixed << showpoint;

	vector<string> sa_varName;
	vector<SParamTriplet> sa_varVal;
	ReadOneWaySAParamValueRange(argFileOnewayParam, sa_varName, sa_varVal);

	int sa_counter = 1;
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios = 10;
	stateType listState[] = { s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr };
	char listGender[] = { 'M','F','M','F','M','F','M','F','M','F' };
	double listAge[] = { 35,35,35,35,35,35,35,35,35,35 };


	_baselineModelParam.ApplyIndiaParameters();

	for (int k = 0; k < nScenarios; k++) {
		cout << "Comparing\t" << _listArm1[argCmpIdx] << "-" << _listArm2[argCmpIdx] << ":\tFibrosis=" << listState[k] << "\tAge=" << listAge[k]
			<< "\tGender=" << listGender[k] << "\tGenotype=" << _listGenotype[argCmpIdx] << "\t" << endl;

		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', _listGenotype[argCmpIdx], 'N', 'A');

		// ===== baseline ======		
		_modelParam = _baselineModelParam;// recover the baseline values

		vector<double> r1_b = EvaluateOneArm_Averaged(_listArm1[argCmpIdx], testCohort);
		vector<double> r2_b = EvaluateOneArm_Averaged(_listArm2[argCmpIdx], testCohort);

		assert(r1_b.size() == r2_b.size());
		outf_low << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
			<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
			<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
			<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
		outf_high << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
			<< r1_b[0] << "\t" << r2_b[0] << "\t" << r2_b[0] - r1_b[0] << "\t"
			<< r1_b[1] << "\t" << r2_b[1] << "\t" << r2_b[1] - r1_b[1] << "\t"
			<< (r2_b[1] - r1_b[1]) / (r2_b[0] - r1_b[0]) << "\t";
		for (int i = 2; i < r1_b.size(); i++) {
			outf_low << r1_b[i] << "\t" << r2_b[i] << "\t";
			outf_high << r1_b[i] << "\t" << r2_b[i] << "\t";
		}
		outf_low << r2_b[9] - r1_b[9] << "\t"; // incremental dalys
		outf_high << r2_b[9] - r1_b[9] << "\t"; // incremental dalys
		outf_low << endl;
		outf_high << endl;

		sa_counter++;

		for (int j = 0; j < sa_varName.size(); j++) {
			cout << " - sensitivity analysis on variable: " << sa_varName[j] << endl;
			// ===== low value ======
			_modelParam = _baselineModelParam;// recover the baseline values
			ChangeValue(sa_varName[j], sa_varVal[j]._lb, _modelParam);
			vector<double> r1_low = EvaluateOneArm_Averaged(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_low = EvaluateOneArm_Averaged(_listArm2[argCmpIdx], testCohort);

			outf_low << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
				<< r1_low[0] << "\t" << r2_low[0] << "\t" << r2_low[0] - r1_low[0] << "\t"
				<< r1_low[1] << "\t" << r2_low[1] << "\t" << r2_low[1] - r1_low[1] << "\t"
				<< (r2_low[1] - r1_low[1]) / (r2_low[0] - r1_low[0]) << "\t";

			for (int i = 2; i < r1_low.size(); i++) {
				outf_low << r1_low[i] << "\t" << r2_low[i] << "\t";
			}

			outf_low << r2_low[9] - r1_low[9] << "\t"; // incremental dalys
			outf_low << endl;

			// ===== high value ======
			_modelParam = _baselineModelParam;// recover the baseline values
			ChangeValue(sa_varName[j], sa_varVal[j]._ub, _modelParam);
			vector<double> r1_high = EvaluateOneArm_Averaged(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_high = EvaluateOneArm_Averaged(_listArm2[argCmpIdx], testCohort);

			outf_high << "F" << listState[k] << " " << listGender[k] << " " << listAge[k] << "\t" << sa_counter << "\t"
				<< r1_high[0] << "\t" << r2_high[0] << "\t" << r2_high[0] - r1_high[0] << "\t"
				<< r1_high[1] << "\t" << r2_high[1] << "\t" << r2_high[1] - r1_high[1] << "\t"
				<< (r2_high[1] - r1_high[1]) / (r2_high[0] - r1_high[0]) << "\t";
			for (int i = 2; i < r1_high.size(); i++) {
				outf_high << r1_high[i] << "\t" << r2_high[i] << "\t";
			}
			outf_high << r2_high[9] - r1_high[9] << "\t"; // incremental dalys
			outf_high << endl;

			sa_counter++;
		}

	}


	_modelParam = _baselineModelParam;//restore the nominal values



	outf_low.close(); outf_high.close();
	return 0;
}

//int Project_India::ReadPSASampledValues( ifstream & inf, modelParamType & argModelParam )
//{
//	//update the model parameters:
//	//Fibrosis score (F0-F4)	sex	age	
//	// q_F0	q_F1	q_F2	q_F3	q_CoCirr	q_DeCirr	q_HCC	q_LivTr	q_PostLivT	q_SVR	q_Anemia	q_TX_oSOC	q_Tx_DAA	
//	inf>>argModelParam._qolData.q_F0;
//	inf>>argModelParam._qolData.q_F1;
//	inf>>argModelParam._qolData.q_F2;
//	inf>>argModelParam._qolData.q_F3;
//	inf>>argModelParam._qolData.q_CoCirr;
//	inf>>argModelParam._qolData.q_DeCirr;
//	inf>>argModelParam._qolData.q_HCC;
//	inf>>argModelParam._qolData.q_LivTr;
//	inf>>argModelParam._qolData.q_PostLivTr;
//	inf>>argModelParam._qolData.q_SVR;
//	inf>>argModelParam._qolData.q_Dec_Anemia;
//	inf>>argModelParam._qolData.q_TX_oSOC;
//	inf>>argModelParam._qolData.q_TX_DAA;
//
//	//c_F0	c_F1	c_F2	c_F3	c_CoCirr	c_DeCirr	c_DeCirr1yrPlus	c_HCC	c_LivTr	c_PostLivTr	
//	inf>>argModelParam._costData.c_F0;
//	inf>>argModelParam._costData.c_F1;
//	inf>>argModelParam._costData.c_F2;
//	inf>>argModelParam._costData.c_F3;
//	inf>>argModelParam._costData.c_CoCirr;
//	inf>>argModelParam._costData.c_DeCirr;
//	inf>>argModelParam._costData.c_DeCirr1yrPlus;
//	inf>>argModelParam._costData.c_HCC;
//	inf>>argModelParam._costData.c_LivTr;
//	inf>>argModelParam._costData.c_PostLivTr;
//
//	//pF0_F1_SA	pF1_F2_SA	pF2_F3_SA	pF3_F4_SA	pF4_DC_SA	pF4_HCC_SA	pDC_HCC_SA	pDC_Liv_Transpl_SA	
//	inf>>argModelParam._transData.pr_F0_F1;
//	inf>>argModelParam._transData.pr_F1_F2;
//	inf>>argModelParam._transData.pr_F2_F3;
//	inf>>argModelParam._transData.pr_F3_CoCirr;
//	inf>>argModelParam._transData.pr_CoCirr_DeCirr;
//	inf>>argModelParam._transData.pr_CoCirr_HCC;
//	inf>>argModelParam._transData.pr_DeCirr_HCC;
//	inf>>argModelParam._transData.pr_DeCirr_LivTr;
//
//	//pMort_dc_cyc_1_SA	pMort_dc_cyc_2_SA	pHCC_Liv_Transpl_SA	pMort_hcc_cyc_SA	pMort_Liv_Transpl_SA	pMort_Post_Liv_Transpl_SA	
//	inf>>argModelParam._transData.pr_DeCirr_DeathLiv;
//	inf>>argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv;
//	inf>>argModelParam._transData.pr_HCC_LivTr;
//	inf>>argModelParam._transData.pr_HCC_DeathLiv;
//	inf>>argModelParam._transData.pr_LivTr_DeathLiv;
//	inf>>argModelParam._transData.pr_LivTr1yrPlus_DeathLiv;
//
//	//pr_SVR_CoCirr_DeCirr	pr_SVR_CoCirr_HCC	pSVR_Delta_oSOC	pSVR_Delta_DAA
//	inf>>argModelParam._transData.pr_SVR_CoCirr_DeCirr;
//	inf>>argModelParam._transData.pr_SVR_CoCirr_HCC;
//	inf>>argModelParam._transData.pr_SVR_Delta_oSOC;
//	inf>>argModelParam._transData.pr_SVR_Delta_DAA;
//
//
//
//	return 0;
//}

int Project_India::CEA_PSA(int argIdx, int argNBatches) {
	ofstream outf_psa;
	outf_psa.open("output_project_india/PSA_" + basicToStr(argIdx) + ".txt");
	outf_psa << fixed << showpoint;

	psaDistrType psaDistr;
	psaDistr.ReadDistrForPSA("project_india_intput_PSA_distr.txt");

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
		outf_psa << argIdx * nIterOfThisBatch + n<<"\t";
		vecResult = GetAggregatedResults();
		for (int k = 0; k < vecResult.size(); k++) {
			outf_psa << vecResult[k] << "\t";
		}
		outf_psa << endl;



	}
	outf_psa.close();

	return 0;


}

//	ofstream outf_psa;
//	outf_psa.open("output_project_doublecheck/PSA_"+_listCmp[argCmpIdx]+"_"+basicToStr(argBatch)+".txt");
//	outf_psa<< fixed << showpoint;
//
//	ifstream inf;
//	inf.open("project_doublecheck_input_PSA.txt");
//	string line;
//	getline(inf,line);//skip the firstline
//	
//	//skip the first argBatch*1000 lines
//	for(int i=0; i<argBatchSize*argBatch; i++) getline(inf,line);
//
//	int sa_counter=1+argBatchSize*argBatch;
//
//	int fib;
//	char gender;
//	double initialAge;
//	for(int i=0; i<argBatchSize; i++){
//		inf>>fib;
//		inf>>gender;
//		inf>>initialAge;
//		// create cohort
//		baseCohortType testCohort(fib,initialAge,gender, 'W', _listGenotype[argCmpIdx], 'N', 'A');	// first three charateristics (fib, age, gender) will be updated in the next line
//		// read sampled parameter values
//		ReadPSASampledValues(inf,_modelParam);
//	
//
//		vector<double> r1=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
//		vector<double> r2=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);
//
//		assert(r1.size()==r2.size());
//		outf_psa<<"F"<<fib<<" "<<gender<<" "<<initialAge<<"\t"<<sa_counter<<"\t"
//			<<r1[0]<<"\t"<<r2[0]<<"\t"<<r2[0]-r1[0]<<"\t"
//			<<r1[1]<<"\t"<<r2[1]<<"\t"<<r2[1]-r1[1]<<"\t"
//			<<(r2[1]-r1[1])/(r2[0]-r1[0])<<"\t";		
//		for(int i=2;i<r1.size(); i++){
//			outf_psa<<r1[i]<<"\t"<<r2[i]<<"\t";
//		}
//		outf_psa<<endl;
//
//		sa_counter++;
//	}
//
//	inf.close();
//	outf_psa.close();
//
//	return 0;
//}

int Project_India::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
	if ("q_F0" == varName) {
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

	else if ("dw_F0" == varName) {
		argModelParam._dalyData._dw_f0 = argVal;
	}
	else if ("dw_F1" == varName) {
		argModelParam._dalyData._dw_f1 = argVal;
	}
	else if ("dw_F2" == varName) {
		argModelParam._dalyData._dw_f2 = argVal;
	}
	else if ("dw_F3" == varName) {
		argModelParam._dalyData._dw_f3 = argVal;
	}
	else if ("dw_CoCirr" == varName) {
		argModelParam._dalyData._dw_CoCirr = argVal;
	}
	else if ("dw_DeCirr" == varName) {
		argModelParam._dalyData._dw_DeCirr = argVal;
		argModelParam._dalyData._dw_DeCirr1yrPlus = argVal;
	}
	else if ("dw_HCC" == varName) {
		argModelParam._dalyData._dw_HCC = argVal;
	}
	else if ("dw_LivTr" == varName) {
		argModelParam._dalyData._dw_LivTr = argVal;
	}
	else if ("dw_PostLivT" == varName) {
		argModelParam._dalyData._dw_LivTr1yrPlus = argVal;
	}

	else if ("c_testing_preTx" == varName) {
		argModelParam._costData.c_testing_preTx = argVal;
	}
	else if ("c_testing_postTx" == varName) {
		argModelParam._costData.c_testing_postTx = argVal;
	}
	else {
		ExitWithMsg("[Error] Project_Doublecheck::ChangeValue(string varName, double argVal): Unknown parameter " + varName);
	}
	return 0;
}

Project_India::Project_India()
{
	
	_age_initial = 35.0; // base case average age.
	_drugCostReduction = 0;
	// read all arms for comparisons
	ifstream inf;
	inf.open("project_india_input_comparators.txt");
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

int Project_India::CEA_VaryTimeHorizon(int argIdx)
{
	double listHorz[] = { 1,2,3,4, 5,10,15,20,25,30,35,40,50 };
	TIME_HORIZON = listHorz[argIdx];
	CEA_BaseResults();
	return 0;
}

int Project_India::CEA_VaryAgeAndTimeHorizon(int argIdx)
{
	ofstream outf;
	outf.open("output_project_india/out_change_age_and_horizon" + basicToStr(argIdx) + ".txt");
	outf << fixed << showpoint;

	double listInitialAge[] = {20,30,40,50,60,70};
	double listHorz[] = {1,2,3,4,5,10,15,20,25,30,35,40,50};
	TIME_HORIZON = listHorz[argIdx];

	for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {
		
		SetInitialAge(listInitialAge[k]);
		vector<double> rs = GetAggregatedResults();

		if (argIdx == 0 && k == 0) {
		}

		outf << listInitialAge[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;
		
	}

	outf.close();
	return 0;
}



int Project_India::CEA_VaryAge(int argIdx)
{

	//double listInitialAge[] = { 20,25,30,35,40,45,50,55,60,65,70 };
	//for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {

	//	SetInitialAge(listInitialAge[k]);

	//	CEA_BaseResults();
	//}

	//return 0;


	ofstream outf;
	outf.open("output_project_india/out_change_age.txt");
	outf << fixed << showpoint;

	double listInitialAge[] = { 20,25,30,35,40,45,50,55,60,65,70 };


	for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {

		SetInitialAge(listInitialAge[k]);
		vector<double> rs = GetAggregatedResults();



		outf << listInitialAge[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

	}

	outf.close();
	return 0;
}


int Project_India::CEA_VaryDrugCost(int argIdx)
{

	//double listDrugCost[] = { 2,5,8 };

	//for (int k = 0; k < sizeof(listDrugCost) / sizeof(double); k++) {

	//	SetDrugCostReduction(listDrugCost[k]);

	//	CEA_BaseResults();
	//}

	//return 0;


	ofstream outf;
	outf.open("output_project_india/out_change_drugcost.txt");
	outf << fixed << showpoint;

	double listDrugCost[] = {0,2,5,8};


	for (int k = 0; k < sizeof(listDrugCost) / sizeof(double); k++) {

		SetDrugCostReduction(listDrugCost[k]);

		
		vector<double> rs = GetAggregatedResults();



		outf <<	_modelParam._costData.c_SOF * (1.0 + listDrugCost[k]) << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

		_modelParam = _baselineModelParam; //restore the nominal values
	}

	outf.close();
	return 0;
}

int Project_India::CEA_SA_BackgroundCost()
{
	ofstream outf;
	outf.open("output_project_india/out_change_background_cost.txt");
	outf << fixed << showpoint;

	double listDrugCost[] = {0,1000,2000,3000,5000,7000,10000};


	for (int k = 0; k < sizeof(listDrugCost) / sizeof(double); k++) {

		_modelParam._costData.c_background = listDrugCost[k];
		vector<double> rs = GetAggregatedResults();

		outf <<	listDrugCost[k] << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

		_modelParam = _baselineModelParam; //restore the nominal values
	}

	outf.close();
	return 0;

}



int Project_India::CEA_R1_LiverRelatedMortality()
{
	ofstream outf;
	outf.open("output_project_india/out_R1_change_lrdMortality.txt");
	outf << fixed << showpoint;

	// list of relative risks (RR)
	double listRR[] = {1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0};


	for (int k = 0; k < sizeof(listRR) / sizeof(double); k++) {

		_modelParam._transData.pr_DeCirr_DeathLiv = 1.0 - pow(1.0 -_modelParam._transData.pr_DeCirr_DeathLiv, listRR[k]);
		_modelParam._transData.pr_DeCirr1yrPlus_DeathLiv = 1.0 - pow(1.0 - _modelParam._transData.pr_DeCirr1yrPlus_DeathLiv, listRR[k]);
		_modelParam._transData.pr_HCC_DeathLiv = 1.0 - pow(1.0 - _modelParam._transData.pr_HCC_DeathLiv, listRR[k]);
		//_modelParam._transData.pr_DeCirr_DeathLiv =  _modelParam._transData.pr_DeCirr_DeathLiv * listRR[k];
		//_modelParam._transData.pr_DeCirr1yrPlus_DeathLiv = _modelParam._transData.pr_DeCirr1yrPlus_DeathLiv * listRR[k];
		//_modelParam._transData.pr_HCC_DeathLiv = _modelParam._transData.pr_HCC_DeathLiv * listRR[k];

		vector<double> rs = GetAggregatedResults();

		outf <<	listRR[k] << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

		_modelParam = _baselineModelParam; //restore the nominal values
	}

	outf.close();
	return 0;

}



int Project_India::CEA_R1_LifeExpectancy_w_wo_HCV()
{

	FILE_background_mortality = "Input_mortality_male_female_India_correct.in";
	ofstream outf;
	outf.open("output_project_india/out_R1_lifeexpectancy.txt");
	outf << fixed << showpoint;

	double listInitialAge[] = {30, 60};

	DISCOUNT_FACT_C = 0; //annual discount factor for costs
	DISCOUNT_FACT_Q = 0; //anual discount factor for QALYs
	EXPECTED_LIFE_ONLY = 1;

	for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {

		SetInitialAge(listInitialAge[k]);
		SetOverridePerfectTreatment(true);
		vector<double> rs = GetAggregatedResults();



		outf << listInitialAge[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

		_modelParam = _baselineModelParam;
	}

	outf.close();

	DISCOUNT_FACT_C = 0.03; //annual discount factor for costs
	DISCOUNT_FACT_Q = 0.03; //anual discount factor for QALYs
	EXPECTED_LIFE_ONLY = 0;
	SetOverridePerfectTreatment(false);

	return 0;
}



int Project_India::CEA_R1_discount()
{

	
	ofstream outf;
	outf.open("output_project_india/out_R1_discount.txt");
	outf << fixed << showpoint;

	double listDiscountFactor[] = {0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2};



	for (int k = 0; k < sizeof(listDiscountFactor) / sizeof(double); k++) {

		DISCOUNT_FACT_C = listDiscountFactor[k]; //annual discount factor for costs
		DISCOUNT_FACT_Q = listDiscountFactor[k]; //anual discount factor for QALYs
				
		
		vector<double> rs = GetAggregatedResults();


		outf << listDiscountFactor[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

		_modelParam = _baselineModelParam;
	}

	outf.close();


	return 0;
}


int Project_India::CEA_R1_progRisk()
{
	SetOutFile("output_project_india/out_change_progRisk_allruns.txt");

	ofstream outf;
	outf.open("output_project_india/out_change_progRisk_summary.txt");
	outf << fixed << showpoint;

	double listHR[] = {0.6,0.8,1.0,1.2};


	for (int k = 0; k < sizeof(listHR) / sizeof(double); k++) {

		_modelParam._transData.pr_F0_F1 = 1.0 - pow(1.0 - _modelParam._transData.pr_F0_F1, listHR[k]);
		_modelParam._transData.pr_F1_F2 = 1.0 - pow(1.0 - _modelParam._transData.pr_F1_F2, listHR[k]);
		_modelParam._transData.pr_F2_F3 = 1.0 - pow(1.0 - _modelParam._transData.pr_F2_F3, listHR[k]);
		_modelParam._transData.pr_F3_CoCirr = 1.0 - pow(1.0 - _modelParam._transData.pr_F3_CoCirr, listHR[k]);

		vector<double> rs = GetAggregatedResults();



		outf << listHR[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

		_modelParam = _baselineModelParam; //restore the nominal values
	}

	outf.close();
	return 0;
}


int Project_India::CEA_R1_progRisk_and_Horizon()
{

	double listHR[] = { 0.6,0.8,1.0,1.2 };
	double listHorz[] = { 1,2,3,4,5,10,15};
	

	for (int k = 0; k < sizeof(listHR) / sizeof(double); k++) {
		
		_modelParam._transData.pr_F0_F1 = 1.0 - pow(1.0 - _modelParam._transData.pr_F0_F1, listHR[k]);
		_modelParam._transData.pr_F1_F2 = 1.0 - pow(1.0 - _modelParam._transData.pr_F1_F2, listHR[k]);
		_modelParam._transData.pr_F2_F3 = 1.0 - pow(1.0 - _modelParam._transData.pr_F2_F3, listHR[k]);
		_modelParam._transData.pr_F3_CoCirr = 1.0 - pow(1.0 - _modelParam._transData.pr_F3_CoCirr, listHR[k]);

		SetOutFile("output_project_india/out_change_hr_and_horizon_allruns_" + basicToStr(k) + ".txt");
		ofstream outf;
		outf.open("output_project_india/out_change_hr_and_horizon_" + basicToStr(k) + ".txt");
		outf << fixed << showpoint;

		for (int j = 0; j < sizeof(listHorz) / sizeof(double); j++) {
			TIME_HORIZON = listHorz[j];

			vector<double> rs = GetAggregatedResults();
			
			outf << listHR[k] << "\t" << TIME_HORIZON << "\t";
			for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
			outf << endl;

			
		}
		outf.close();


		_modelParam = _baselineModelParam; //restore the nominal values
	}

	return 0;
}
