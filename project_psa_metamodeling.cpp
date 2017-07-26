#include "project_psa_metamodeling.h"


int Project_PSA_MetaModeling::RunBaseCase()
{

	
	_outf.open("output_project_metamodel/project_psa_basecase.txt");
	if (_outf.fail()) {
		ExitWithMsg("Error @ RunBaseCase(): Fail to open file output_PSA_metamodel/project_psa_basecase.txt");
	}
	PrintHeader(_outf);
	cout << "["<< currentDateTime() <<"] Running base case ..." << endl;
	// ======== version 1 ===============================	
	//NUM_PATIENTS = 1;	
	//for (int iterMC = 0; iterMC < 10000; iterMC++) {
	//	_seedOfOneMCRun = SIM_SEED + iterMC;
	//	Compare(iterMC, 0);
	//}
	// ======== version 2 ===============================
	NUM_PATIENTS = 10000;
	_seedOfOneMCRun = SIM_SEED;
	Compare(0, 0);	
	// ==================================================

	_outf.close();
	return 0;
}


int Project_PSA_MetaModeling::Compare(int argIdxMC, int argIdxPSA)
{



	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios = 1;
	stateType listState[] = { s_F0 };
	char listGender[] = { 'M' };
	double listAge[] = { 55.0 };
	int listGenotype[] = { 1 };

	//stateType listState[] = { s_CoCirr };
	//char listGender[] = { 'M' };
	//double listAge[] = { 35.0 };
	//int listGenotype[] = { 1 };
	//_modelParam._transData.pr_DeCirr_LivTr = 0;
	//_modelParam._transData.pr_HCC_LivTr = 0;

	string listTxtxArm[] = { "TN_TOL_G1_BOC","TN_G1_LDP_SOF12" };
	int k = 0;

	// define cohort
	// create patient profile: for example [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
	baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', listGenotype[k], 'N', 'A');

	vector<double> r1 = EvaluateOneArm(listTxtxArm[0], testCohort);
	vector<double> r2 = EvaluateOneArm(listTxtxArm[1], testCohort);
	assert(r1.size() == r2.size());

	// r1 and r2 has the outcomes in the following order:
	// [0]=QALY, [1]=Cost, [2]=#DeCirr, [3]=#HCC, [4]=#LivTr, [5]=#DeathLiv, [6]=Treatment Cost		

	_outf << argIdxPSA << "\t" << argIdxMC << "\t" << r1[0] << "\t" << r2[0] << "\t" << r2[0] - r1[0] << "\t"
		<< r1[1] << "\t" << r2[1] << "\t" << r2[1] - r1[1] << "\t"
		<< (r2[1] - r1[1]) / (r2[0] - r1[0]) << "\t";

	for (int k = 2; k < r1.size(); k++) {
		_outf << r1[k] << "\t" << r2[k] << "\t";
	}
	


	return 0;
}

vector<double> Project_PSA_MetaModeling::EvaluateOneArm(string strArm, const baseCohortType & testCohort)
{
	modelParamType arm = _modelParam;
	arm._armName = strArm;
	arm._cohortData = testCohort;
	arm.ReduceDrugCost(_drugCostReduction);

	HepCSim * mySim;
	mySim = new HepCSim;
	mySim->SetRandomSeed(_seedOfOneMCRun);
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

	FreeMem(mySim);
	return r;

}

int Project_PSA_MetaModeling::RunPSA_Multivariate_OutputSingleMCRun(long argNumPSAIter, int argNumRplc)
{
	ReadDistrForPSA("SA_PSA_distr.in");
	cout << "[" << currentDateTime() << "] Running PSA ..." << endl;
	string strFile = "output_project_metamodel/project_psa_output_" + basicToStr(argNumPSAIter) + "PSA_" + basicToStr(argNumRplc) + "MC_SEED" + basicToStr(PSA_SEED) + ".txt";
	_outf.open(strFile);
	if (_outf.fail()) {
		ExitWithMsg("Error @ RunBaseCase(): Fail to open file "+strFile);
	}
	
	//PrintHeader(_outf);	
	

	HepCSim paramSampler;		//used for sampling
	paramSampler.PSA_initializeSampler(_distrPSA_varName, _distrPSA_distrName, _distrPSA_param1, _distrPSA_param2);

	int fileCounter=0;
	for(int iterPSA=0; iterPSA < argNumPSAIter; iterPSA++){
		clock_t start=clock();

		if (argNumPSAIter >= 10) {
			if ((iterPSA + 1) % (argNumPSAIter / 10) == 0) {
				cout << "[" << currentDateTime() << "] " << (iterPSA + 1) << " PSA iterations finished...(" << (iterPSA + 1) / (1.0 * argNumPSAIter) * 100 << "%)" << endl;
			}
		}
		
		// sampling the value
		paramSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
		

		for (int iterMC = 0; iterMC < argNumRplc; iterMC++) {
			_seedOfOneMCRun = SIM_SEED + iterMC +iterPSA;
			Compare(iterMC, iterPSA);

			if (_flagPrintParameterValues) {

				_outf << _modelParam._qolData.q_F0 << "\t"
					<< _modelParam._qolData.q_F1 << "\t"
					<< _modelParam._qolData.q_F2 << "\t"
					<< _modelParam._qolData.q_F3 << "\t"
					<< _modelParam._qolData.q_CoCirr << "\t"
					<< _modelParam._qolData.q_DeCirr << "\t"
					<< _modelParam._qolData.q_HCC << "\t"
					<< _modelParam._qolData.q_LivTr << "\t"
					<< _modelParam._qolData.q_PostLivTr << "\t"
					<< _modelParam._qolData.q_SVR << "\t"
					<< _modelParam._qolData.q_Dec_Anemia << "\t"
					<< _modelParam._qolData.q_TX_oSOC << "\t"
					<< _modelParam._qolData.q_TX_DAA << "\t"

					<< _modelParam._costData.c_F0 << "\t"
					<< _modelParam._costData.c_F1 << "\t"
					<< _modelParam._costData.c_F2 << "\t"
					<< _modelParam._costData.c_F3 << "\t"
					<< _modelParam._costData.c_CoCirr << "\t"
					<< _modelParam._costData.c_DeCirr << "\t"
					<< _modelParam._costData.c_DeCirr1yrPlus << "\t"
					<< _modelParam._costData.c_HCC << "\t"
					<< _modelParam._costData.c_LivTr << "\t"
					<< _modelParam._costData.c_PostLivTr << "\t"

					<< _modelParam._transData.pr_F0_F1 << "\t"
					<< _modelParam._transData.pr_F1_F2 << "\t"
					<< _modelParam._transData.pr_F2_F3 << "\t"
					<< _modelParam._transData.pr_F3_CoCirr << "\t"
					<< _modelParam._transData.pr_CoCirr_DeCirr << "\t"
					<< _modelParam._transData.pr_CoCirr_HCC << "\t"
					<< _modelParam._transData.pr_DeCirr_HCC << "\t"
					<< _modelParam._transData.pr_DeCirr_LivTr << "\t"
					<< _modelParam._transData.pr_DeCirr_DeathLiv << "\t"
					<< _modelParam._transData.pr_DeCirr1yrPlus_DeathLiv << "\t"
					<< _modelParam._transData.pr_HCC_LivTr << "\t"
					<< _modelParam._transData.pr_HCC_DeathLiv << "\t"
					<< _modelParam._transData.pr_LivTr_DeathLiv << "\t"
					<< _modelParam._transData.pr_LivTr1yrPlus_DeathLiv << "\t"
					<< _modelParam._transData.pr_SVR_CoCirr_DeCirr << "\t"
					<< _modelParam._transData.pr_SVR_CoCirr_HCC << "\t"
					<< _modelParam._transData.pr_SVR_Delta_oSOC << "\t"
					<< _modelParam._transData.pr_SVR_Delta_DAA;
			}

				_outf<<endl;
		}

		
	}

	_outf.close();
	return 0;
}



int Project_PSA_MetaModeling::ReadDistrForPSA(string argFileStr)
{
	cout<<"Reading parameter distributions for PSA ..."<<endl;
	_distrPSA_varName.clear();
	_distrPSA_distrName.clear();
	_distrPSA_param1.clear();
	_distrPSA_param2.clear();

	ifstream inf;
	inf.open(argFileStr);

	string varName, line;
	while(inf>>varName){
		if("//" == varName ) {
			getline(inf,line);
			continue;
		}
		string distrName;
		inf>>distrName;
		double p1,p2;
		inf>>p1;
		inf>>p2;
		_distrPSA_varName.push_back(varName);
		_distrPSA_distrName.push_back(distrName);
		_distrPSA_param1.push_back(p1);
		_distrPSA_param2.push_back(p2);
	}
	assert(_distrPSA_distrName.size() == _distrPSA_param1.size() && _distrPSA_param1.size() == _distrPSA_param2.size());

	inf.close();

	return 0;
}

void Project_PSA_MetaModeling::PrintHeader( ofstream & outf )
{
	// *************************************************************
	// MAKE SURE THE HEADER IS MATCHING THE ACTUAL SIMULATION OUTPUT
	// *************************************************************
	// r1 and r2 has the outcomes in the following order:
	// [0]=QALY, [1]=Cost, [2]=#DeCirr, [3]=#HCC, [4]=#LivTr, [5]=#DeathLiv, [6]=Treatment Cost		
	//	_outf << argIdxPSA << "\t" << argIdxMC << "\t" << r1[0] << "\t" << r2[0] << "\t" << r2[0] - r1[0] << "\t"
	//	<< r1[1] << "\t" << r2[1] << "\t" << r2[1] - r1[1] << "\t"
	//	<< (r2[1] - r1[1]) / (r2[0] - r1[0]) << "\t";
	//	for (int k = 2; k < r1.size(); k++) {
	//		_outf << r1[k] << "\t" << r2[k] << "\t";
	//	}

	outf<< "IdxPSA\tIdxMC\tQALY1\tQALY2\tDeltaQ\tCost1\tCost2\tDeltaC\tICER\t" ;

	outf << "DC1\tDC2\tHCC1\tHCC2\tLT1\tLT2\tDeathLiv1\tDeathLiv2\tTxCost1\tTxCost2\t" ;

	if (_flagPrintParameterValues) {
		outf<< "q_F0"<<"\t"
			<<"q_F1"<<"\t"
			<<"q_F2"<<"\t"
			<<"q_F3"<<"\t"
			<<"q_CoCirr"<<"\t"
			<<"q_DeCirr"<<"\t"
			<<"q_HCC"<<"\t"
			<<"q_LivTr"<<"\t"
			<<"q_PostLivT"<<"\t"
			<<"q_SVR"<<"\t"
			<<"q_Anemia"<<"\t"
			<<"q_TX_oSOC"<<"\t"
			<<"q_Tx_DAA"<<"\t"
			<<"c_F0"<<"\t"
			<<"c_F1"<<"\t"
			<<"c_F2"<<"\t"
			<<"c_F3"<<"\t"
			<<"c_CoCirr"<<"\t"
			<<"c_DeCirr"<<"\t"
			<<"c_DeCirr1yrPlus"<<"\t"
			<<"c_HCC"<<"\t"
			<<"c_LivTr"<<"\t"
			<<"c_PostLivTr"<<"\t"
			<<"pF0_F1_SA"<<"\t"
			<<"pF1_F2_SA"<<"\t"
			<<"pF2_F3_SA"<<"\t"
			<<"pF3_F4_SA"<<"\t"
			<<"pF4_DC_SA"<<"\t"
			<<"pF4_HCC_SA"<<"\t"
			<<"pDC_HCC_SA"<<"\t"
			<<"pDC_Liv_Transpl_SA"<<"\t"
			<<"pMort_dc_cyc_1_SA"<<"\t"
			<<"pMort_dc_cyc_2_SA"<<"\t"
			<<"pHCC_Liv_Transpl_SA"<<"\t"
			<<"pMort_hcc_cyc_SA"<<"\t"
			<<"pMort_Liv_Transpl_SA"<<"\t"
			<<"pMort_Post_Liv_Transpl_SA"<<"\t"
			<<"pr_SVR_CoCirr_DeCirr"<<"\t"
			<<"pr_SVR_CoCirr_HCC"<<"\t"
			<<"pSVR_Delta_oSOC"<<"\t"
			<<"pSVR_Delta_DAA"<<"\t";

	}
	outf << endl;
	outf<<fixed<<setprecision(3);
}

int Project_PSA_MetaModeling::RunParameterSmplingTest( int seed, int num )
{
	PSA_SEED=seed;

	ofstream outf;
	outf.open("project_psa_output_sampled_param_values_testset.txt");
	outf<<"baseAge"<<"\t"<<"baseState"<<"\t"<<"sex"<<"\t"
		<<"q_F0"<<"\t"
		<<"q_F1"<<"\t"
		<<"q_F2"<<"\t"
		<<"q_F3"<<"\t"
		<<"q_CoCirr"<<"\t"
		<<"q_DeCirr"<<"\t"
		<<"q_HCC"<<"\t"
		<<"q_LivTr"<<"\t"
		<<"q_PostLivT"<<"\t"
		<<"q_SVR"<<"\t"
		<<"q_Anemia"<<"\t"
		<<"q_TX_oSOC"<<"\t"
		<<"q_Tx_DAA"<<"\t"
		<<"c_F0"<<"\t"
		<<"c_F1"<<"\t"
		<<"c_F2"<<"\t"
		<<"c_F3"<<"\t"
		<<"c_CoCirr"<<"\t"
		<<"c_DeCirr"<<"\t"
		<<"c_DeCirr1yrPlus"<<"\t"
		<<"c_HCC"<<"\t"
		<<"c_LivTr"<<"\t"
		<<"c_PostLivTr"<<"\t"
		<<"pF0_F1_SA"<<"\t"
		<<"pF1_F2_SA"<<"\t"
		<<"pF2_F3_SA"<<"\t"
		<<"pF3_F4_SA"<<"\t"
		<<"pF4_DC_SA"<<"\t"
		<<"pF4_HCC_SA"<<"\t"
		<<"pDC_HCC_SA"<<"\t"
		<<"pDC_Liv_Transpl_SA"<<"\t"
		<<"pMort_dc_cyc_1_SA"<<"\t"
		<<"pMort_dc_cyc_2_SA"<<"\t"
		<<"pHCC_Liv_Transpl_SA"<<"\t"
		<<"pMort_hcc_cyc_SA"<<"\t"
		<<"pMort_Liv_Transpl_SA"<<"\t"
		<<"pMort_Post_Liv_Transpl_SA"<<"\t"
		<<"pr_SVR_CoCirr_DeCirr"<<"\t"
		<<"pr_SVR_CoCirr_HCC"<<"\t"
		<<"pSVR_Delta_oSOC"<<"\t"
		<<"pSVR_Delta_DAA"<<"\t"
		<<endl;

	ReadDistrForPSA("SA_PSA_distr.in");
	HepCSim mySim;		//used for sampling
	mySim.PSA_initializeSampler(_distrPSA_varName, _distrPSA_distrName, _distrPSA_param1, _distrPSA_param2);
	for(int n=0; n<num; n++){
		modelParamType txArm;
		baseCohortType testCohort(s_F0,49.9,'M', 'W', 2, 'N', 'A');	
		mySim.PSA_sampleModelParamValue(txArm);	// sample cost/probability/quality value
		mySim.PSA_sampleCohort(testCohort);	// sample the state/age/gender value



		outf<<testCohort.baseAge<<"\t"<<testCohort.baseState<<"\t"<<testCohort.baseGender<<"\t"
			<<txArm._qolData.q_F0<<"\t"
			<<txArm._qolData.q_F1<<"\t"
			<<txArm._qolData.q_F2<<"\t"
			<<txArm._qolData.q_F3<<"\t"
			<<txArm._qolData.q_CoCirr<<"\t"
			<<txArm._qolData.q_DeCirr<<"\t"
			<<txArm._qolData.q_HCC<<"\t"
			<<txArm._qolData.q_LivTr<<"\t"
			<<txArm._qolData.q_PostLivTr<<"\t"
			<<txArm._qolData.q_SVR<<"\t"
			<<txArm._qolData.q_Dec_Anemia<<"\t"
			<<txArm._qolData.q_TX_oSOC<<"\t"
			<<txArm._qolData.q_TX_DAA<<"\t"

			<<txArm._costData.c_F0<<"\t"
			<<txArm._costData.c_F1<<"\t"
			<<txArm._costData.c_F2<<"\t"
			<<txArm._costData.c_F3<<"\t"
			<<txArm._costData.c_CoCirr<<"\t"
			<<txArm._costData.c_DeCirr<<"\t"
			<<txArm._costData.c_DeCirr1yrPlus<<"\t"
			<<txArm._costData.c_HCC<<"\t"
			<<txArm._costData.c_LivTr<<"\t"
			<<txArm._costData.c_PostLivTr<<"\t"

			<<txArm._transData.pr_F0_F1<<"\t"
			<<txArm._transData.pr_F1_F2<<"\t"
			<<txArm._transData.pr_F2_F3<<"\t"
			<<txArm._transData.pr_F3_CoCirr<<"\t"
			<<txArm._transData.pr_CoCirr_DeCirr<<"\t"
			<<txArm._transData.pr_CoCirr_HCC<<"\t"
			<<txArm._transData.pr_DeCirr_HCC<<"\t"
			<<txArm._transData.pr_DeCirr_LivTr<<"\t"
			<<txArm._transData.pr_DeCirr_DeathLiv<<"\t"
			<<txArm._transData.pr_DeCirr1yrPlus_DeathLiv<<"\t"
			<<txArm._transData.pr_HCC_LivTr<<"\t"
			<<txArm._transData.pr_HCC_DeathLiv<<"\t"
			<<txArm._transData.pr_LivTr_DeathLiv<<"\t"
			<<txArm._transData.pr_LivTr1yrPlus_DeathLiv<<"\t"
			<<txArm._transData.pr_SVR_CoCirr_DeCirr<<"\t"
			<<txArm._transData.pr_SVR_CoCirr_HCC<<"\t"
			<<txArm._transData.pr_SVR_Delta_oSOC<<"\t"
			<<txArm._transData.pr_SVR_Delta_DAA
			<<endl;


	}
	outf.close();
	return 0;
}

int Project_PSA_MetaModeling::RunModelWithDiffParamValue()
{
	_outf.open("output_project_metamodel/project_psa_diff_value.txt", ios::app);
	if (_outf.fail()) {
		ExitWithMsg("Error @ RunBaseCase(): Fail to open file output_PSA_metamodel/project_psa_diff_value.txt");
	}

	//string change_var = "c_F0"; double change_val = 500; // default value = 768 in data_type.cpp
	 string change_var = "c_F0"; double change_val = 1000; // default value = 768 in data_type.cpp
	//string change_var = "c_0"; double change_val = 1000; // WRONG parameter name
	
	ChangeValue(change_var, change_val);	

	_outf << "param\tvalue\t";
	PrintHeader(_outf);
	_outf << change_var<<"\t" << change_val << "\t";

	
	cout << "[" << currentDateTime() << "] Running base case ..." << endl;
	NUM_PATIENTS = 10000;
	Compare(0, 0);

	_outf << endl;
	_outf.close();

	return 0;

}



int Project_PSA_MetaModeling::ChangeValue(string varName, double argVal) {
	if ("q_F0" == varName) {
		_modelParam._qolData.q_F0 = argVal;
	}
	else if ("q_F1" == varName) {
		_modelParam._qolData.q_F1 = argVal;
	}
	else if ("q_F2" == varName) {
		_modelParam._qolData.q_F2 = argVal;
	}
	else if ("q_F3" == varName) {
		_modelParam._qolData.q_F3 = argVal;
	}
	else if ("q_CoCirr" == varName) {
		_modelParam._qolData.q_CoCirr = argVal;
	}
	else if ("q_DeCirr" == varName) {
		_modelParam._qolData.q_DeCirr = argVal;
	}
	else if ("q_HCC" == varName) {
		_modelParam._qolData.q_HCC = argVal;
	}
	else if ("q_LivTr" == varName) {
		_modelParam._qolData.q_LivTr = argVal;
	}
	else if ("q_PostLivT" == varName) {
		_modelParam._qolData.q_PostLivTr = argVal;
	}
	else if ("q_SVR" == varName) {
		_modelParam._qolData.q_SVR = argVal;
	}
	else if ("q_Anemia" == varName) {
		_modelParam._qolData.q_Dec_Anemia = argVal;
	}
	else if ("q_TX_oSOC" == varName) {
		_modelParam._qolData.q_TX_oSOC = argVal;
	}
	else if ("q_Tx_DAA" == varName) {
		_modelParam._qolData.q_TX_DAA = argVal;
	}
	else if ("c_F0" == varName) {
		_modelParam._costData.c_F0 = argVal;
	}
	else if ("c_F1" == varName) {
		_modelParam._costData.c_F1 = argVal;
	}
	else if ("c_F2" == varName) {
		_modelParam._costData.c_F2 = argVal;
	}
	else if ("c_F3" == varName) {
		_modelParam._costData.c_F3 = argVal;
	}
	else if ("c_CoCirr" == varName) {
		_modelParam._costData.c_CoCirr = argVal;
	}
	else if ("c_DeCirr" == varName) {
		_modelParam._costData.c_DeCirr = argVal;
	}
	else if ("c_DeCirr1yrPlus" == varName) {
		_modelParam._costData.c_DeCirr1yrPlus = argVal;
	}
	else if ("c_HCC" == varName) {
		_modelParam._costData.c_HCC = argVal;
	}
	else if ("c_LivTr" == varName) {
		_modelParam._costData.c_LivTr = argVal;
	}
	else if ("c_PostLivTr" == varName) {
		_modelParam._costData.c_PostLivTr = argVal;
	}
	else if ("pF0_F1_SA" == varName) {
		_modelParam._transData.pr_F0_F1 = argVal;
	}
	else if ("pF1_F2_SA" == varName) {
		_modelParam._transData.pr_F1_F2 = argVal;
	}
	else if ("pF2_F3_SA" == varName) {
		_modelParam._transData.pr_F2_F3 = argVal;
	}
	else if ("pF3_F4_SA" == varName) {
		_modelParam._transData.pr_F3_CoCirr = argVal;
	}
	else if ("pF4_DC_SA" == varName) {
		_modelParam._transData.pr_CoCirr_DeCirr = argVal;
	}
	else if ("pF4_HCC_SA" == varName) {
		_modelParam._transData.pr_CoCirr_HCC = argVal;
	}
	else if ("pDC_HCC_SA" == varName) {
		_modelParam._transData.pr_DeCirr_HCC = argVal;
	}
	else if ("pDC_Liv_Transpl_SA" == varName) {
		_modelParam._transData.pr_DeCirr_LivTr = argVal;
	}
	else if ("pMort_dc_cyc_1_SA" == varName) {
		_modelParam._transData.pr_DeCirr_DeathLiv = argVal;
	}
	else if ("pMort_dc_cyc_2_SA" == varName) {
		_modelParam._transData.pr_DeCirr1yrPlus_DeathLiv = argVal;
	}
	else if ("pHCC_Liv_Transpl_SA" == varName) {
		_modelParam._transData.pr_HCC_LivTr = argVal;
	}
	else if ("pMort_hcc_cyc_SA" == varName) {
		_modelParam._transData.pr_HCC_DeathLiv = argVal;
	}
	else if ("pMort_Liv_Transpl_SA" == varName) {
		_modelParam._transData.pr_LivTr_DeathLiv = argVal;
	}
	else if ("pMort_Post_Liv_Transpl_SA" == varName) {
		_modelParam._transData.pr_LivTr1yrPlus_DeathLiv = argVal;
	}
	else if ("pr_SVR_CoCirr_DeCirr" == varName) {
		_modelParam._transData.pr_SVR_CoCirr_DeCirr = argVal;
	}
	else if ("pr_SVR_CoCirr_HCC" == varName) {
		_modelParam._transData.pr_SVR_CoCirr_HCC = argVal;
	}
	else if ("pSVR_Delta_oSOC" == varName) {
		_modelParam._transData.pr_SVR_Delta_oSOC = argVal;
	}
	else if ("pSVR_Delta_DAA" == varName) {
		_modelParam._transData.pr_SVR_Delta_DAA = argVal;


	}
	else {
		ExitWithMsg("[Error] Project_PSA_MetaModeling::ChangeValue(string varName, double argVal): Unknown parameter " + varName);
	}
	return 0;
}
