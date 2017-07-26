#include "project_hcvcalculator.h"
using namespace std;



int Project_HCVCalculator::RunBaseCase()
{
	_strOutputFile="output_hcvcalc/project_hcvcalculator_output_age"+basicToStr(_initialAge)+"_0.txt";
	_outFileSummary.open(_strOutputFile);
	_outFileSummary	<< fixed << showpoint;

	// print the table header:
	//for(vector<string>::const_iterator it=_vecVarOfInterest.begin(); it!=_vecVarOfInterest.end(); it++){
	//	_outFileSummary<<(*it)<<"\t";
	//}
	//_outFileSummary<<"Fib\tAge\tGender\tGenotype\tTxExp\t";
	//_outFileSummary<<"QALY1\tQALY2\tIncrQALY\tCost1\tCost2\tIncrCost\tICER\t";
	//_outFileSummary<<"DC1\tDC2\tHCC1\tHCC2\tLT1\tLT2\tLRD1\tLRD2\tTxCost1\tTxCost2\t";
	
	
	CompareAllCases();

	_outFileSummary<<endl;
	_outFileSummary.close();

	return 0;
}


int Project_HCVCalculator::RunPerturbations(string argPertValFile, int idxPertub)
{
	ReadPerturbedValues(argPertValFile);
	_baselineModelParam=_modelParam;// backup the baseline model values;

	if(idxPertub == 0){
		_outFileSummary.open("output_hcvcalc/project_hcvcalculator_output_age"+basicToStr(_initialAge)+"_"+basicToStr(idxPertub)+".txt",ios::app);
	}else{
		_outFileSummary.open("output_hcvcalc/project_hcvcalculator_output_age"+basicToStr(_initialAge)+"_"+basicToStr(idxPertub)+".txt");
	}
	_outFileSummary	<< fixed << showpoint;


	//for(int k=0; k<_vecPertVarName.size(); k++){
	int k=idxPertub;

		_modelParam=_baselineModelParam; // reset the model parameters;
		ChangeValue(_vecPertVarName[k],_vecPertVarVal[k],_modelParam); // change the value;

		cout<<"----- Change variable: "<<_vecPertVarName[k]<<"="<<_vecPertVarVal[k]<<" ---------"<<endl;
		CompareAllCases();
	//}

	_outFileSummary.close();
	return 0;
}

int Project_HCVCalculator::CompareAllCases(){

	// ========================================================================
	//							Before Aug 2015
	// ========================================================================

	////TN TOL
	////Compare("TN_TOL_G1_BOC","TN_TOL_G1_SOF_PEG_RBV12",1);
	////Compare("TN_TOL_G1_TEL","TN_TOL_G1_SOF_PEG_RBV12",1);

	////Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF8",1);
	////Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF12",1);
	////Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF",1);
	////Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF8",1);
	////Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF12",1);
	////Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF",1);

	//_strTxExperience="TN_TOL";
	//Compare("TN_TOL_G1_BOC_TEL","TN_G1_LDP_SOF",1);

	//Compare("TN_TOL_G2_PEG_RBV24","TN_TOL_G2_SOF_RBV12",2);
	//Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_RBV24",3);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_PEG_RBV12",4);

	//// TN NOT-TOL

	////Compare("NoTx","TN_NOT_G1_SOF_SMV12",1);

	////Compare("NoTx","TN_G1_LDP_SOF8",1);
	////Compare("NoTx","TN_G1_LDP_SOF12",1);
	//_strTxExperience="TN_INTOL";
	//Compare("NoTx","TN_G1_LDP_SOF",1);
	//Compare("NoTx","TN_NOT_G2_SOF_RBV12",2);
	//Compare("NoTx","TN_NOT_G3_SOF_RBV24",3);
	//Compare("NoTx","TN_NOT_G4_SOF_RBV24",4);

	//// TE
	////Compare("TE_G1_BOC","TE_G1_SOF_SMV12",1);
	////Compare("TE_G1_TEL","TE_G1_SOF_SMV12",1);

	////Compare("TE_G1_BOC","TE_G1_LDP_SOF12_24",1);
	////Compare("TE_G1_TEL","TE_G1_LDP_SOF12_24",1);
	//_strTxExperience="TE";
	//Compare("TE_G1_BOC_TEL","TE_G1_LDP_SOF12_24",1);
	//Compare("TE_G2_PEG_RBV24","TE_G2_SOF_RBV12",2);
	//Compare("TE_G3_PEG_RBV24","TE_G3_SOF_RBV24",3);
	//Compare("TE_G4_PEG_RBV48","TE_G4_SOF_PEG_RBV12",4);



	// ========================================================================
	// Updated on Aug 27, 2015
	// genotype 11 = 1a
	// genotype 12 = 1b
	// - 1a and 1b share the same treatment effectiveness of TN_G1_LDP_SOF and TE_G1_LDV_SOF12_24
	// ========================================================================
	_strTxExperience="TN_TOL";
	Compare("TN_TOL_G1_BOC_TEL","TN_G1_LDP_SOF",11);
	Compare("TN_TOL_G1_BOC_TEL","TN_TOL_G1a_DCV_SOF12_24",11);
	Compare("TN_TOL_G1_BOC_TEL","TN_TOL_G1a_PrOD_RBV12_24",11);

	Compare("TN_TOL_G1_BOC_TEL","TN_G1_LDP_SOF",12);
	Compare("TN_TOL_G1_BOC_TEL","TN_TOL_G1b_DCV_SOF12_24",12);
	Compare("TN_TOL_G1_BOC_TEL","TN_TOL_G1b_PrOD12",12);
	Compare("TN_TOL_G1_BOC_TEL","TN_TOL_G1b_SOF_SMV12",12);

	Compare("TN_TOL_G2_PEG_RBV24","TN_TOL_G2_SOF_RBV12_16",2);

	Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_RBV24",3);
	Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_DCV_SOF12_24",3);
	Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_PEG_RBV12",3);

	Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_LDV_SOF12",4);
	Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_PrOD_RBV12",4);
	Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_RBV24",4);
	Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_PEG_RBV12",4);


	//_strTxExperience="TN_INTOL";
	Compare("NoTx","TN_G1_LDP_SOF",11);
	Compare("NoTx","TN_TOL_G1a_DCV_SOF12_24",11);
	Compare("NoTx","TN_TOL_G1a_PrOD_RBV12_24",11);

	Compare("NoTx","TN_G1_LDP_SOF",12);
	Compare("NoTx","TN_TOL_G1b_DCV_SOF12_24",12);
	Compare("NoTx","TN_TOL_G1b_PrOD12",12);
	Compare("NoTx","TN_TOL_G1b_SOF_SMV12",12);

	Compare("NoTx","TN_TOL_G2_SOF_RBV12_16",2);

	Compare("NoTx","TN_TOL_G3_SOF_RBV24",3);
	Compare("NoTx","TN_TOL_G3_DCV_SOF12_24",3);
	Compare("NoTx","TN_TOL_G3_SOF_PEG_RBV12",3);


	Compare("NoTx","TN_TOL_G4_LDV_SOF12",4);
	Compare("NoTx","TN_TOL_G4_PrOD_RBV12",4);
	Compare("NoTx","TN_TOL_G4_SOF_RBV24",4);
	Compare("NoTx","TN_TOL_G4_SOF_PEG_RBV12",4);

	_strTxExperience="TE";
	Compare("TE_G1_BOC_TEL","TE_G1_LDP_SOF12_24",11);
	Compare("TE_G1_BOC_TEL","TE_G1a_DCV_SOF12_24",11);
	Compare("TE_G1_BOC_TEL","TE_G1a_PrOD_RBV12_24",11);
	Compare("TE_G1_BOC_TEL","TE_G1a_SMV_SOF12",11);

	Compare("TE_G1_BOC_TEL","TE_G1_LDP_SOF12_24",12);
	Compare("TE_G1_BOC_TEL","TE_G1b_DCV_SOF12_24",12);
	Compare("TE_G1_BOC_TEL","TE_G1b_PrOD12",12);
	Compare("TE_G1_BOC_TEL","TE_G1b_SMV_SOF12",12);
	
	Compare("TE_G2_PEG_RBV24","TE_G2_LDV_RBV12_16",2);
	Compare("TE_G2_PEG_RBV24","TE_G2_SOF_PEG_RBV12",2);
	Compare("TE_G2_PEG_RBV24","TE_G2_SMV_SOF12",2);
	
	Compare("TE_G3_PEG_RBV24","TE_G3_DCV_SOF12_24",3);
	Compare("TE_G3_PEG_RBV24","TE_G3_SOF_PEG_RBV12",3);
	
	
	Compare("TE_G4_PEG_RBV48","TE_G4_LDV_SOF12",4);
	Compare("TE_G4_PEG_RBV48","TE_G4_PrOD_RBV12",4);
	Compare("TE_G4_PEG_RBV48","TE_G4_SOF_RBV24",4);
	Compare("TE_G4_PEG_RBV48","TE_G4_SOF_PEG_RBV12",4);
	return 0;
}

int Project_HCVCalculator::Compare(string strArm1, string strArm2, int argGenotype)
{





	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios=10;
	stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F','M','F','M','F'};
	double listAge[]={49.9,49.9,56.2,56.2,57.6,57.6,58.3,58.3,58.7,58.7};
	if(_initialAge - 0 > EPSILON){// if specific initial age is assigned
		for(int ageK=0; ageK<10; ageK++){
			listAge[ageK]=_initialAge;
		}
	}
	//int nScenarios=30;
	//stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr,
	//					   s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr,
	//					   s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	//char listGender[]={'M','F','M','F','M','F','M','F','M','F',
	//				   'M','F','M','F','M','F','M','F','M','F',
	//				   'M','F','M','F','M','F','M','F','M','F'};
	//double listAge[]={40,40,40,40,40,40,40,40,40,40,
	//				  55,55,55,55,55,55,55,55,55,55,
	//				  70,70,70,70,70,70,70,70,70,70};


	for(int k=0; k<nScenarios; k++){
		cout<<"Comparing\t"<<strArm1<<"-"<<strArm2<<":\tFibrosis="<<listState[k]<<"\tAge="<<listAge[k]
		<<"\tGender="<<listGender[k]<<"\tGenotype="<<argGenotype<<"\tTxExp="<<_strTxExperience<<endl;
		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]

		int genotype_noSubtype = (argGenotype == 11 || argGenotype == 12)?1:argGenotype;
		baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', genotype_noSubtype, 'N', 'A');	

		//vector<double> r1=EvaluateOneArm(strArm1, testCohort);
		//vector<double> r2=EvaluateOneArm(strArm2, testCohort);
		vector<double> r1=EvaluateOneArm_Averaged(strArm1, testCohort);
		vector<double> r2=EvaluateOneArm_Averaged(strArm2, testCohort);

		assert(r1.size()==r2.size());

		for(vector<string>::const_iterator it=_vecVarOfInterest.begin(); it!=_vecVarOfInterest.end(); it++){
			_outFileSummary<<GetValue(*it)<<"\t";
		}
		_outFileSummary<<strArm1<<"\t"<<strArm2<<"\t";
		_outFileSummary<<listState[k]<<"\t"<<listAge[k]<<"\t"<<listGender[k]<<"\t"<<argGenotype<<"\t"<<_strTxExperience<<"\t";
		_outFileSummary<<r1[0]<<"\t"<<r2[0]<<"\t"<<r2[0]-r1[0]<<"\t"
		<<r1[1]<<"\t"<<r2[1]<<"\t"<<r2[1]-r1[1]<<"\t"
			<<(r2[1]-r1[1])/(r2[0]-r1[0])<<"\t";

		for(int k=2;k<r1.size(); k++){
			_outFileSummary<<r1[k]<<"\t"<<r2[k]<<"\t";
		}
		_outFileSummary<<endl;
	}


	return 0;
}

vector<double> Project_HCVCalculator::EvaluateOneArm( string strArm, const baseCohortType & testCohort )
{
	modelParamType arm=_modelParam;
	arm._armName=strArm;
	arm._cohortData=testCohort;
	arm.ReduceDrugCost(_drugCostReduction);

	HepCSim * mySim;
	mySim=new HepCSim;
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

	FreeMem(mySim);
	return r;
	
}

vector<double> Project_HCVCalculator::EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort )
{
	vector<double> r;
	if("TN_TOL_G1_BOC_TEL" == txArm){
		vector<double> r_boc=EvaluateOneArm("TN_TOL_G1_BOC", testCohort);
		vector<double> r_tel=EvaluateOneArm("TN_TOL_G1_TEL", testCohort);
		r=LinearCombTwoVec(r_boc, r_tel, _modelParam._transData._pBOC);

	}else if ("TE_G1_BOC_TEL" == txArm){
		vector<double> r_boc=EvaluateOneArm("TE_G1_BOC", testCohort);
		vector<double> r_tel=EvaluateOneArm("TE_G1_TEL", testCohort);
		r=LinearCombTwoVec(r_boc, r_tel, _modelParam._transData._pBOC);

	}else if("TN_G1_LDP_SOF" == txArm){		
		vector<double> r_8week=EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
		vector<double> r_12week=EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
		r=LinearCombTwoVec(r_8week, r_12week, _modelParam._transData._p8Week_LDV);
		//assert(r_8week.size()==r_12week.size());
		//for(int l=0; l<r_8week.size(); l++){
		//	r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1-_modelParam._transData._p8Week_LDV)*r_12week[l]);
		//}
	}else{
		r=EvaluateOneArm(txArm, testCohort);
	}
	return r;
}



int Project_HCVCalculator::ChangeValue(string varName, double argVal, modelParamType & argModelParam){
	if ("q_F0"== varName){argModelParam._qolData.q_F0 = argVal;
	}else if ("q_F1"== varName){argModelParam._qolData.q_F1 = argVal;
	}else if ("q_F2"== varName){argModelParam._qolData.q_F2 = argVal;
	}else if ("q_F3"== varName){argModelParam._qolData.q_F3 = argVal;
	}else if ("q_CoCirr"== varName){argModelParam._qolData.q_CoCirr = argVal;
	}else if ("q_DeCirr"== varName){argModelParam._qolData.q_DeCirr = argVal;
	}else if ("q_HCC"== varName){argModelParam._qolData.q_HCC = argVal;
	}else if ("q_LivTr"== varName){argModelParam._qolData.q_LivTr = argVal;
	}else if ("q_PostLivT"== varName){argModelParam._qolData.q_PostLivTr = argVal;
	}else if ("q_SVR"== varName){argModelParam._qolData.q_SVR = argVal;
	}else if ("q_Anemia"== varName){argModelParam._qolData.q_Dec_Anemia = argVal;
	}else if ("q_TX_oSOC"== varName){argModelParam._qolData.q_TX_oSOC = argVal;
	}else if ("q_Tx_DAA"== varName){argModelParam._qolData.q_TX_DAA = argVal;
	}else if ("c_F0"== varName){argModelParam._costData.c_F0 = argVal;
	}else if ("c_F1"== varName){argModelParam._costData.c_F1 = argVal;
	}else if ("c_F2"== varName){argModelParam._costData.c_F2 = argVal;
	}else if ("c_F3"== varName){argModelParam._costData.c_F3 = argVal;
	}else if ("c_CoCirr"== varName){argModelParam._costData.c_CoCirr = argVal;
	}else if ("c_DeCirr"== varName){argModelParam._costData.c_DeCirr = argVal;
	}else if ("c_DeCirr1yrPlus"== varName){argModelParam._costData.c_DeCirr1yrPlus = argVal;
	}else if ("c_HCC"== varName){argModelParam._costData.c_HCC = argVal;
	}else if ("c_LivTr"== varName){argModelParam._costData.c_LivTr = argVal;
	}else if ("c_PostLivTr"== varName){argModelParam._costData.c_PostLivTr = argVal;
	}else if ("pF0_F1_SA"== varName){argModelParam._transData.pr_F0_F1 = argVal;
	}else if ("pF1_F2_SA"== varName){argModelParam._transData.pr_F1_F2 = argVal;
	}else if ("pF2_F3_SA"== varName){argModelParam._transData.pr_F2_F3 = argVal;
	}else if ("pF3_F4_SA"== varName){argModelParam._transData.pr_F3_CoCirr = argVal;
	}else if ("pF4_DC_SA"== varName){argModelParam._transData.pr_CoCirr_DeCirr = argVal;
	}else if ("pF4_HCC_SA"== varName){argModelParam._transData.pr_CoCirr_HCC = argVal;
	}else if ("pDC_HCC_SA"== varName){argModelParam._transData.pr_DeCirr_HCC = argVal;
	}else if ("pDC_Liv_Transpl_SA"== varName){argModelParam._transData.pr_DeCirr_LivTr = argVal;
	}else if ("pMort_dc_cyc_1_SA"== varName){argModelParam._transData.pr_DeCirr_DeathLiv = argVal;
	}else if ("pMort_dc_cyc_2_SA"== varName){argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = argVal;
	}else if ("pHCC_Liv_Transpl_SA"== varName){argModelParam._transData.pr_HCC_LivTr = argVal;
	}else if ("pMort_hcc_cyc_SA"== varName){argModelParam._transData.pr_HCC_DeathLiv = argVal;
	}else if ("pMort_Liv_Transpl_SA"== varName){argModelParam._transData.pr_LivTr_DeathLiv = argVal;
	}else if ("pMort_Post_Liv_Transpl_SA"== varName){argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = argVal;
	}else if ("pr_SVR_CoCirr_DeCirr"== varName){argModelParam._transData.pr_SVR_CoCirr_DeCirr = argVal;
	}else if ("pr_SVR_CoCirr_HCC"== varName){argModelParam._transData.pr_SVR_CoCirr_HCC = argVal;
	}else if ("pSVR_Delta_oSOC"== varName){argModelParam._transData.pr_SVR_Delta_oSOC = argVal;
	}else if ("pSVR_Delta_DAA"== varName){argModelParam._transData.pr_SVR_Delta_DAA = argVal;
	}
	// ------ added 6/26/2015, QC, for HCV calculator	----------
	else if ("cost_TEL"== varName){argModelParam._costData.c_TEL = argVal;
	}else if ("cost_BOC"== varName){argModelParam._costData.c_BOC = argVal;
	}else if ("cost_RBV"== varName){argModelParam._costData.c_RBV = argVal;
	}else if ("cost_PEG"== varName){argModelParam._costData.c_PEG = argVal;
	}else if ("cost_SMV"== varName){argModelParam._costData.c_SMV = argVal;
	}else if ("cost_SOF"== varName){argModelParam._costData.c_SOF = argVal;
	}else if ("cost_LDV"== varName){argModelParam._costData.c_LDV = argVal;
	}else if ("cost_DCV"== varName){argModelParam._costData.c_DCV = argVal;
	}else if ("cost_PrOD"== varName){argModelParam._costData.c_PrOD = argVal;

	// -----------------------------------------------------------
	
	}else{
		ExitWithMsg("[Error] Project_Doublecheck::ChangeValue(string varName, double argVal): Unknown parameter "+varName);
	}
	return 0;
}

Project_HCVCalculator::Project_HCVCalculator()
{
	_initialAge=0;
	_drugCostReduction=0;
	// read all arms for comparisons
	ifstream inf;
	inf.open("project_hcvcalculator_input_varOfInterest.txt");
	ReadList(_vecVarOfInterest,inf);
	inf.close();

	_outFileSummary.open(_strOutputFile);
	_outFileSummary.close();

}

int Project_HCVCalculator::ReadPerturbedValues( string argPertValFile )
{
	ifstream inf;
	inf.open(argPertValFile);
	string word;
	while(inf>>word){
		_vecPertVarName.push_back(word);
		double v;
		inf>>v;
		_vecPertVarVal.push_back(v);
	}

	inf.close();
	return 0;
}



double Project_HCVCalculator::GetValue( string varName )
{
	if ("cost_TEL"== varName){return _modelParam._costData.c_TEL;
	}else if ("cost_BOC"== varName){return _modelParam._costData.c_BOC;
	}else if ("cost_RBV"== varName){return _modelParam._costData.c_RBV;
	}else if ("cost_PEG"== varName){return _modelParam._costData.c_PEG;
	}else if ("cost_SMV"== varName){return _modelParam._costData.c_SMV;
	}else if ("cost_SOF"== varName){return _modelParam._costData.c_SOF;
	}else if ("cost_LDV"== varName){return _modelParam._costData.c_LDV;
	}else if ("cost_DCV"== varName){return _modelParam._costData.c_DCV;
	}else if ("cost_PrOD"== varName){return _modelParam._costData.c_PrOD;
	// -----------------------------------------------------------
	}else{
		ExitWithMsg("[Error] Project_Doublecheck::ChangeValue(string varName, double argVal): Unknown parameter "+varName);
		return 1;
	}
}




