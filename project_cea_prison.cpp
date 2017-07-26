#include "project_cea_prison.h"
using namespace std;




int Project_CEA_Prison::BaseCase(){
	ifstream inFileLHS;
	ofstream outFileSummary;


	outFileSummary.open("project_cea_prison_output_basecase.txt");
	outFileSummary<<endl<<"===== "<<currentDateTime()<<" ====="<<endl
		<< fixed << showpoint;
	// print the table header:
	outFileSummary << "Scenario" << "\t" << "Treat"<<"\t"
		<< "QALY" << "\t" << "COST" << "\t" << "COUNT_DECIRR" << "\t"
		<< "COUNT_HCC" << "\t" << "COUNT_LIVTR" << "\t" << "COUNT_DEATHLIV" 
		<< "\t" << "COST_TX" << endl;


	// ==== prepare the model parameters ====
	vector<modelParamType> vecArms;
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios=10;
	stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F','M','F','M','F'};
	double listAge[]={49.9,49.9,56.2,56.2,57.6,57.6,58.3,58.3,58.7,58.7};

	// === arms ===
	string listTxArm[]={"TN_G1_LDP_SOF8",
		"TN_G1_LDP_SOF12",
		"TN_G1_TOL_BOC",
		"TN_G1_TOL_TEL",
		"TE_G1_LDP_SOF12_24",
		"TE_G1_BOC",
		"TE_G1_TEL"};
	int listGenotype[]={1,1,1,1,1,1,1};


	for(int a=0; a<sizeof(listGenotype)/sizeof(int); a++){
		//if(a!=2) continue; // only run arm 6

		for(int k=0; k<nScenarios; k++){
			//		if(k!=2) continue; // only run scenario1

			if(listTxArm[a]=="TN_G1_LDP_SOF8" && listState[k]==s_CoCirr) continue;	

			// define cohort
			// create patient profile:
			// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
			baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', listGenotype[a], 'N', 'A');	

			modelParamType armX(listTxArm[a],testCohort); // arm 10 in the original code
			armX.ReduceDrugCost(_drugCostReduction);

			//HepCSim * mySim;
			//mySim=new HepCSim;
			//mySim->SetRandomSeed(SIM_SEED);
			//mySim->Run(armX);
			//// output one line summary of this arm
			//outFileSummary<<armX._armName<<"\t"<<TREATMENT<<"\t";
			//mySim->OutputBaseCase(outFileSummary);
			//FreeMem(mySim);

			HepCSim mySim;
			mySim.SetRandomSeed(SIM_SEED);
			mySim.Run(armX);
			// output one line summary of this arm
			outFileSummary<<armX._armName<<"\t"<<TREATMENT<<"\t";
			mySim.OutputBaseCase(outFileSummary);

		}


	}


	/// ==== no treatment for TN-INTOL
	TREATMENT = false;
	for(int k=0; k<nScenarios; k++){

		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', 1, 'N', 'A');	

		modelParamType armX("NoTx",testCohort); // arm 10 in the original code
		armX.ReduceDrugCost(_drugCostReduction);

		HepCSim mySim;
		mySim.SetRandomSeed(SIM_SEED);
		mySim.Run(armX);
		// output one line summary of this arm
		outFileSummary<<armX._armName<<"\t"<<TREATMENT<<"\t";
		mySim.OutputBaseCase(outFileSummary);

	}

	outFileSummary.close();
	return 0;
}

int Project_CEA_Prison::PSA(){
	// read parameter values from external file:
	ifstream inFileLHS;
	inFileLHS.open(FILE_LHS.c_str());
	if (!inFileLHS) //exit the program if any of the output files not found
	{
		cout << "Cannot open onput file LHSinput.in"<<endl<< "The program terminates." << endl;
		return 1;
	}


	inFileLHS.close();

	return 0;
}


int Project_CEA_Prison::RunDifferentMCReplc(){
	ofstream outFileSummary;


	outFileSummary.open("project_cea_prison_output_mcreplc.txt");
	outFileSummary<<endl<<"===== "<<currentDateTime()<<" ====="<<endl
		<< fixed << setprecision(4);
	// print the table header:


	// ==== prepare the model parameters ====
	vector<modelParamType> vecArms;
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	int nScenarios=10;
	stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F','M','F','M','F'};
	double listAge[]={49.9,49.9,56.2,56.2,57.6,57.6,58.3,58.3,58.7,58.7};
	string listTxArm[]={"TN_TOL_G1_SOF_PEG_RBV12",
		"TN_TOL_G2_SOF_RBV12",
		"TN_TOL_G2_PEG_RBV24"};
	int listGenotype[]={1,2,2,3,4,1,2,3,4,1,2,3,4};

	for(int a=0; a<3; a++){
		if(a!=2) continue; // only run arm 0

		for(int k=0; k<nScenarios; k++){
			if(k!=0) continue; // only run scenario 1

			// define cohort
			// create patient profile:
			// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
			baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', listGenotype[a], 'N', 'A');	

			modelParamType armX(listTxArm[a],testCohort); // arm 10 in the original code
			armX.ReduceDrugCost(_drugCostReduction);
			NUM_PATIENTS=1000000;

			HepCSim mySim;
			mySim.SetPrintResultForEveryPatient(true);
			mySim.SetRandomSeed(SIM_SEED);
			mySim.Run(armX);
			// output one line summary of this arm

		}
	}
	outFileSummary.close();
	return 0;
}



int Project_CEA_Prison::CEA_BaseResults()
{
	//Compare("FBOP","TN_G1_LDP_SOF",1);
	Compare("FBOP","TN_G1_LDP_SOF12",1);

	


	return 0;
}

int Project_CEA_Prison::Compare(string strArm1, string strArm2, int argGenotype)
{


	_outFileSummary.open("project_cea_prison_output_compare.txt",ios::app);
	_outFileSummary<<endl<<"===== "<<currentDateTime()<<" ===== "<<strArm1<<" vs "<<strArm2<<" ==="<<endl
		<< fixed << showpoint;
	// print the table header:



	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios=6;
	stateType listState[]={s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F'};
	double listAge[]={50.7,50.7,50.7,50.7,50.6,50.6};


	for(int k=0; k<nScenarios; k++){
		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', argGenotype, 'N', 'A');	

		//vector<double> r1=EvaluateOneArm(strArm1, testCohort);
		//vector<double> r2=EvaluateOneArm(strArm2, testCohort);
		vector<double> r1=EvaluateOneArm_CombineLDV8Or12WkOnline(strArm1, testCohort);
		vector<double> r2=EvaluateOneArm_CombineLDV8Or12WkOnline(strArm2, testCohort);

		assert(r1.size()==r2.size());
		_outFileSummary<<r1[0]<<"\t"<<r2[0]<<"\t"<<r2[0]-r1[0]<<"\t"
		<<r1[1]<<"\t"<<r2[1]<<"\t"<<r2[1]-r1[1]<<"\t"
			<<(r2[1]-r1[1])/(r2[0]-r1[0])<<"\t";

		for(int k=2;k<r1.size(); k++){
			_outFileSummary<<r1[k]<<"\t"<<r2[k]<<"\t";
		}
		_outFileSummary<<endl;
	}

	_outFileSummary<<endl;
	_outFileSummary.close();
	return 0;
}

vector<double> Project_CEA_Prison::EvaluateOneArm( string strArm, const baseCohortType & testCohort )
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

vector<double> Project_CEA_Prison::EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort )
{

	if("TN_G1_LDP_SOF" == txArm){
		vector<double> r;
		vector<double> r_8week=EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
		vector<double> r_12week=EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
		assert(r_8week.size()==r_12week.size());
		for(int l=0; l<r_8week.size(); l++){
			r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1-_modelParam._transData._p8Week_LDV)*r_12week[l]);
		}
		return r;
	}else{
		return EvaluateOneArm(txArm, testCohort);
	}

}

int Project_CEA_Prison::CEA_OneWay(int argCmpIdx)
{

	

	ofstream outf_low, outf_high;
	outf_low.open("output_project_cea_prison/SA_LOW_"+_listCmp[argCmpIdx]+".txt");
	outf_high.open("output_project_cea_prison/SA_HIGH_"+_listCmp[argCmpIdx]+".txt");
	outf_low<< fixed << showpoint;
	outf_high<< fixed << showpoint;

	vector<string> sa_varName;
	vector<SParamTriplet> sa_varVal;
	ReadOneWaySAParamValueRange("project_cea_prison_input_OneWay.txt",sa_varName,sa_varVal);

	int sa_counter=1;
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios=10;
	stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F','M','F','M','F'};
	double listAge[]={49.9,49.9,56.2,56.2,57.6,57.6,58.3,58.3,58.7,58.7};


	for(int k=0; k<nScenarios; k++){
		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', _listGenotype[argCmpIdx], 'N', 'A');	

		// ===== baseline ======
		
		_modelParam=_baselineModelParam;// recover the baseline values

		vector<double> r1_b=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
		vector<double> r2_b=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

		assert(r1_b.size()==r2_b.size());
		outf_low<<"F"<<listState[k]<<" "<<listGender[k]<<" "<<listAge[k]<<"\t"<<sa_counter<<"\t"
			<<r1_b[0]<<"\t"<<r2_b[0]<<"\t"<<r2_b[0]-r1_b[0]<<"\t"
			<<r1_b[1]<<"\t"<<r2_b[1]<<"\t"<<r2_b[1]-r1_b[1]<<"\t"
			<<(r2_b[1]-r1_b[1])/(r2_b[0]-r1_b[0])<<"\t";
		outf_high<<"F"<<listState[k]<<" "<<listGender[k]<<" "<<listAge[k]<<"\t"<<sa_counter<<"\t"
			<<r1_b[0]<<"\t"<<r2_b[0]<<"\t"<<r2_b[0]-r1_b[0]<<"\t"
			<<r1_b[1]<<"\t"<<r2_b[1]<<"\t"<<r2_b[1]-r1_b[1]<<"\t"
			<<(r2_b[1]-r1_b[1])/(r2_b[0]-r1_b[0])<<"\t";
		for(int i=2;i<r1_b.size(); i++){
			outf_low<<r1_b[i]<<"\t"<<r2_b[i]<<"\t";
			outf_high<<r1_b[i]<<"\t"<<r2_b[i]<<"\t";
		}
		outf_low<<endl;
		outf_high<<endl;

		sa_counter++;

		for(int j=0; j<sa_varName.size(); j++){

		// ===== low value ======
			_modelParam=_baselineModelParam;// recover the baseline values
			ChangeValue(sa_varName[j],sa_varVal[j]._lb,_modelParam);
			vector<double> r1_low=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_low=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

			outf_low<<"F"<<listState[k]<<" "<<listGender[k]<<" "<<listAge[k]<<"\t"<<sa_counter<<"\t"
				<<r1_low[0]<<"\t"<<r2_low[0]<<"\t"<<r2_low[0]-r1_low[0]<<"\t"
				<<r1_low[1]<<"\t"<<r2_low[1]<<"\t"<<r2_low[1]-r1_low[1]<<"\t"
				<<(r2_low[1]-r1_low[1])/(r2_low[0]-r1_low[0])<<"\t";
			
			for(int i=2; i<r1_low.size(); i++){
				outf_low<<r1_low[i]<<"\t"<<r2_low[i]<<"\t";
			}
			outf_low<<endl;

		// ===== high value ======
			_modelParam=_baselineModelParam;// recover the baseline values
			ChangeValue(sa_varName[j],sa_varVal[j]._ub,_modelParam);
			vector<double> r1_high=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
			vector<double> r2_high=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

			outf_high<<"F"<<listState[k]<<" "<<listGender[k]<<" "<<listAge[k]<<"\t"<<sa_counter<<"\t"
				<<r1_high[0]<<"\t"<<r2_high[0]<<"\t"<<r2_high[0]-r1_high[0]<<"\t"
				<<r1_high[1]<<"\t"<<r2_high[1]<<"\t"<<r2_high[1]-r1_high[1]<<"\t"
				<<(r2_high[1]-r1_high[1])/(r2_high[0]-r1_high[0])<<"\t";
			for(int i=2;i<r1_high.size(); i++){
				outf_high<<r1_high[i]<<"\t"<<r2_high[i]<<"\t";
			}
			outf_high<<endl;

			sa_counter++;
		}
		
	}


	_modelParam=_baselineModelParam;//restore the nominal values



	outf_low.close(); outf_high.close();
	return 0;
}

int Project_CEA_Prison::ReadPSASampledValues( ifstream & inf, modelParamType & argModelParam )
{
	//update the model parameters:
	//Fibrosis score (F0-F4)	sex	age	
	// q_F0	q_F1	q_F2	q_F3	q_CoCirr	q_DeCirr	q_HCC	q_LivTr	q_PostLivT	q_SVR	q_Anemia	q_TX_oSOC	q_Tx_DAA	
	inf>>argModelParam._qolData.q_F0;
	inf>>argModelParam._qolData.q_F1;
	inf>>argModelParam._qolData.q_F2;
	inf>>argModelParam._qolData.q_F3;
	inf>>argModelParam._qolData.q_CoCirr;
	inf>>argModelParam._qolData.q_DeCirr;
	inf>>argModelParam._qolData.q_HCC;
	inf>>argModelParam._qolData.q_LivTr;
	inf>>argModelParam._qolData.q_PostLivTr;
	inf>>argModelParam._qolData.q_SVR;
	inf>>argModelParam._qolData.q_Dec_Anemia;
	inf>>argModelParam._qolData.q_TX_oSOC;
	inf>>argModelParam._qolData.q_TX_DAA;

	//c_F0	c_F1	c_F2	c_F3	c_CoCirr	c_DeCirr	c_DeCirr1yrPlus	c_HCC	c_LivTr	c_PostLivTr	
	inf>>argModelParam._costData.c_F0;
	inf>>argModelParam._costData.c_F1;
	inf>>argModelParam._costData.c_F2;
	inf>>argModelParam._costData.c_F3;
	inf>>argModelParam._costData.c_CoCirr;
	inf>>argModelParam._costData.c_DeCirr;
	inf>>argModelParam._costData.c_DeCirr1yrPlus;
	inf>>argModelParam._costData.c_HCC;
	inf>>argModelParam._costData.c_LivTr;
	inf>>argModelParam._costData.c_PostLivTr;

	//pF0_F1_SA	pF1_F2_SA	pF2_F3_SA	pF3_F4_SA	pF4_DC_SA	pF4_HCC_SA	pDC_HCC_SA	pDC_Liv_Transpl_SA	
	inf>>argModelParam._transData.pr_F0_F1;
	inf>>argModelParam._transData.pr_F1_F2;
	inf>>argModelParam._transData.pr_F2_F3;
	inf>>argModelParam._transData.pr_F3_CoCirr;
	inf>>argModelParam._transData.pr_CoCirr_DeCirr;
	inf>>argModelParam._transData.pr_CoCirr_HCC;
	inf>>argModelParam._transData.pr_DeCirr_HCC;
	inf>>argModelParam._transData.pr_DeCirr_LivTr;

	//pMort_dc_cyc_1_SA	pMort_dc_cyc_2_SA	pHCC_Liv_Transpl_SA	pMort_hcc_cyc_SA	pMort_Liv_Transpl_SA	pMort_Post_Liv_Transpl_SA	
	inf>>argModelParam._transData.pr_DeCirr_DeathLiv;
	inf>>argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv;
	inf>>argModelParam._transData.pr_HCC_LivTr;
	inf>>argModelParam._transData.pr_HCC_DeathLiv;
	inf>>argModelParam._transData.pr_LivTr_DeathLiv;
	inf>>argModelParam._transData.pr_LivTr1yrPlus_DeathLiv;

	//pr_SVR_CoCirr_DeCirr	pr_SVR_CoCirr_HCC	pSVR_Delta_oSOC	pSVR_Delta_DAA
	inf>>argModelParam._transData.pr_SVR_CoCirr_DeCirr;
	inf>>argModelParam._transData.pr_SVR_CoCirr_HCC;
	inf>>argModelParam._transData.pr_SVR_Delta_oSOC;
	inf>>argModelParam._transData.pr_SVR_Delta_DAA;



	return 0;
}

int Project_CEA_Prison::CEA_PSA(int argCmpIdx, int argBatch, int argBatchSize)
{
	ofstream outf_psa;
	outf_psa.open("output_project_cea_prison/PSA_"+_listCmp[argCmpIdx]+"_"+basicToStr(argBatch)+".txt");
	outf_psa<< fixed << showpoint;

	ifstream inf;
	inf.open("project_cea_prison_input_PSA.txt");
	string line;
	getline(inf,line);//skip the firstline
	
	//skip the first argBatch*1000 lines
	for(int i=0; i<argBatchSize*argBatch; i++) getline(inf,line);

	int sa_counter=1+argBatchSize*argBatch;

	int fib;
	char gender;
	double initialAge;
	for(int i=0; i<argBatchSize; i++){
		inf>>fib;
		inf>>gender;
		inf>>initialAge;
		// create cohort
		baseCohortType testCohort(fib,initialAge,gender, 'W', _listGenotype[argCmpIdx], 'N', 'A');	// first three charateristics (fib, age, gender) will be updated in the next line
		// read sampled parameter values
		ReadPSASampledValues(inf,_modelParam);
	

		vector<double> r1=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
		vector<double> r2=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);

		assert(r1.size()==r2.size());
		outf_psa<<"F"<<fib<<" "<<gender<<" "<<initialAge<<"\t"<<sa_counter<<"\t"
			<<r1[0]<<"\t"<<r2[0]<<"\t"<<r2[0]-r1[0]<<"\t"
			<<r1[1]<<"\t"<<r2[1]<<"\t"<<r2[1]-r1[1]<<"\t"
			<<(r2[1]-r1[1])/(r2[0]-r1[0])<<"\t";		
		for(int i=2;i<r1.size(); i++){
			outf_psa<<r1[i]<<"\t"<<r2[i]<<"\t";
		}
		outf_psa<<endl;

		sa_counter++;
	}

	inf.close();
	outf_psa.close();

	return 0;
}

int Project_CEA_Prison::ChangeValue(string varName, double argVal, modelParamType & argModelParam){
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


	}else{
		ExitWithMsg("[Error] Project_CEA_Prison::ChangeValue(string varName, double argVal): Unknown parameter "+varName);
	}
	return 0;
}

Project_CEA_Prison::Project_CEA_Prison()
{
	_drugCostReduction=0;
	// read all arms for comparisons
	ifstream inf;
	inf.open("project_cea_prison_input_comparators.txt");
	string word;
	while(inf>>word){
		_listCmp.push_back(word);
		inf>>word; _listArm1.push_back(word);
		inf>>word; _listArm2.push_back(word);
		int g;
		inf>>g;	_listGenotype.push_back(g);
	}
	inf.close();

}




