#include "project_tfs.h"
using namespace std;




int Project_TFS::BaseCase(){
	ifstream inFileLHS;
	ofstream outFileSummary;


	outFileSummary.open("project_tfs_output_basecase.txt");
	outFileSummary<<endl<<"===== "<<currentDateTime()<<" ====="<<endl
		<< fixed << showpoint;
	// print the table header:
	outFileSummary << "Scenario" << "\t" << "Genotype" << "\t" << "Treat"<<"\t"
		<< "QALY" << "\t" << "COST" << "\t" << "COUNT_DECIRR" << "\t"
		<< "COUNT_HCC" << "\t" << "COUNT_LIVTR" << "\t" << "COUNT_DEATHLIV" 
		<< "\t" << "COST_TX" << endl;

	ofstream oTfs;
	oTfs.open("project_tfs_output_transpl_free_survival.txt");
	oTfs<<"fibState\tgender\tgenotype\tcounter\t";
	int horz = (int)(15 * CYCLES_PER_YEAR);
	for(int k=0; k<horz; k++) oTfs<<k+1<<"\t";
	oTfs<<endl;

	// ==== prepare the model parameters ====
	vector<modelParamType> vecArms;
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios=10;
	stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F','M','F','M','F'};
	//double listAge[]={49.9,49.9,56.2,56.2,57.6,57.6,58.3,58.3,58.7,58.7};
	double listAge[]={54.5,54.5,54.5,54.5,54.5,54.5,54.5,54.5,54.5,54.5};

	// === arms ===
	string listTxArm[]={"NoTx","NoTx"};
	int listGenotype[]={1,2};


	for(int a=0; a<sizeof(listGenotype)/sizeof(int); a++){

		for(int s=0; s<nScenarios; s++){
			
			// define cohort
			// create patient profile:
			// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
			baseCohortType testCohort(listState[s],listAge[s],listGender[s], 'W', listGenotype[a], 'N', 'A');	

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

			cout<<armX._armName<<"\t"<<listState[s]<<"\t"<<listGender[s]<<"\t"<<listGenotype[a]<<"\t"<<TREATMENT<<"..."<<endl;

			HepCSim mySim;
			mySim.SetRandomSeed(SIM_SEED);
			mySim.Run(armX);
			// output one line summary of this arm
			outFileSummary<<armX._armName<<"\t"<<listGenotype[a]<<"\t"<<TREATMENT<<"\t";
			mySim.OutputBaseCase(outFileSummary);

			
		// 			prevStateCountF0.push_back(0);
		// prevStateCountF1.push_back(0);
		// prevStateCountF2.push_back(0);
		// prevStateCountF3.push_back(0);
		// prevStateCountCoCirr.push_back(0);
		// prevStateCountDeCirr.push_back(0);
		// prevStateCountDeCirr1yrPlus.push_back(0);
		// prevStateCountHCC.push_back(0);

			vector<int> vF0=mySim.GetCounter().prevStateCountF0;
			vector<int> vF1=mySim.GetCounter().prevStateCountF1;
			vector<int> vF2=mySim.GetCounter().prevStateCountF2;
			vector<int> vF3=mySim.GetCounter().prevStateCountF3;
			vector<int> vF4=mySim.GetCounter().prevStateCountCoCirr;
			vector<int> vDC1=mySim.GetCounter().prevStateCountDeCirr;
			vector<int> vDC2=mySim.GetCounter().prevStateCountDeCirr1yrPlus;
			vector<int> vHCC=mySim.GetCounter().prevStateCountHCC;
			

			oTfs<<listState[s]<<"\t"<<listGender[s]<<"\t"<<listGenotype[a]<<"\t"<<"atRisk\t";
			for(int k=0; k<horz; k++) oTfs<<vF0[k]+vF1[k]+vF2[k]+vF3[k]+vF4[k]+vDC1[k]+vDC2[k]+vHCC[k]<<"\t";
			oTfs<<endl;

			vector<int> vTr=mySim.GetCounter().incidentStateCountLivTr;
			oTfs<<listState[s]<<"\t"<<listGender[s]<<"\t"<<listGenotype[a]<<"\t"<<"trspl\t";
			for(int k=0; k<horz; k++) oTfs<<vTr[k]<<"\t";
			oTfs<<endl;

			vector<int> vTr1plus=mySim.GetCounter().incidentStateCountLivTr1yrPlus;
			oTfs<<listState[s]<<"\t"<<listGender[s]<<"\t"<<listGenotype[a]<<"\t"<<"tr1plus\t";
			for(int k=0; k<horz; k++) oTfs<<vTr1plus[k]<<"\t";
			oTfs<<endl;
			
			

		}


	}



	outFileSummary.close();
	return 0;
}





int Project_TFS::GetTFS()
{
	//TN TOL
	//Compare("TN_TOL_G1_BOC","TN_TOL_G1_SOF_PEG_RBV12",1);
	//Compare("TN_TOL_G1_TEL","TN_TOL_G1_SOF_PEG_RBV12",1);

	//Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF8",1);
	//Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF12",1);
	Compare("TN_TOL_G1_BOC","TN_G1_LDP_SOF",1);
	//Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF8",1);
	//Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF12",1);
	Compare("TN_TOL_G1_TEL","TN_G1_LDP_SOF",1);
	Compare("TN_TOL_G2_PEG_RBV24","TN_TOL_G2_SOF_RBV12",2);
	Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_RBV24",3);
	Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_PEG_RBV12",4);

	// TN NOT-TOL
	//Compare("NoTx","TN_NOT_G1_SOF_SMV12",1);

	//Compare("NoTx","TN_G1_LDP_SOF8",1);
	//Compare("NoTx","TN_G1_LDP_SOF12",1);
	Compare("NoTx","TN_G1_LDP_SOF",1);
	Compare("NoTx","TN_NOT_G2_SOF_RBV12",2);
	Compare("NoTx","TN_NOT_G3_SOF_RBV24",3);
	Compare("NoTx","TN_NOT_G4_SOF_RBV24",4);

	// TE
	//Compare("TE_G1_BOC","TE_G1_SOF_SMV12",1);
	//Compare("TE_G1_TEL","TE_G1_SOF_SMV12",1);

	Compare("TE_G1_BOC","TE_G1_LDP_SOF12_24",1);
	Compare("TE_G1_TEL","TE_G1_LDP_SOF12_24",1);
	Compare("TE_G2_PEG_RBV24","TE_G2_SOF_RBV12",2);
	Compare("TE_G3_PEG_RBV24","TE_G3_SOF_RBV24",3);
	Compare("TE_G4_PEG_RBV48","TE_G4_SOF_PEG_RBV12",4);

	return 0;
}

int Project_TFS::Compare(string strArm1, string strArm2, int argGenotype)
{


	_outFileSummary.open("project_doublecheck_output_compare.txt",ios::app);
	_outFileSummary<<endl<<"===== "<<currentDateTime()<<" ===== "<<strArm1<<" vs "<<strArm2<<" ==="<<endl
		<< fixed << showpoint;
	// print the table header:



	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios=10;
	stateType listState[]={s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr};
	char listGender[]={'M','F','M','F','M','F','M','F','M','F'};
	double listAge[]={49.9,49.9,56.2,56.2,57.6,57.6,58.3,58.3,58.7,58.7};

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

vector<double> Project_TFS::EvaluateOneArm( string strArm, const baseCohortType & testCohort )
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

vector<double> Project_TFS::EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort )
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


Project_TFS::Project_TFS()
{
	_drugCostReduction=0;
	// read all arms for comparisons
	ifstream inf;
	inf.open("project_doublecheck_input_comparators.txt");
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




