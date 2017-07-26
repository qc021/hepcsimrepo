#include "project_ageNoTx.h"



int Project_AgeNoTx::BaseCase_detailed()
{
	clock_t start = clock();
	/**************** SET-UPs *****************************************/
	string FILE_Basecase_AgeNoTx = "project_ageNoTx_output_base_result.txt";
	string FILE_scenario_input="project_ageNoTx_input_scenarios.txt";
	stateType listInitialState[]={s_F0,s_F1,s_F2,s_F3,s_CoCirr};
	char listGender[]={'M','F'};
	double ageSeq[]={30.0,2.0, 85.0};	// {start, interval, end}
	vector<double> listAge;
	for(double k=ageSeq[0]; k<ageSeq[2]; k=k+ageSeq[1])	listAge.push_back(k);
	//listAge.push_back(49.9);

	/******************************************************************/
	_outf.open(FILE_Basecase_AgeNoTx);
	_outf<<fixed << setprecision(4);

	_outf_detail.open("project_ageNoTx_output_detail.txt", ios::app);
	_outf_detail<<endl<<"===== "<<currentDateTime()<<" ====="<<endl<< fixed << setprecision(4);

	//// print the table header:
	//outf<<endl<<"===== "<<currentDateTime()<<" ====="<<endl<< fixed << setprecision(4);
	//outf << "TN/TE"<<"\t"<<"IFN-Tol"<<"\t"<<"Genotype"<<"\t"<<"Treatment"<<"\t"<<"InitialFibState"<<"\t"<<"Gender"<<"\t"
	//	<<"InitialAge"<<"\t"<<"ICER"<<"\t"<<"QALY_NoTX"<<"\t"<<"Cost_noTX"<<"\t"<<"QALY_TX"<<"\t"<<"Cost_TX"<<"\t"<<"Delta_QALY"<<"\t"<<"Delta_Cost"<<endl;

	ifstream inFileScenario;
	inFileScenario.open(FILE_scenario_input);
	string txExp; // treatment experience
	while(inFileScenario>>txExp){
		string ifnTol,  txArmName;
		int genotype;
		inFileScenario>>ifnTol>>genotype>>txArmName;

		for(int stateIdx=0; stateIdx<sizeof(listInitialState)/sizeof(stateType); stateIdx++){
			for(int genderIdx=0; genderIdx<sizeof(listGender)/sizeof(char); genderIdx++){
				for(int ageIdx=0; ageIdx<listAge.size(); ageIdx++){

					cout<<endl<<endl<<"===  Running scenario: "<< txExp<<"\t"<<ifnTol<<"\t"<<genotype<<"\t"<<txArmName<<"\t"<<listInitialState[stateIdx]<<"\t"<<listGender[genderIdx]<<"\t"<<listAge[ageIdx]<<"\t["<<Time(start)<<" sec]"<<endl;
					// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
					baseCohortType testCohort(listInitialState[stateIdx],listAge[ageIdx],listGender[genderIdx], 'W', genotype, 'N', 'A');

					vector<double> r=CompareTxAndNoTx(txArmName, testCohort);
					double q1=r[0];
					double c1=r[1];
					double q2=r[2];
					double c2=r[3];
					//------------------ output result ----------------------					
					_outf<<txExp<<"\t"<<ifnTol<<"\t"<<genotype<<"\t"<<txArmName<<"\t"<<listInitialState[stateIdx]<<"\t"<<listGender[genderIdx]<<"\t"<<listAge[ageIdx]<<"\t"
						<<(c2-c1)/(q2-q1)<<"\t"<<q1<<"\t"<<q2<<"\t"<<c1<<"\t"<<c2<<"\t"<<(q2-q1)<<"\t"<<(c2-c1)<<endl;
				}// end of for(int ageIdx...)
			}// end of for(int genderIdx...)
		}// end of for(int stateidx...)
	}//end of while()
	_outf.close(); _outf_detail.close(); inFileScenario.close(); 

	return 0;
}

int Project_AgeNoTx::BaseCase_aggregated(double argAgeStart, double argAgeEnd, double argAgeIncr)
{
	//if(_outputDetail){
	//	_outf_detail.open("project_ageNoTx_output_detail.txt");
	//	_outf_detail<<endl<<"===== "<<currentDateTime()<<" ====="<<endl<< fixed << setprecision(4);
	//}
	GetCutoffAge(argAgeStart, argAgeEnd,argAgeIncr,true);
	//if(_outputDetail){_outf_detail.close();}
	return 0;
}

vector<double> Project_AgeNoTx::GetCutoffAge(double argAgeStart, double argAgeEnd, double argAgeIncr, bool argOutput)
{
	string label[3]={"50K","100K","150K"};
	double benchvalue[3]={50000.0, 100000.0, 150000.0};
	int nFib=8;
	int nCut=3;

	if(argOutput){
		_outf.open("project_ageNoTx_output_ICER_age_fib_disc"+basicToStr(int(100*abs(_drugCostDisc)))+".txt");
		_outf<<fixed << setprecision(4);
		// print the table header:
		_outf << "age"<<"\t"<<"F0"<<"\t"<<"F1"<<"\t"<<"F2"<<"\t"<<"F3"<<"\t"<<"F4"<<"\t"<<"F0-F2"<<"\t"<<"F0-F3"<<"\t"<<"F3-F4"<<endl;
	}

	// [benchvalue][fibrosis]: 3 * 8
	vector<double> vec_age_cut;
	for(int k=0; k<24; k++)vec_age_cut.push_back(0);

	vector<double> icer_prev;
	int counter=0;
	while(argAgeStart+counter * argAgeIncr <= argAgeEnd){

		vector<double> icer_age_fib=EstimateICER_age(argAgeStart+counter * argAgeIncr);
		assert(icer_age_fib.size() == nFib);

		// output ICER for each age-fib
		if(argOutput){
			_outf<<argAgeStart+counter * argAgeIncr<<"\t";
			for(int k=0; k<icer_age_fib.size(); k++){
				_outf<<icer_age_fib[k]<<"\t";
			}
			_outf<<endl;
		}

		// compute age threshold
		if(counter == 0){
			for(int i=0; i<sizeof(benchvalue)/sizeof(double); i++){
				for(int j=0; j<nFib; j++){
					if(icer_age_fib[j] > benchvalue[i]){
						vec_age_cut[i*nFib + j] = argAgeStart; // the age in the last iteration, which is still cost-effective
					}
				}
			}
		}

		// for each fibrosis state, compare with benchvalues
		if(counter >= 1){
			for(int i=0; i<sizeof(benchvalue)/sizeof(double); i++){
				for(int j=0; j<nFib; j++){
					if(icer_prev[j] <= benchvalue[i] && icer_age_fib[j] > benchvalue[i]){
						vec_age_cut[i*nFib + j]=argAgeStart + (counter - 1 )*argAgeIncr; // the age in the last iteration, which is still cost-effective
					}
				}
			}
		}
		icer_prev=icer_age_fib;
		counter++;
	}
	// for those with ICER <= benchvalue in the last iteration,
	for(int i=0; i<sizeof(benchvalue)/sizeof(double); i++){
		for(int j=0; j<nFib; j++){
			if(icer_prev[j] <= benchvalue[i]){
				vec_age_cut[i*nFib + j]= argAgeStart + (counter - 1 )*argAgeIncr; // the age in the last iteration, which is still cost-effective
			}
		}
	}


	// output age threshold
	if(argOutput){
		_outf<<endl<<endl;
		for(int i=0; i<nCut; i++){
			_outf<<"CutOff"<<label[i]<<"\t";
			for(int k=0; k<nFib; k++){
				_outf<<vec_age_cut[i*nFib+k]<<"\t";
			}
			_outf<<endl;
		}		
		_outf.close();
	}
	return vec_age_cut;
}

vector<double> Project_AgeNoTx::CompareTxAndNoTx(string argTxArmName, const baseCohortType & argCohort)
{
	vector<double> rs;

	// ----------------  No Treatment arm: ----------------------
	TREATMENT = false;
	HepCSim * mySim1 = new HepCSim;
	modelParamType arm_NoTX;//(argTxArmName, argCohort);
	arm_NoTX=_modelParam;	// copy the sampled paramter values;	[9/26]
	arm_NoTX._armName=argTxArmName;	// reassign the treatment arm name;	[9/26]
	arm_NoTX._cohortData=argCohort;	// reassign the cohort data;	[9/26]
	mySim1->SetRandomSeed(SIM_SEED);
	mySim1->Run(arm_NoTX);
	//rs.push_back(mySim1->GetAvgQALY());
	//rs.push_back(mySim1->GetAvgCost());
	rs.push_back(Round(mySim1->GetAvgQALY(),4));
	rs.push_back(Round(mySim1->GetAvgCost(),4));

	FreeMem(mySim1);

	// ----------------  Treatment arm: ----------------------
	TREATMENT = true;
	HepCSim * mySim2 = new HepCSim;
	modelParamType arm_TX;//(argTxArmName, argCohort);
	arm_TX=_modelParam;	// copy the sampled paramter values;	[9/26]
	arm_TX._armName=argTxArmName;	// reassign the treatment arm name;	[9/26]
	arm_TX._cohortData=argCohort;	// reassign the cohort data;	[9/26]
	mySim2->SetRandomSeed(SIM_SEED);
	mySim2->Run(arm_TX);
	//rs.push_back(mySim2->GetAvgQALY());
	//rs.push_back(mySim2->GetAvgCost());
	rs.push_back(Round(mySim2->GetAvgQALY(),4));
	rs.push_back(Round(mySim2->GetAvgCost(),4));

	FreeMem(mySim2);




	return rs;
}



vector<double> Project_AgeNoTx::EstimateICER_age_fib(double argAge, stateType argFib)
{
	cout<<"Project_AgeNoTx::EstimateICER_age_fib: age "<<argAge<<", fibrosis stage "<<argFib<<", drug cost discount "<<_drugCostDisc<<endl;

	double icer=0;
	double q1=0;
	double q2=0;
	double c1=0;
	double c2=0;
	double sumWt=0;


	char listGender[]={'M','F'};

	ifstream inFileScenario;
	inFileScenario.open(FILE_scenario_input);

	string txExp; // treatment experience
	while(inFileScenario>>txExp){
		string ifnTol,  txArmName;
		int genotype;
		inFileScenario>>ifnTol>>genotype>>txArmName;

		//cout<<"Estimate ICER (Tx vs NoTx): "<<argFib<<" "<<txExp<<" "<<ifnTol<<" G"<<genotype<<" "<<txArmName<<endl;
		for(int genderIdx=0; genderIdx<sizeof(listGender)/sizeof(char); genderIdx++){

			// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
			baseCohortType testCohort(argFib,argAge,listGender[genderIdx], 'W', genotype, 'N', 'A');
			vector<double> r;

			if(txExp=="TN" && genotype == 1 && "TN_G1_LDP_SOF" == txArmName){
				vector<double> r_8week=CompareTxAndNoTx("TN_G1_LDP_SOF8", testCohort);
				vector<double> r_12week=CompareTxAndNoTx("TN_G1_LDP_SOF12", testCohort);
				assert(r_8week.size()==r_12week.size());
				for(int l=0; l<r_8week.size(); l++){
					r.push_back(_p8Week*r_8week[l] + (1-_p8Week)*r_12week[l]);
				}
			}else{
				r=CompareTxAndNoTx(txArmName, testCohort);
			}

			// output the detailed ICER for each scenario
			// 
			if(_outputDetail){
				_outf_detail.open("project_ageNoTx_output_detail.txt",ios::app);
				_outf<<fixed<<setprecision(6);
				_outf_detail<<argFib<<"\t"<<argAge<<"\t"<<listGender[genderIdx]<<"\t"<<genotype<<"\t"
					<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<"\t"<<r[3]<<"\t"<<(r[3]-r[1])/(r[2]-r[0])<<endl;
				_outf_detail.close();
			}

			double wt=GetWeight(txExp,ifnTol,genotype,listGender[genderIdx]);

			q1=q1+r[0] * wt;
			c1=c1+r[1] * wt;
			q2=q2+r[2] * wt;
			c2=c2+r[3] * wt;

			sumWt = sumWt + wt;
		}// end of for(int genderIdx...)
	}//end of while()
	assert(abs(sumWt-1)<EPSILON);
	inFileScenario.close(); 
	icer=(c2-c1)/(q2-q1);

	vector<double> vecRs;
	vecRs.push_back(q1);
	vecRs.push_back(c1);
	vecRs.push_back(q2);
	vecRs.push_back(c2);
	vecRs.push_back(icer);

	return vecRs;
}

vector<double> Project_AgeNoTx::EstimateICER_age(double argAge)
{


	vector<double> vecICER;
	// in following rs_XX:  double[5]={q1, c1, q2, c2, ICER}
	vector<double> rs_f0=EstimateICER_age_fib(argAge,s_F0);
	vector<double> rs_f1=EstimateICER_age_fib(argAge,s_F1);
	vector<double> rs_f2=EstimateICER_age_fib(argAge,s_F2);
	vector<double> rs_f3=EstimateICER_age_fib(argAge,s_F3);
	vector<double> rs_f4=EstimateICER_age_fib(argAge,s_CoCirr);
	vecICER.push_back(rs_f0[4]);
	vecICER.push_back(rs_f1[4]);
	vecICER.push_back(rs_f2[4]);
	vecICER.push_back(rs_f3[4]);
	vecICER.push_back(rs_f4[4]);

	vector<vector<double>> mtxQandC = ZeroVec2d(4,5); //[q1,c1,q2,c2] x [f0,f1,...,f4]
	for(int k=0; k<4; k++){
		mtxQandC[k][0]=rs_f0[k];
		mtxQandC[k][1]=rs_f1[k];
		mtxQandC[k][2]=rs_f2[k];
		mtxQandC[k][3]=rs_f3[k];
		mtxQandC[k][4]=rs_f4[k];
	}
	double q1, q2, c1, c2;
	int beg,end;
	// F0-F2
	beg=0;
	end=2;
	q1=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[0],beg,end));
	c1=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[1],beg,end));
	q2=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[2],beg,end));
	c2=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[3],beg,end));
	vecICER.push_back((c2-c1)/(q2-q1));
	//cout<<"============"<<q1<<"\t"<<c1<<"\t"<<q2<<"\t"<<c2<<endl;
	// F0-F3
	beg=0;
	end=3;
	q1=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[0],beg,end));
	c1=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[1],beg,end));
	q2=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[2],beg,end));
	c2=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[3],beg,end));
	vecICER.push_back((c2-c1)/(q2-q1));
	//cout<<"============"<<q1<<"\t"<<c1<<"\t"<<q2<<"\t"<<c2<<endl;
	// F3-F4
	beg=3;
	end=4;
	q1=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[0],beg,end));
	c1=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[1],beg,end));
	q2=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[2],beg,end));
	c2=InnerProd(Renormalize(_fibDistr,beg,end),SubVector(mtxQandC[3],beg,end));
	vecICER.push_back((c2-c1)/(q2-q1));
	//cout<<"============"<<q1<<"\t"<<c1<<"\t"<<q2<<"\t"<<c2<<endl;

	return vecICER;
}

Project_AgeNoTx::Project_AgeNoTx()
{
	_fibDistr.push_back(0.0975);
	_fibDistr.push_back(0.2543);
	_fibDistr.push_back(0.2066);
	_fibDistr.push_back(0.1960);
	_fibDistr.push_back(0.2457);

	_genoDistr.push_back(0.796);
	_genoDistr.push_back(0.13);
	_genoDistr.push_back(0.063);
	_genoDistr.push_back(0.011);

	_pMale=0.64;
	_pIFN_intol=0.23;
	_pTE=0.39;
	_p8Week=0.57;

	FILE_scenario_input="project_ageNoTx_input_scenarios.txt";
	_outputDetail=false;
}

double Project_AgeNoTx::GetWeight(string argTxExp, string argIFNTol, int argGenotype, char argGender)
{
	double w=1;

	// genotype: 1,2,3,4
	w=w*_genoDistr[argGenotype-1];
	// gender
	if(argGender == 'M'){
		w=w*_pMale;
	}else if (argGender == 'F'){
		w=w*(1-_pMale);
	}else{
		ExitWithMsg("Error: GetWeight(): unknown argGender");
	}
	// TE / TN
	if(argTxExp == "TE"){
		w=w*_pTE;
	}else if (argTxExp == "TN"){
		if(argIFNTol == "TOL"){
			w=w * (1-_pTE) * (1-_pIFN_intol);
		}else if (argIFNTol == "NOT"){
			w=w * (1-_pTE) * _pIFN_intol;
		}else{
			ExitWithMsg("Error: GetWeight(): unknown argTxExp: "+argIFNTol);
		}
	}else{
		ExitWithMsg("Error: GetWeight(): unknown argTxExp: "+argTxExp);
	}


	return w;
}

int Project_AgeNoTx::PSA_ICER_for_each_age(double argAge,int argBatchIdx, long argSeed)
{
	if(argBatchIdx != -1 && argSeed != -1 ){
		NUM_RUNS_PSA=NUM_RUNS_PSA_ONE_BATCH;
		PSA_SEED=argSeed;
	}

	string label[3]={"50K","100K","150K"};
	double benchvalue[3]={50000.0, 100000.0, 150000.0};
	double icer_conf[3][8]={0.0};	//[icer standard][fibrosis state]



	cout<<endl<<"---- PSA_ICER for age "<<argAge<<endl;
	psaDistrType psaDistr;
	psaDistr.ReadDistrForPSA("SA_PSA_distr.in");

	HepCSim mySim_psaSampler;		//used for sampling
	mySim_psaSampler.PSA_initializeSampler(psaDistr);

	for(int n=0; n<NUM_RUNS_PSA; n++){
		if((n+1) % (NUM_RUNS_PSA/10) == 0 ){
			cout<<(n+1)<<" PSA iterations finished...("<<(n+1)*1.0 / (1.0 * NUM_RUNS_PSA) * 100<<"%)"<<endl;
		}
		// sample
		mySim_psaSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
		// solve
		vector<double> icer_age_fib=EstimateICER_age(argAge);

		ofstream icerOut;
		icerOut.open("output_PSA_tx_for_age/ICER_samples_age_"+basicToStr(argAge)+".txt", ios::app);
		icerOut<<fixed<<setprecision(3);
		int bidx= (argBatchIdx==-1?0:argBatchIdx);
		icerOut<<bidx * NUM_RUNS_PSA_ONE_BATCH + n<<"\t"<<fixed<<setprecision(2);
		for(int idx=0; idx<icer_age_fib.size(); idx++){
			icerOut<<icer_age_fib[idx]<<"\t";
		}
		icerOut<<endl;
		icerOut.close();


		// count the confidence (i.e., percentage)
		for(int i=0; i<sizeof(benchvalue)/sizeof(double); i++){
			for(int j=0; j<sizeof(icer_conf[i])/sizeof(double); j++){
				if(icer_age_fib[j]<=benchvalue[i]){
					icer_conf[i][j] = icer_conf[i][j] + 1.0/NUM_RUNS_PSA;
				}
			}
		}
	}
	for(int i=0; i<sizeof(benchvalue)/sizeof(double); i++){
		// output one line in output-file:
		_psaOutf_icer_conf.open("output_PSA_tx_for_age/ICER_conf_by_age_fib_"+label[i]+".txt",ios::app);
		_psaOutf_icer_conf<<fixed<<setprecision(6);
		_psaOutf_icer_conf<<argAge<<"\t";
		for(int j=0; j<sizeof(icer_conf[i])/sizeof(double); j++){
			_psaOutf_icer_conf<<icer_conf[i][j]<<"\t";
		}
		_psaOutf_icer_conf<<endl;
		_psaOutf_icer_conf.close();

	}

	return 0;
}




int Project_AgeNoTx::PSA_age_threshold(int argNumPSARunsPerBatch, int idx,double argAgeStart, double argAgeEnd, double argAgeIncr)
{
	PSA_SEED=3*(2+idx);

	_psaOutf_threshold.open("output_PSA_tx_for_age/age_threhold_"+basicToStr(idx)+".txt");


	psaDistrType psaDistr;
	psaDistr.ReadDistrForPSA("SA_PSA_distr.in");

	HepCSim mySim_psaSampler;		//used for sampling
	mySim_psaSampler.PSA_initializeSampler(psaDistr);

	for(int n=0; n<argNumPSARunsPerBatch; n++){
		if((n+1) % (argNumPSARunsPerBatch/10) == 0 ){
			cout<<(n+1)<<" PSA iterations finished...("<<(n+1)*1.0 / (1.0 * argNumPSARunsPerBatch) * 100<<"%)"<<endl;
		}

		// sample
		mySim_psaSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
		// solve
		vector<double> thre=GetCutoffAge(argAgeStart, argAgeEnd, argAgeIncr);

		// output one line for this PSA sample
		_psaOutf_threshold<< idx*argNumPSARunsPerBatch + n<<"\t"<<fixed<<setprecision(1);
		for(int j=0; j<thre.size(); j++){
			_psaOutf_threshold<<thre[j]<<"\t";
		}
		_psaOutf_threshold<<endl;

	}
	_psaOutf_threshold.close();
	return 0;
}


int Project_AgeNoTx::ReduceDrugCost( double reductPerct )
{
	_drugCostDisc=reductPerct;
	_modelParam=_modelParam_backup;

	_modelParam._costData.c_BOC=_modelParam._costData.c_BOC * (1 + reductPerct);
	_modelParam._costData.c_RBV=_modelParam._costData.c_RBV * (1 + reductPerct);
	_modelParam._costData.c_PEG=_modelParam._costData.c_PEG * (1 + reductPerct);
	_modelParam._costData.c_SMV=_modelParam._costData.c_SMV * (1 + reductPerct);
	_modelParam._costData.c_SOF=_modelParam._costData.c_SOF * (1 + reductPerct);
	_modelParam._costData.c_TEL=_modelParam._costData.c_TEL * (1 + reductPerct);
	_modelParam._costData.c_LDV=_modelParam._costData.c_LDV * (1 + reductPerct);
	_modelParam._costData.c_DCV=_modelParam._costData.c_DCV * (1 + reductPerct);
	_modelParam._costData.c_PrOD=_modelParam._costData.c_PrOD * (1 + reductPerct);

	return 0;

}

int Project_AgeNoTx::DSA_OneWay_ICER_for_age(double argAge)
{
	// backup the model parameters;
	modelParamType backupValue=_modelParam;


	ofstream outfSA;
	outfSA.open("project_ageNoTx_output_SA_ICER_disc"+basicToStr(int(abs(_drugCostDisc)*100))+"_age_"+basicToStr(argAge)+".txt");
	outfSA<<fixed<<setprecision(3);
	map<string,SParamTriplet> ranges=ReadOneWaySAParamValueRange("SA_OneWay_param_ranges.txt");

	// output baseline value
	outfSA<<"Baseline\t\t\t\t";

	vector<double> icer_age_fib_base=EstimateICER_age(argAge);
	for(int k=0; k<icer_age_fib_base.size(); k++){
		outfSA<<icer_age_fib_base[k]<<"\t";
	}
	outfSA<<endl;

	// oneway sensitivity analysis
	for(map<string,SParamTriplet>::const_iterator it=ranges.begin(); it!=ranges.end(); it++){
		outfSA<<it->first<<"\t"<<it->second._base<<"\t"<<it->second._lb<<"\t"<<it->second._ub<<"\t";
		// try low value
		_modelParam=backupValue;
		DSA_OneWay_ChangeValue(it->first, it->second._lb); //update the _modelParam
		vector<double> icer_age_fib_low=EstimateICER_age(argAge);
		for(int k=0; k<icer_age_fib_low.size(); k++){
			outfSA<<icer_age_fib_low[k]<<"\t";
		}

		// try high value
		_modelParam=backupValue;
		DSA_OneWay_ChangeValue(it->first, it->second._ub); //update the _modelParam
		vector<double> icer_age_fib_high=EstimateICER_age(argAge);
		for(int k=0; k<icer_age_fib_high.size(); k++){
			outfSA<<icer_age_fib_high[k]<<"\t";
		}
		outfSA<<endl;

	}
	outfSA.close();
	_modelParam=backupValue;//restore the nominal values
	return 0;
}

int Project_AgeNoTx::DSA_OneWay_ChangeValue(string varName, double argVal){
	if ("q_F0"== varName){_modelParam._qolData.q_F0 = argVal;
	}else if ("q_F1"== varName){_modelParam._qolData.q_F1 = argVal;
	}else if ("q_F2"== varName){_modelParam._qolData.q_F2 = argVal;
	}else if ("q_F3"== varName){_modelParam._qolData.q_F3 = argVal;
	}else if ("q_CoCirr"== varName){_modelParam._qolData.q_CoCirr = argVal;
	}else if ("q_DeCirr"== varName){_modelParam._qolData.q_DeCirr = argVal;
	}else if ("q_HCC"== varName){_modelParam._qolData.q_HCC = argVal;
	}else if ("q_LivTr"== varName){_modelParam._qolData.q_LivTr = argVal;
	}else if ("q_PostLivT"== varName){_modelParam._qolData.q_PostLivTr = argVal;
	}else if ("q_SVR"== varName){_modelParam._qolData.q_SVR = argVal;
	}else if ("q_Anemia"== varName){_modelParam._qolData.q_Dec_Anemia = argVal;
	}else if ("q_TX_oSOC"== varName){_modelParam._qolData.q_TX_oSOC = argVal;
	}else if ("q_Tx_DAA"== varName){_modelParam._qolData.q_TX_DAA = argVal;
	}else if ("c_F0"== varName){_modelParam._costData.c_F0 = argVal;
	}else if ("c_F1"== varName){_modelParam._costData.c_F1 = argVal;
	}else if ("c_F2"== varName){_modelParam._costData.c_F2 = argVal;
	}else if ("c_F3"== varName){_modelParam._costData.c_F3 = argVal;
	}else if ("c_CoCirr"== varName){_modelParam._costData.c_CoCirr = argVal;
	}else if ("c_DeCirr"== varName){_modelParam._costData.c_DeCirr = argVal;
	}else if ("c_DeCirr1yrPlus"== varName){_modelParam._costData.c_DeCirr1yrPlus = argVal;
	}else if ("c_HCC"== varName){_modelParam._costData.c_HCC = argVal;
	}else if ("c_LivTr"== varName){_modelParam._costData.c_LivTr = argVal;
	}else if ("c_PostLivTr"== varName){_modelParam._costData.c_PostLivTr = argVal;
	}else if ("pF0_F1_SA"== varName){_modelParam._transData.pr_F0_F1 = argVal;
	}else if ("pF1_F2_SA"== varName){_modelParam._transData.pr_F1_F2 = argVal;
	}else if ("pF2_F3_SA"== varName){_modelParam._transData.pr_F2_F3 = argVal;
	}else if ("pF3_F4_SA"== varName){_modelParam._transData.pr_F3_CoCirr = argVal;
	}else if ("pF4_DC_SA"== varName){_modelParam._transData.pr_CoCirr_DeCirr = argVal;
	}else if ("pF4_HCC_SA"== varName){_modelParam._transData.pr_CoCirr_HCC = argVal;
	}else if ("pDC_HCC_SA"== varName){_modelParam._transData.pr_DeCirr_HCC = argVal;
	}else if ("pDC_Liv_Transpl_SA"== varName){_modelParam._transData.pr_DeCirr_LivTr = argVal;
	}else if ("pMort_dc_cyc_1_SA"== varName){_modelParam._transData.pr_DeCirr_DeathLiv = argVal;
	}else if ("pMort_dc_cyc_2_SA"== varName){_modelParam._transData.pr_DeCirr1yrPlus_DeathLiv = argVal;
	}else if ("pHCC_Liv_Transpl_SA"== varName){_modelParam._transData.pr_HCC_LivTr = argVal;
	}else if ("pMort_hcc_cyc_SA"== varName){_modelParam._transData.pr_HCC_DeathLiv = argVal;
	}else if ("pMort_Liv_Transpl_SA"== varName){_modelParam._transData.pr_LivTr_DeathLiv = argVal;
	}else if ("pMort_Post_Liv_Transpl_SA"== varName){_modelParam._transData.pr_LivTr1yrPlus_DeathLiv = argVal;
	}else if ("pr_SVR_CoCirr_DeCirr"== varName){_modelParam._transData.pr_SVR_CoCirr_DeCirr = argVal;
	}else if ("pr_SVR_CoCirr_HCC"== varName){_modelParam._transData.pr_SVR_CoCirr_HCC = argVal;
	}else if ("pSVR_Delta_oSOC"== varName){_modelParam._transData.pr_SVR_Delta_oSOC = argVal;
	}else if ("pSVR_Delta_DAA"== varName){_modelParam._transData.pr_SVR_Delta_DAA = argVal;


	}else{
		ExitWithMsg("[Error] Project_AgeNoTx::DSA_OneWay_ChangeValue(string varName, double argVal): Unknown parameter "+varName);
	}
	return 0;
}

int Project_AgeNoTx::BaseCase_distributed( double argAgeStart,double argAgeEnd, double argAgeIncr )
{
	int nFib=8;
	_outf.open("project_ageNoTx_output_ICER_age_fib_distr_"+basicToStr(int(abs(_drugCostDisc)*100))+".txt",ios::app);
	_outf<<fixed << setprecision(4);
	int counter=0;
	while(argAgeStart+counter * argAgeIncr < argAgeEnd){

		vector<double> icer_age_fib=EstimateICER_age(argAgeStart+counter * argAgeIncr);
		assert(icer_age_fib.size() == nFib);

		// output ICER for each age-fib
		_outf<<argAgeStart+counter * argAgeIncr<<"\t";
		for(int k=0; k<icer_age_fib.size(); k++){
			_outf<<icer_age_fib[k]<<"\t";
		}
		_outf<<endl;
		counter++;
	}
	_outf.close();
	return 0;
}





