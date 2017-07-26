#include"project_priortzTx.h"


Project_priortzTx::Project_priortzTx()
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

}

vector<double> Project_priortzTx::GetLossInQALYs(double argAge, stateType argFib,txHistryTolType argTxHistrTol,char argGender, string argArmName)
{
	double q_noTX, q_TxNow, q_delay;

	rs_debug0.clear();
	rs_debug1.clear();
	rs_debug2.clear();
	vector<double> rs;

	// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
	baseCohortType testCohort(argFib, argAge, argGender, 'W', 1, 'N', 'A');	
	modelParamType modelParam;//(argTxArmName, argCohort);
	modelParam=_modelParam;	// copy the sampled paramter values;	[9/30]
	modelParam._cohortData=testCohort;
	modelParam._armName=argArmName;//GetArmName(argFib,argTxHistrTol);	// reassign the treatment arm name;	[9/30]




	//// ============== No treatment ==============================================================
	TREATMENT=false;
	HepCSim * mySim0 = new HepCSim;
	mySim0->SetRandomSeed(SIM_SEED);	
	mySim0->Run(modelParam);
	//rs.push_back(mySim1->GetAvgQALY());
	//rs.push_back(mySim1->GetAvgCost());	
	q_noTX=Round(mySim0->GetAvgQALY(),4);
	rs.push_back(q_noTX);

	
	FreeMem(mySim0);


	TREATMENT=true;
	// ============== treat immediately ==============================================================
	TREATMENT=true;
	HepCSim * mySim1 = new HepCSim;
	mySim1->SetRandomSeed(SIM_SEED);	
	mySim1->Run(modelParam);
	//rs.push_back(mySim1->GetAvgQALY());
	//rs.push_back(mySim1->GetAvgCost());	
	q_TxNow=Round(mySim1->GetAvgQALY(),4);
	rs.push_back(q_TxNow);
	
	FreeMem(mySim1);

	// ============== Delay treatment for one year ======================================================
	TREATMENT=true;
	HepCSim * mySim2 = new HepCSim;
	mySim2->SetRandomSeed(SIM_SEED);
	mySim2->SetTxDelayedWeeks(_delayedWeeks);	// delay treatment for 1 year (~ 52 weeks)
	mySim2->Run(modelParam);
	//rs.push_back(mySim1->GetAvgQALY());
	//rs.push_back(mySim1->GetAvgCost());
	q_delay=Round(mySim2->GetAvgQALY(),4);
	rs.push_back(q_delay);
	
	FreeMem(mySim2);




	rs.push_back(q_TxNow-q_delay);
	return rs;

}


vector<double> Project_priortzTx::GetLossInQALYs_Combine8and12Weeks( double argAge, stateType argFib,txHistryTolType argTxHistrTol,char argGender )
{
	
	if(argFib!=s_CoCirr && argTxHistrTol != TE_NA){
		vector<double> r;
		vector<double> r_8week=GetLossInQALYs(argAge,argFib,argTxHistrTol,argGender,"TN_G1_LDP_SOF8");
		vector<double> r_12week=GetLossInQALYs(argAge,argFib,argTxHistrTol,argGender,"TN_G1_LDP_SOF12");
		r=LinearCombTwoVec(r_8week,r_12week,_p8Week);
		return r;
	}else{
		return GetLossInQALYs(argAge,argFib,argTxHistrTol,argGender,GetArmName(argFib,argTxHistrTol));
	}
}



std::string Project_priortzTx::GetArmName(stateType argFib, txHistryTolType argTxHistr)
{
	string s;

	if(argFib != s_CoCirr && argTxHistr == TN_TOL){
		s="TN_G1_LDP_SOF12";
	}else if(argFib != s_CoCirr && argTxHistr == TN_INT){
		s="TN_G1_LDP_SOF12";
	}else if(argFib != s_CoCirr && argTxHistr == TE_NA){
		s="TE_G1_LDP_SOF12_24";
	}else if(argFib == s_CoCirr && argTxHistr == TN_TOL){
		s="TN_G1_LDP_SOF12";
	}else if(argFib == s_CoCirr && argTxHistr == TN_INT){
		s="TN_G1_LDP_SOF12";
	}else if(argFib == s_CoCirr && argTxHistr == TE_NA){
		s="TE_G1_LDP_SOF12_24";
	}else{
		ExitWithMsg("[Error]: GetArmName(): unknown argFib or argTxHistry");
	}

	return s;
}

int Project_priortzTx::BaseCase()
{

	int listDelayedWeeks[]={52};
	double listAge[]={25,30,35,40,45,50,55,60,65,70,75};
	stateType listFib[]={s_F0,s_F1,s_F2,s_F3,s_CoCirr};
	txHistryTolType listTxHstrTol[]={TN_TOL,TE_NA};
	char listGender[]={'M','F'};

	//int listDelayedWeeks[]={52};
	////double listAge[]={55+52/CYCLES_PER_YEAR,60+52/CYCLES_PER_YEAR,65+52/CYCLES_PER_YEAR,70+52/CYCLES_PER_YEAR,75+52/CYCLES_PER_YEAR};
	//double listAge[]={55,60,65,70,75};
	//stateType listFib[]={s_CoCirr}; //s_F0,s_F1,s_F2,s_F3,
	////stateType listFib[]={s_CoCirr}; 
	//txHistryTolType listTxHstrTol[]={TE_NA};
	//char listGender[]={'M'};

	_outf.open("project_priortzTx_output_basecase.txt");
	_outf<<"age\tfibState\tTxHistory\tGender\tDelayedWeeks\tQ_NoTx\tQ_TxNow\tQ_Delayed\tLossInQALYs"<<endl;
	for(int a=0; a<sizeof(listAge)/sizeof(double); a++){
		for(int b=0; b<sizeof(listFib)/sizeof(stateType); b++){
			for(int c=0; c<sizeof(listTxHstrTol)/sizeof(txHistryType); c++){
				for(int d=0; d<sizeof(listGender)/sizeof(char); d++){
					for(int f=0; f<sizeof(listDelayedWeeks)/sizeof(int); f++){
						_outf<<listAge[a]<<"\t"
							<<listFib[b]<<"\t"
							<<listTxHstrTol[c]<<"\t"
							<<listGender[d]<<"\t"
							<<listDelayedWeeks[f]<<"\t";


						_delayedWeeks=listDelayedWeeks[f];
						vector<double> r;
						

						//r=GetLossInQALYs(listAge[a],listFib[b],listTxHstrTol[c],listGender[d],
						//								GetArmName(listFib[b],listTxHstrTol[c]));


						// Note: GetLossInQALYs_Combined8and12Weeks: rs_debug[] are the results only for 12-week case
						r=GetLossInQALYs_Combine8and12Weeks(listAge[a],listFib[b],listTxHstrTol[c],listGender[d]);
						
						for(int k=0; k<r.size(); k++){
							_outf<<r[k]<<"\t";
						}
						//for(int k=0; k<rs_debug0.size(); k++){
						//	_outf<<rs_debug0[k]<<"\t";
						//}
						//for(int k=0; k<rs_debug1.size(); k++){
						//	_outf<<rs_debug1[k]<<"\t";
						//}
						for(int k=0; k<rs_debug2.size(); k++){
							_outf<<rs_debug2[k]<<"\t";
						}
						_outf<<endl;
					}
				}
			}
		}
	}

	_outf.close();
	return 0;
}


int Project_priortzTx::OneWaySA(double argAge, int argDelayedWeeks, txHistryTolType argHistryTol, stateType argFib)
{
	//_baselineModelParam=_modelParam; // QC: no need to have this as _baselineModelParam is initialized with default values
	double theAge=argAge;
	int delayedWeeks=argDelayedWeeks;	
	txHistryTolType theHistryTol=argHistryTol;

	stateType thefibState=argFib;


	ofstream outf_onewaySA;
	outf_onewaySA.open("project_priortzTx_output_SA_"
		+basicToStr(theAge)+"_F"+basicToStr(thefibState)+"_TxHis"+basicToStr(theHistryTol)
		+".txt");
	outf_onewaySA<< fixed << showpoint;

	vector<string> sa_varName;
	vector<SParamTriplet> sa_varVal;
	ReadOneWaySAParamValueRange("SA_OneWay_param_ranges.txt",sa_varName,sa_varVal);
	_delayedWeeks=delayedWeeks;


	// === Get baseline result for comparison ====
	_modelParam=_baselineModelParam;// recover the baseline values
	vector<double> r_male=GetLossInQALYs_Combine8and12Weeks(theAge,thefibState,theHistryTol,'M');
	vector<double> r_female=GetLossInQALYs_Combine8and12Weeks(theAge,thefibState,theHistryTol,'F');		

	assert(r_male.size()==r_female.size());
	outf_onewaySA<<"\t\tbaseline\t";
	for(int l=0; l<r_male.size(); l++){outf_onewaySA<<r_male[l]*_pMale + r_female[l]*(1-_pMale)<<"\t";}
	for(int l=0; l<r_male.size(); l++){outf_onewaySA<<r_male[l]*_pMale + r_female[l]*(1-_pMale)<<"\t";}// repeat once
	outf_onewaySA<<endl;

	// === start one-way SA ====
	for(int j=0; j<sa_varName.size(); j++){

		outf_onewaySA<<sa_varVal[j]._lb<<"\t"<<sa_varVal[j]._ub<<"\t"<<sa_varName[j]<<"\t";

		// ===== low value ======
		_modelParam=_baselineModelParam;// recover the baseline values
		ChangeValue(sa_varName[j],sa_varVal[j]._lb,_modelParam);

		vector<double> r_m_low=GetLossInQALYs_Combine8and12Weeks(theAge,thefibState,theHistryTol,'M');
		vector<double> r_f_low=GetLossInQALYs_Combine8and12Weeks(theAge,thefibState,theHistryTol,'F');		
		assert(r_m_low.size()==r_f_low.size());
		for(int l=0; l<r_m_low.size(); l++){outf_onewaySA<<r_m_low[l]*_pMale + r_f_low[l]*(1-_pMale)<<"\t";}

		// ===== high value ======
		_modelParam=_baselineModelParam;// recover the baseline values
		ChangeValue(sa_varName[j],sa_varVal[j]._ub,_modelParam);

		vector<double> r_m_high=GetLossInQALYs_Combine8and12Weeks(theAge,thefibState,theHistryTol,'M');
		vector<double> r_f_high=GetLossInQALYs_Combine8and12Weeks(theAge,thefibState,theHistryTol,'F');		
		assert(r_m_high.size()==r_f_high.size());
		for(int l=0; l<r_m_high.size(); l++){outf_onewaySA<<r_m_high[l]*_pMale + r_f_high[l]*(1-_pMale)<<"\t";}

		outf_onewaySA<<endl;
	}

	_modelParam=_baselineModelParam;//restore the nominal values
	outf_onewaySA.close();



	return 0;
}

int Project_priortzTx::ChangeValue( string varName, double argVal, modelParamType & _modelParam )
{
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

