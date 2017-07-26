#include"trials.h"



txProfileType GetTx(const modelParamType & argParam, const patientType & argPat)
{
	if("TN_TOL_G1_SOF_PEG_RBV12" == argParam._armName)	return GetTx_TN_TOL_G1_SOF_PEG_RBV12(argParam, argPat);
	else if("TN_NOT_G1_SOF_SMV12" == argParam._armName)	return GetTx_TN_NOT_G1_SOF_SMV12(argParam, argPat);
	else if("TE_G1_SOF_SMV12" == argParam._armName)	return GetTx_TE_G1_SOF_SMV12(argParam, argPat);
	
	/************ Updated 11/28/2016 for Acute HCV treatment CEA analysis ****************/
	else if ("Acute_G1" == argParam._armName)	return GetNoTx();
	//else if ("Acute_G2" == argParam._armName)	return GetNoTx();
	//else if ("Acute_G3" == argParam._armName)	return GetNoTx();
	//else if ("Acute_G4" == argParam._armName)	return GetNoTx();
	else if ("F0_G1" == argParam._armName)	{
		txProfileType theTx =  GetTx_G1_SOF_LDV8_NoDiscontNoAE(argParam, argPat);
		if(ONEWAYSA_OPTION_PROJ_ACUTE == 1){
			theTx._pr_SVRgivenETR = argParam._transData._pr_svr_for_SA_F0 / theTx._pr_TxComplete;		
		}
		return theTx;
	} 
	//else if ("F0_G2" == argParam._armName)	return  GetTx_TN_TOL_G2_SOF_RBV12_16(argParam, argPat);
	//else if ("F0_G3" == argParam._armName) return GetTx_G3_SOF_DCV12_24(argParam, argPat);
	//else if ("F0_G4" == argParam._armName) return GetTx_G4_SOF_LDV12(argParam, argPat);
	else if ("F0_G1_ChangableDuration" == argParam._armName) {
		txProfileType theTx = GetTx_G1_SOF_LDV12_NoDiscontNoAE_ChangableDuration(argParam, argPat);
		if (ONEWAYSA_OPTION_PROJ_ACUTE == 1) {
			theTx._pr_SVRgivenETR = argParam._transData._pr_svr_for_SA_F0 / theTx._pr_TxComplete;
		}
		return theTx;
	}


	/************ Updated 7/13/2016 for India CEA analysis ****************/
	else if ("G1_SOF_LDV12" == argParam._armName) return GetTx_G1_SOF_LDV12(argParam,argPat);
	else if ("G1_SOF_DCV12_24" == argParam._armName) return GetTx_G1_SOF_DCV12_24(argParam,argPat);
	else if ("G3_SOF_DCV12_24" == argParam._armName) return GetTx_G3_SOF_DCV12_24(argParam,argPat);
	else if ("G4_SOF_LDV12" == argParam._armName) return GetTx_G4_SOF_LDV12(argParam,argPat);

	



	/************ Updated 12/12/2015 ****************/
	else if ("TN_TOL_G1b_SOF_LDV_ASV3" == argParam._armName) return GetTx_TN_TOL_G1b_SOF_LDV_ASV3(argParam,argPat);
	else if ("TN_TOL_G1b_SOF_DCV_SMV3" == argParam._armName) return GetTx_TN_TOL_G1b_SOF_DCV_SMV3(argParam,argPat);
	else if ("TN_TOL_G1b_SOF_DCV_ASV3" == argParam._armName) return GetTx_TN_TOL_G1b_SOF_DCV_ASV3(argParam,argPat);
	
	
	/************ Updated 8/27/2015 ****************/
	// existed:  "TN_TOL_G1_SOF_LDV8" and "TN_TOL_G1_SOF_LDV12"
	else if("TN_TOL_G1a_DCV_SOF12_24"== argParam._armName) return  GetTx_TN_TOL_G1a_DCV_SOF12_24(argParam, argPat);
	else if("TN_TOL_G1a_PrOD_RBV12_24"== argParam._armName) return  GetTx_TN_TOL_G1a_PrOD_RBV12_24(argParam, argPat);

	else if("TN_TOL_G1b_DCV_SOF12_24"== argParam._armName) return  GetTx_TN_TOL_G1b_DCV_SOF12_24(argParam, argPat);
	else if("TN_TOL_G1b_PrOD12"== argParam._armName) return  GetTx_TN_TOL_G1b_PrOD12(argParam, argPat);
	else if("TN_TOL_G1b_SOF_SMV12"== argParam._armName) return  GetTx_TN_TOL_G1b_SOF_SMV12(argParam, argPat);

	else if("TN_TOL_G2_SOF_RBV12_16"== argParam._armName) return  GetTx_TN_TOL_G2_SOF_RBV12_16(argParam, argPat);


	else if("TN_TOL_G3_DCV_SOF12_24"== argParam._armName) return  GetTx_TN_TOL_G3_DCV_SOF12_24(argParam, argPat);
	else if("TN_TOL_G3_SOF_PEG_RBV12"== argParam._armName) return  GetTx_TN_TOL_G3_SOF_PEG_RBV12(argParam, argPat);

	else if("TN_TOL_G4_LDV_SOF12"== argParam._armName) return  GetTx_TN_TOL_G4_LDV_SOF12(argParam, argPat);
	else if("TN_TOL_G4_PrOD_RBV12"== argParam._armName) return  GetTx_TN_TOL_G4_PrOD_RBV12(argParam, argPat);
	else if("TN_TOL_G4_SOF_RBV24"== argParam._armName) return  GetTx_TN_TOL_G4_SOF_RBV24(argParam, argPat);





	// existed: 'TE_G1_LDP_SOF12_24"

	else if("TE_G1a_DCV_SOF12_24"== argParam._armName) return  GetTx_TE_G1a_DCV_SOF12_24(argParam, argPat);
	else if("TE_G1a_PrOD_RBV12_24"== argParam._armName) return  GetTx_TE_G1a_PrOD_RBV12_24(argParam, argPat);
	else if("TE_G1a_SMV_SOF12"== argParam._armName) return  GetTx_TE_G1a_SMV_SOF12(argParam, argPat);


	else if("TE_G1b_DCV_SOF12_24"== argParam._armName) return  GetTx_TE_G1b_DCV_SOF12_24(argParam, argPat);
	else if("TE_G1b_PrOD12"== argParam._armName) return  GetTx_TE_G1b_PrOD12(argParam, argPat);
	else if("TE_G1b_SMV_SOF12"== argParam._armName) return  GetTx_TE_G1b_SMV_SOF12(argParam, argPat);

	else if("TE_G2_LDV_RBV12_16"== argParam._armName) return  GetTx_TE_G2_LDV_RBV12_16(argParam, argPat);
	else if("TE_G2_SOF_PEG_RBV12"== argParam._armName) return  GetTx_TE_G2_SOF_PEG_RBV12(argParam, argPat);
	else if("TE_G2_SMV_SOF12"== argParam._armName) return  GetTx_TE_G2_SMV_SOF12(argParam, argPat);

	else if("TE_G3_DCV_SOF12_24"== argParam._armName) return  GetTx_TE_G3_DCV_SOF12_24(argParam, argPat);
	else if("TE_G3_SOF_PEG_RBV12"== argParam._armName) return  GetTx_TE_G3_SOF_PEG_RBV12(argParam, argPat);


	else if("TE_G4_LDV_SOF12"== argParam._armName) return  GetTx_TE_G4_LDV_SOF12(argParam, argPat);
	else if("TE_G4_PrOD_RBV12"== argParam._armName) return  GetTx_TE_G4_PrOD_RBV12(argParam, argPat);
	else if("TE_G4_SOF_RBV24"== argParam._armName) return  GetTx_TE_G4_SOF_RBV24(argParam, argPat);









	/************ new treatment ****************/
	else if("TN_G1_LDP_SOF8" == argParam._armName)	return  GetTx_TN_G1_LDP_SOF8(argParam, argPat);
	else if("TN_G1_LDP_SOF12" == argParam._armName)	return GetTx_TN_G1_LDP_SOF12(argParam, argPat);
	else if("TN_TOL_G2_SOF_RBV12" == argParam._armName)	return GetTx_TN_TOL_G2_SOF_RBV12(argParam, argPat);	
	else if("TN_TOL_G3_SOF_RBV24" == argParam._armName)	return GetTx_TN_TOL_G3_SOF_RBV24(argParam, argPat);
	else if("TN_TOL_G4_SOF_PEG_RBV12" == argParam._armName)	return GetTx_TN_TOL_G4_SOF_PEG_RBV12(argParam, argPat);

	// TN_NOT_G1 same as TN_TOL_G1, skip here
	else if("TN_NOT_G2_SOF_RBV12" == argParam._armName)	return GetTx_TN_NOT_G2_SOF_RBV12(argParam, argPat);
	else if("TN_NOT_G3_SOF_RBV24" == argParam._armName)	return GetTx_TN_NOT_G3_SOF_RBV24(argParam, argPat);
	else if("TN_NOT_G4_SOF_RBV24" == argParam._armName)	return GetTx_TN_NOT_G4_SOF_RBV24(argParam, argPat);

	else if("TE_G1_LDP_SOF12_24" == argParam._armName) return GetTx_TE_G1_LDP_SOF12_24(argParam, argPat);
	else if("TE_G2_SOF_RBV12" == argParam._armName)	return GetTx_TE_G2_SOF_RBV12(argParam, argPat);
	else if("TE_G3_SOF_RBV24" == argParam._armName)	return GetTx_TE_G3_SOF_RBV24(argParam, argPat);
	else if("TE_G4_SOF_PEG_RBV12" == argParam._armName)	return GetTx_TE_G4_SOF_PEG_RBV12(argParam, argPat);

	/************ old SOC ****************/
	else if("TN_TOL_G1_BOC" == argParam._armName)	return GetTx_TN_TOL_G1_BOC(argParam,argPat);
	else if("TN_TOL_G1_TEL" == argParam._armName)	return GetTx_TN_TOL_G1_TEL(argParam,argPat);
	else if("TN_TOL_G2_PEG_RBV24" == argParam._armName) return GetTx_TN_TOL_G2_PEG_RBV24(argParam, argPat);
	else if("TN_TOL_G3_PEG_RBV24" == argParam._armName) return GetTx_TN_TOL_G3_PEG_RBV24(argParam, argPat);
	else if("TN_TOL_G4_PEG_RBV48" == argParam._armName) return GetTx_TN_TOL_G4_PEG_RBV48(argParam, argPat);


	else if ("TE_G1_BOC" == argParam._armName)	return GetTx_TE_G1_BOC(argParam,argPat);
	else if ("TE_G1_TEL" == argParam._armName)	return GetTx_TE_G1_TEL(argParam,argPat);
	else if("TE_G2_PEG_RBV24" == argParam._armName) return GetTx_TE_G2_PEG_RBV24(argParam, argPat);
	else if("TE_G3_PEG_RBV24" == argParam._armName) return GetTx_TE_G3_PEG_RBV24(argParam, argPat);
	else if("TE_G4_PEG_RBV48" == argParam._armName) return GetTx_TE_G4_PEG_RBV48(argParam, argPat);


	else if("FBOP" == argParam._armName) return GetTx_FBOP(argParam, argPat);





	///************ for project_priortzTx ****************/
	//else if("TN_G1_F0F3_LDP_SOF8" == argParam._armName)	return GetTx_TN_G1_F0F3_LDP_SOF8(argParam, argPat);
	//else if("TN_G1_F4_LDP_SOF12" == argParam._armName)	return GetTx_TN_G1_F4_LDP_SOF12(argParam, argPat);
	//else if("TE_G1_F0F3_LDP_SOF12" == argParam._armName)	return GetTx_TE_G1_F0F3_LDP_SOF12(argParam, argPat);
	//else if("TE_G1_F4_LDP_SOF24" == argParam._armName)	return GetTx_TE_G1_F4_LDP_SOF24(argParam, argPat);
	//else if("TN_G1_LDP_SOF" == argParam._armName)	return GetTx_TN_G1_LDP_SOF(argParam, argPat);
	//else if("TE_G1_LDP_SOF" == argParam._armName)	return GetTx_TE_G1_LDP_SOF(argParam, argPat);
	///**************************************************/


	else if("NoTx" == argParam._armName) return GetNoTx();
	else{
		ExitWithMsg("[Error] GetTx(): Unknown trial arm: "+argParam._armName);
		// below lines won't be executed
		txProfileType tx;
		return tx;
	}
}




txProfileType	GetTx_TN_TOL_G1a_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.96/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=2;

	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.76/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=2;
	}


	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}


txProfileType	GetTx_TN_TOL_G1a_PrOD_RBV12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;

	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=7;
	}

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PrOD, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType	GetTx_TN_TOL_G1b_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=2;

	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=2;
	}

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType	GetTx_TN_TOL_G1b_PrOD12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=12;

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	// probability of anemia
	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PrOD, pair<int,int>(1,theTx._txDuration)));


	return theTx;
}

txProfileType	GetTx_TN_TOL_G1b_SOF_SMV12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=12;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment

	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  0.88/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SMV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType	GetTx_TN_TOL_G2_SOF_RBV12_16	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.94/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.08; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;

	}else{
		theTx._txDuration=16;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.94/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.08; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=6;
	}

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}		



txProfileType	GetTx_TN_TOL_G3_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;

	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.88/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;
	}

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}
txProfileType	GetTx_TN_TOL_G3_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.96/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  0.91/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}


	// probability of anemia
	theTx._pr_AEAnemia =  0.21; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G4_LDV_SOF12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=12;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	// probability of anemia
	theTx._pr_AEAnemia =  0.02; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G4_PrOD_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=12;

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	// probability of anemia
	theTx._pr_AEAnemia =  0.13; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PrOD, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G4_SOF_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=24;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	theTx._pr_SVRgivenETR =  0.92/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	// probability of anemia
	theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=7;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType	GetTx_TE_G1a_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=2;

		// set the treatment timeline: [drug,(first week, last week)]
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));


	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=7;

		// set the treatment timeline: [drug,(first week, last week)]
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	}


	return theTx;
}
txProfileType	GetTx_TE_G1a_PrOD_RBV12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.96/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;

	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=7;
	}

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PrOD, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}
txProfileType	GetTx_TE_G1a_SMV_SOF12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	// probability of anemia
	theTx._pr_AEAnemia =  0.02; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SMV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}


txProfileType	GetTx_TE_G1b_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat){
	return GetTx_TE_G1a_DCV_SOF12_24(argModelParam,argPat);
}

txProfileType	GetTx_TE_G1b_PrOD12	(const modelParamType & argModelParam,const patientType & argPat){
	return GetTx_TN_TOL_G1b_PrOD12(argModelParam,argPat);
}

txProfileType	GetTx_TE_G1b_SMV_SOF12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	// probability of anemia
	theTx._pr_AEAnemia =  0.02; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SMV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}



txProfileType	GetTx_TE_G2_LDV_RBV12_16	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;
		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.87/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

		// probability of anemia
		theTx._pr_AEAnemia =  0.13; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;
	}else{
		theTx._txDuration=16;
		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

		// probability of anemia
		theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=6;
	}


	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType	GetTx_TE_G2_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){

		theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

	}else{
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.13; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}
txProfileType	GetTx_TE_G2_SMV_SOF12	(const modelParamType & argModelParam,const patientType & argPat){
	return GetTx_TE_G1b_SMV_SOF12(argModelParam,argPat);
}

txProfileType	GetTx_TE_G3_DCV_SOF12_24	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;
		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;
	}else{
		theTx._txDuration=24;
		// discountinuation rate		
		theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;
	}


	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType	GetTx_TE_G3_SOF_PEG_RBV12	(const modelParamType & argModelParam,const patientType & argPat){
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){

		theTx._pr_SVRgivenETR =  0.94/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

	}else{
		theTx._pr_SVRgivenETR =  0.86/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 

	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.21; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}



txProfileType GetTx_TE_G4_LDV_SOF12( const modelParamType & argModelParam,const patientType & argPat )
{
	return GetTx_TN_TOL_G4_LDV_SOF12(argModelParam,argPat);
}

txProfileType GetTx_TE_G4_PrOD_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=12;

	// discountinuation rate		
	theTx._pr_TxComplete = 1.00; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	theTx._pr_SVRgivenETR =  1.00/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	// probability of anemia
	theTx._pr_AEAnemia =  0.13; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PrOD, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TE_G4_SOF_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=24;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  0.67/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	// probability of anemia
	theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=7;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}














// ====================================================================================================================
// ====================================================================================================================
// ====================================================================================================================

txProfileType GetTx_TN_TOL_G1_SOF_PEG_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;	
	// discountinuation rate		
	theTx._pr_TxComplete = 1-0.02; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.998317857; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.915/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.789/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia = 0.208; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G2_SOF_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;	
	// discountinuation rate		
	theTx._pr_TxComplete = 250.0/253.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999006446; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=249.0/250.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.967/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.833/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia =  20.0/253.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G2_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;



	// treatment duration: 24 weeks
	theTx._txDuration=24;	
	// discountinuation rate		
	theTx._pr_TxComplete = 217.0/243.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.995295929; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=217.0/217.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.815/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.615/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia =  28.0/243.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=8;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G3_SOF_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 24 weeks
	theTx._txDuration=24;	
	// discountinuation rate	
	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999158575; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.92/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia = 11.0/100.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=7;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G4_SOF_PEG_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;


	// ============== OLD VERSION ====================
	// treatment duration: 12 weeks
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment
	//theTx._pr_ContinueTx = 0.998317857; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.987/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.851/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// ============== END - OLD VERSION ====================

	//// treatment duration: 12 weeks
	//theTx._txDuration=12;
	//// discountinuation rate		
	//theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	//theTx._pr_ContinueTx = 1.0; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	//if(argPat.state != s_CoCirr ){
	//	theTx._pr_SVRgivenETR =  0.96/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	//}
	//else if(argPat.state == s_CoCirr ){
	//	theTx._pr_SVRgivenETR =  0.96/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	//}

	// probability of anemia
	theTx._pr_AEAnemia =  0.208; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}




txProfileType GetTx_TN_NOT_G1_SOF_SMV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;	
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999162823; // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia = 0.02; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SMV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_NOT_G2_SOF_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	// original code: POSITRON_SOF12

	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;

	// discountinuation rate		
	theTx._pr_TxComplete = 203.0/207.0; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.998375253; // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	// treatment response
	theTx._pr_ETRgivenEOT=202.0/203.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.92/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.94/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}

	// probability of anemia
	theTx._pr_AEAnemia =  27.0/207.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}




txProfileType GetTx_TN_NOT_G3_SOF_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 24 weeks
	theTx._txDuration=24;
	// discountinuation rate	

	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999158575; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.92/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia = 11.0/100.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=7;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_NOT_G4_SOF_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 24 weeks
	theTx._txDuration=24;	
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999581324; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia =  0.11; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=7;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}




txProfileType GetTx_TE_G1_SOF_SMV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;	
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999162823; // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.93/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia = 0.02; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SMV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TE_G2_SOF_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;	
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999162823; // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.962/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.60/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia =  11.0/100.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}




txProfileType GetTx_TE_G3_SOF_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 24 weeks
	theTx._txDuration=24;	
	// discountinuation rate	
	theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.999158575; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.85/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.60/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia = 11.0/100.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=7;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TE_G4_SOF_PEG_RBV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;



	// treatment duration: 12 weeks
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.90; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.991258389; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.69/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.69/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	// probability of anemia
	theTx._pr_AEAnemia =  0.208; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=4;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}




//txProfileType GetTx_TN_G1_F0F3_LDP_SOF8(const modelParamType & argModelParam,const patientType & argPat)
//{
//	txProfileType theTx;
//	theTx._oSOC = false;
//
//	// treatment duration: 12 weeks
//	theTx._txDuration=8;
//	// discountinuation rate		
//	theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
//	theTx._pr_ContinueTx = 1; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment
//
//	// treatment response
//	theTx._pr_ETRgivenEOT=1.0;
//	theTx._pr_SVRgivenETR =  0.94/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 0.94
//	// probability of anemia
//	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
//	// duration of anemia (weeks)
//	theTx._assignedEpoUse=1;
//
//	// set the treatment timeline: [drug,(first week, last week)]
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
//
//	return theTx;
//}
//
//txProfileType GetTx_TN_G1_F4_LDP_SOF12(const modelParamType & argModelParam,const patientType & argPat)
//{
//	txProfileType theTx;
//	theTx._oSOC = false;
//
//	// treatment duration: 12 weeks
//	theTx._txDuration=12;
//	// discountinuation rate		
//	theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
//	theTx._pr_ContinueTx = 1; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment
//
//	// treatment response
//	theTx._pr_ETRgivenEOT=1.0;
//	theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //0.97
//	// probability of anemia
//	theTx._pr_AEAnemia =  0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
//	// duration of anemia (weeks)
//	theTx._assignedEpoUse=2;
//
//	// set the treatment timeline: [drug,(first week, last week)]
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
//
//	return theTx;
//}
//
//
//
//txProfileType GetTx_TN_G1_LDP_SOF( const modelParamType & argModelParam,const patientType & argPat )
//{
//
//	txProfileType theTx;
//
//	if(argPat.state != s_CoCirr){
//		return GetTx_TN_G1_F0F3_LDP_SOF8(argModelParam,argPat);
//	}else if(argPat.state == s_CoCirr ){
//		return GetTx_TN_G1_F4_LDP_SOF12(argModelParam,argPat);
//	}else{
//		return theTx;
//	}
//	
//}
//txProfileType GetTx_TE_G1_F0F3_LDP_SOF12(const modelParamType & argModelParam,const patientType & argPat)
//{
//	txProfileType theTx;
//	theTx._oSOC = false;
//
//	// treatment duration: 12 weeks
//	theTx._txDuration=12;
//	// discountinuation rate		
//	theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
//	theTx._pr_ContinueTx = 1; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment
//
//	// treatment response
//	theTx._pr_ETRgivenEOT=1.0;
//	theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
//	// probability of anemia
//	theTx._pr_AEAnemia =  0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
//	// duration of anemia (weeks)
//	theTx._assignedEpoUse=2;
//
//	// set the treatment timeline: [drug,(first week, last week)]
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
//
//	return theTx;
//}
//
//txProfileType GetTx_TE_G1_F4_LDP_SOF24(const modelParamType & argModelParam,const patientType & argPat)
//{
//	txProfileType theTx;
//	theTx._oSOC = false;
//
//
//	// treatment duration: 24 weeks
//	theTx._txDuration=24;
//	// discountinuation rate		
//	theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
//	theTx._pr_ContinueTx = 1; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment
//
//	// treatment response
//	theTx._pr_ETRgivenEOT=1.0;
//	theTx._pr_SVRgivenETR = 0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
//	// probability of anemia
//	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
//	// duration of anemia (weeks)
//	theTx._assignedEpoUse=4;
//
//	// set the treatment timeline: [drug,(first week, last week)]
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
//	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
//
//
//
//
//	return theTx;
//}
//
//txProfileType GetTx_TE_G1_LDP_SOF( const modelParamType & argModelParam,const patientType & argPat )
//{
//	txProfileType theTx;
//
//	if(argPat.state != s_CoCirr){
//		return GetTx_TE_G1_F0F3_LDP_SOF12(argModelParam,argPat);
//	}else if(argPat.state == s_CoCirr ){
//		return GetTx_TE_G1_F4_LDP_SOF24(argModelParam,argPat);
//	}else{
//		return theTx;
//	}
//}
//


// ============ updated 10/26/2014 ================

txProfileType GetTx_TN_G1_LDP_SOF8( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment duration: 8 weeks
	theTx._txDuration=8;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.99874; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 0.94G
	}
	else if(argPat.state == s_CoCirr ){
		/*ExitWithMsg("[Error] GetTx_TN_G1_LDP_SOF8: LDP_SOF8 is not for cirrhotic patients.");*/
		return GetTx_TN_G1_LDP_SOF12(argModelParam,argPat);
	}


	// probability of anemia
	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=1;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_G1_LDP_SOF12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment duration: 12 weeks
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.99916; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.96/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}


	// probability of anemia
	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}
txProfileType GetTx_TE_G1_LDP_SOF12_24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// discountinuation rate		
	theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 1; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		// treatment duration: 12 weeks
		theTx._txDuration=12;

		theTx._pr_SVRgivenETR =  0.95/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=0;

		// set the treatment timeline: [drug,(first week, last week)]
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	}
	else if(argPat.state == s_CoCirr ){
		// treatment duration: 24 weeks
		theTx._txDuration=24;

		theTx._pr_SVRgivenETR =  0.99/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;

		// set the treatment timeline: [drug,(first week, last week)]
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	}


	return theTx;
}

txProfileType GetTx_TN_TOL_G1_BOC( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;


	theTx._pr_ETRgivenEOT = 1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration = 32; //average duration of treatment (78% 28 weeks and 22% 48 weeks)
		theTx._pr_TxComplete = 0.72; // = (1 - 165/2312)
		theTx._pr_ContinueTx = 0.989786761; //pow(pr_TxComplete, 1.0/32.0); //weekly probability of continuing treatment
		theTx._pr_AEAnemia = 0.49; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		theTx._pr_SVRgivenETR =  0.67/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
		theTx._assignedEpoUse = 15;

	}
	else if( argPat.state == s_CoCirr ){
		theTx._txDuration = 48;
		theTx._pr_TxComplete = 0.58; // = (1 - 165/2312)
		theTx._pr_ContinueTx = 0.988715668; //pow(pr_TxComplete, 1.0/48.0); //weekly probability of continuing treatment
		theTx._pr_AEAnemia = 0.49; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		theTx._pr_SVRgivenETR =  0.52/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
		theTx._assignedEpoUse = 21;
	}


	// set the treatment timeline
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	if(argPat.state != s_CoCirr){
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(BOC, pair<int,int>(5,30)));
	}else{
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(BOC, pair<int,int>(5,theTx._txDuration)));
	}

	return theTx;
}

txProfileType GetTx_TN_TOL_G1_TEL( const modelParamType & argModelParam,const patientType & argPat )
{

	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration=34;

	theTx._pr_TxComplete = 0.79; // = (1 - 165/2312)
	theTx._pr_ContinueTx = 0.993090968; //pow(pr_TxComplete, 1.0/34.0); //weekly probability of continuing treatment

	theTx._pr_ETRgivenEOT = 1.0;
	if( argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.754/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
	} else if( argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.62/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
	}
	theTx._pr_AEAnemia = 0.37; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=12;

	// set the treatment timeline
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(TEL, pair<int,int>(1,12)));

	return theTx;
}



txProfileType GetTx_TE_G1_BOC(const modelParamType & argModelParam,const patientType & argPat)
{
	//Average length of treatment = (62%)*(36weeks) + (38%)*(48 weeks) = 41 weeks
	txProfileType theTx;
	theTx._oSOC = true;
	theTx._pr_ETRgivenEOT = 1.0;

	if(argPat.state != s_CoCirr ){
		theTx._txDuration = 41; //average duration of treatment (78% 28 weeks and 22% 48 weeks)
		theTx._pr_TxComplete = 89.0/132.0; // = 67.4%
		theTx._pr_ContinueTx = 0.990432271; //pow(pr_TxComplete, 1.0/41.0); //weekly probability of continuing treatment
		theTx._pr_AEAnemia = 0.41; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		theTx._pr_SVRgivenETR =  0.578/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
		theTx._assignedEpoUse = 16;

	}
	else if( argPat.state  == s_CoCirr ){
		theTx._txDuration = 48;
		theTx._pr_TxComplete = 0.727; // = (1 - 165/2312)
		theTx._pr_ContinueTx = 0.993387507; //pow(pr_TxComplete, 1.0/48.0); //weekly probability of continuing treatment
		theTx._pr_AEAnemia = 0.47; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		theTx._pr_SVRgivenETR =  0.522/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
		theTx._assignedEpoUse = 19;
	}

	// set the treatment timeline
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	if(argPat.state != s_CoCirr){
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(BOC, pair<int,int>(5,36)));
	}else{
		theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(BOC, pair<int,int>(5,theTx._txDuration)));
	}
	return theTx;
}


txProfileType GetTx_TE_G1_TEL( const modelParamType & argModelParam,const patientType & argPat )
{
	//Average length of treatment = (62%)*(36weeks) + (38%)*(48 weeks) = 41 weeks
	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration = 48; //average duration of treatment

	theTx._pr_ETRgivenEOT = 1.0;

	if(argPat.state != s_CoCirr ){
		theTx._pr_TxComplete = 0.771;
		theTx._pr_ContinueTx = 0.994596591; //pow(pr_TxComplete, 1.0/41.0); //weekly probability of continuing treatment
		theTx._pr_SVRgivenETR =  0.701/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
	}
	else if( argPat.state  == s_CoCirr ){
		theTx._pr_TxComplete = 0.641;
		theTx._pr_ContinueTx = 0.990777668;  //pow(pr_TxComplete, 1.0/48.0); //weekly probability of continuing treatment
		theTx._pr_SVRgivenETR =  0.583/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_oSOC); //CODEFLAG
	}

	theTx._pr_AEAnemia = 0.30; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	theTx._assignedEpoUse = 12;

	// set the treatment timeline
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(TEL, pair<int,int>(1,12)));


	return theTx;
}










txProfileType GetTx_TN_TOL_G3_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration=24;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.93; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.996980788;  // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	// treatment response
	theTx._pr_ETRgivenEOT=1;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.70/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.49/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.16; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=8;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G4_PEG_RBV48( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration=48;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.93; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.996980788;  // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	// treatment response
	theTx._pr_ETRgivenEOT=1;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.584/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.316/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.16; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=18;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TE_G2_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration=24;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.93; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.996980788;  // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	// treatment response
	theTx._pr_ETRgivenEOT=1;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.64/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.51/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.16; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=8;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TE_G3_PEG_RBV24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration=24;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.93; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.996980788;  // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	// treatment response
	theTx._pr_ETRgivenEOT=1;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.60/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.47/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.16; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=8;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TE_G4_PEG_RBV48( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;

	theTx._txDuration=48;

	// discountinuation rate		
	theTx._pr_TxComplete = 0.60; // pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.989414227;   // pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment


	// treatment response
	theTx._pr_ETRgivenEOT=1;
	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.31/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.24/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0.16; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=14;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_FBOP( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = true;
	// treatment duration: 48 weeks
	theTx._txDuration=48;	

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		// discountinuation rate	
		theTx._pr_TxComplete = 0.39; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = 0.980574; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.31/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG

	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_TxComplete = 0.32; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = 0.976541; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.12/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); //CODEFLAG

	}
	// probability of anemia
	theTx._pr_AEAnemia = 0.04; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=18;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(PEG, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(RBV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}



// ------- updated 12/12/2015 -------

txProfileType GetTx_TN_TOL_G1b_SOF_LDV_ASV3(const modelParamType & argModelParam,const patientType & argPat)
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=3;

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment

	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  1.0/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  1.0/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(ASV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G1b_SOF_DCV_SMV3(const modelParamType & argModelParam,const patientType & argPat)
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=3;

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment

	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  1.0/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  1.0/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SMV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_TN_TOL_G1b_SOF_DCV_ASV3(const modelParamType & argModelParam,const patientType & argPat)
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	theTx._txDuration=3;

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment

	theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  1.0/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}else{
		theTx._pr_SVRgivenETR =  1.0/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}

	// probability of anemia
	theTx._pr_AEAnemia =  0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(ASV, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}


// ------------------- added 7/13/2016 for India CEA -------------------------
txProfileType GetTx_G1_SOF_LDV12( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment duration: 12 weeks
	theTx._txDuration=12;
	// discountinuation rate		
	theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 0.99916; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment
	//theTx._pr_TxComplete = 1; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	//theTx._pr_ContinueTx = 1; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;

	if(argPat.state != s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.981/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		//theTx._pr_SVRgivenETR = 0.99 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	else if(argPat.state == s_CoCirr ){
		theTx._pr_SVRgivenETR =  0.932/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}


	// probability of anemia
	theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	//theTx._pr_AEAnemia = 0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG

	// duration of anemia (weeks)
	theTx._assignedEpoUse=2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(LDP, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));



	// Re-assign SVR value for PSA analysis
	if (PSA_OPTION == 1) {
		theTx._pr_SVRgivenETR = argModelParam._transData._pr_svr_for_SA / theTx._pr_TxComplete;
	}



	return theTx;
}



txProfileType GetTx_G1_SOF_DCV12_24( const modelParamType & argModelParam,const patientType & argPat )
{
	return GetTx_TN_TOL_G1a_DCV_SOF12_24(argModelParam,argPat);
}




txProfileType GetTx_G3_SOF_DCV12_24( const modelParamType & argModelParam,const patientType & argPat )
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment response
	theTx._pr_ETRgivenEOT=1.0;
	if(argPat.state != s_CoCirr ){
		theTx._txDuration=12;

		// discountinuation rate		
		theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment

		theTx._pr_SVRgivenETR =  0.97/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;

	}else{
		theTx._txDuration=24;

		// discountinuation rate		
		theTx._pr_TxComplete = 0.98; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
		theTx._pr_ContinueTx = pow(theTx._pr_TxComplete, 1.0/(double)theTx._txDuration); //weekly probability of continuing treatment


		theTx._pr_SVRgivenETR =  0.86/theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		// probability of anemia
		theTx._pr_AEAnemia =  0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
		// duration of anemia (weeks)
		theTx._assignedEpoUse=4;
	}

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(DCV, pair<int,int>(1,theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType,pair<int,int> >(SOF, pair<int,int>(1,theTx._txDuration)));

	return theTx;
}

txProfileType GetTx_G4_SOF_LDV12( const modelParamType & argModelParam,const patientType & argPat )
{
	return GetTx_TN_TOL_G4_LDV_SOF12(argModelParam,argPat);
}



// -----------------  2016/11/22 added for acute HCV phase --------------------------
txProfileType GetTx_Acute(const modelParamType & argParam, const patientType & argPat)
{
	if ("Acute_G1" == argParam._armName) {
		txProfileType theTx = GetTx_Acute_G1_TreatmentX(argParam, argPat);
		if (ONEWAYSA_OPTION_PROJ_ACUTE == 1) {
			theTx._pr_SVRgivenETR = argParam._transData._pr_svr_for_SA_acute / theTx._pr_TxComplete;
		}
		return theTx;
	}
	else if ("Acute_G1_ChangableDuration" == argParam._armName) {
		return GetTx_Acute_G1_TreatmentX_ChangableDuration(argParam, argPat);
	}
	//else if ("Acute_G2" == argParam._armName)	return GetTx_Acute_G1_TreatmentX(argParam, argPat);
	//else if ("Acute_G3" == argParam._armName)	return GetTx_Acute_G1_TreatmentX(argParam, argPat);
	//else if ("Acute_G4" == argParam._armName)	return GetTx_Acute_G1_TreatmentX(argParam, argPat);
	else if ("F0_G1" == argParam._armName)	return GetNoTx();
	//else if ("F0_G2" == argParam._armName)	return  GetNoTx();
	//else if ("F0_G3" == argParam._armName) return GetNoTx();
	//else if ("F0_G4" == argParam._armName) return GetNoTx();
	else {
		ExitWithMsg("[Error] GetTx_Acute(): Unknown trial arm: " + argParam._armName);
		txProfileType emptyTx;
		return emptyTx;
	}

}

txProfileType GetTx_Acute_G1_TreatmentX(const modelParamType & argModelParam, const patientType & argPat)
{

	txProfileType theTx = GetTx_G1_SOF_LDV8_NoDiscontNoAE(argModelParam, argPat);
	theTx._txDuration = 6;
	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(LDP, pair<int, int>(1, theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(SOF, pair<int, int>(1, theTx._txDuration)));

	return theTx;
}



txProfileType GetTx_G1_SOF_LDV8_NoDiscontNoAE(const modelParamType & argModelParam, const patientType & argPat)
{
	txProfileType theTx;
	theTx._oSOC = false;

	// treatment duration: 12 weeks
	theTx._txDuration = 8; //base case, updated by the new guidelines, 5/4/2017

	// discountinuation rate	
	//theTx._pr_TxComplete = 0.99; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	//theTx._pr_ContinueTx = 0.99916; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 1.0; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

	// treatment response
	theTx._pr_ETRgivenEOT = 1.0;

	if (argPat.state != s_CoCirr) {
		theTx._pr_SVRgivenETR = 0.981 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	else if (argPat.state == s_CoCirr) {
		theTx._pr_SVRgivenETR = 0.932 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}
	
	// probability of anemia
	//theTx._pr_AEAnemia = 0.01; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	theTx._pr_AEAnemia = 0.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG

	// duration of anemia (weeks)
	//theTx._assignedEpoUse = 2;
	theTx._assignedEpoUse = 0.0;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(LDP, pair<int, int>(1, theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(SOF, pair<int, int>(1, theTx._txDuration)));



	// Re-assign SVR value for PSA analysis
	if (PSA_OPTION == 1) {
		theTx._pr_SVRgivenETR = argModelParam._transData._pr_svr_for_SA / theTx._pr_TxComplete;
	}



	return theTx;
}

txProfileType GetTx_G1_SOF_LDV12_NoDiscontNoAE_ChangableDuration(const modelParamType & argModelParam, const patientType & argPat)
{
	txProfileType theTx;
	theTx._oSOC = false;

	// ******************* Changable treatment duration **********************
	theTx._txDuration = argModelParam._f0_tx_duration;
	// ***********************************************************************

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 1.0; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

							  // treatment response
	theTx._pr_ETRgivenEOT = 1.0;

	if (argPat.state != s_CoCirr) {
		theTx._pr_SVRgivenETR = 0.981 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 		
	}
	else if (argPat.state == s_CoCirr) {
		theTx._pr_SVRgivenETR = 0.932 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}

	// probability of anemia
	theTx._pr_AEAnemia = 0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG

							// duration of anemia (weeks)
	theTx._assignedEpoUse = 2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(LDP, pair<int, int>(1, theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(SOF, pair<int, int>(1, theTx._txDuration)));


	// Re-assign SVR value for PSA analysis
	if (PSA_OPTION == 1) {
		theTx._pr_SVRgivenETR = argModelParam._transData._pr_svr_for_SA / theTx._pr_TxComplete;
	}



	return theTx;
}


txProfileType GetTx_Acute_G1_TreatmentX_ChangableDuration(const modelParamType & argModelParam, const patientType & argPat)
{
	txProfileType theTx;
	theTx._oSOC = false;

	// ******************* Changable treatment duration **********************
	theTx._txDuration = argModelParam._acute_tx_duration;
	// ***********************************************************************

	// discountinuation rate		
	theTx._pr_TxComplete = 1.0; //pow(pr_ContinueTx,12); //convert weekly probability to overall probability of completing treatment
	theTx._pr_ContinueTx = 1.0; //pow(theTx._pr_TxComplete, 1.0/theTx._txDuration); //weekly probability of continuing treatment

								// treatment response
	theTx._pr_ETRgivenEOT = 1.0;

	if (argPat.state != s_CoCirr) {
		theTx._pr_SVRgivenETR = 0.981 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
		
	}
	else if (argPat.state == s_CoCirr) {
		theTx._pr_SVRgivenETR = 0.932 / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA); // 
	}

	if (ONEWAYSA_OPTION_PROJ_ACUTE == 1) {
		theTx._pr_SVRgivenETR = argModelParam._transData._pr_svr_for_SA_acute / theTx._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA);
	}

	// probability of anemia
	theTx._pr_AEAnemia = 0.0; //pat.transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
							  // duration of anemia (weeks)
	theTx._assignedEpoUse = 2;

	// set the treatment timeline: [drug,(first week, last week)]
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(LDP, pair<int, int>(1, theTx._txDuration)));
	theTx._tx_timeline.insert(pair<drugType, pair<int, int> >(SOF, pair<int, int>(1, theTx._txDuration)));

	return theTx;
}


txProfileType GetTx_Retreatment(const modelParamType & argParam, const patientType & argPat) {
	return GetReTx_TreatmentNameX(argParam, argPat);
}


txProfileType GetReTx_TreatmentNameX(const modelParamType & argModelParam, const patientType & argPat) {
	// ****************** USING SOF/LDV FOR NOW ***********************
	return GetTx_G1_SOF_LDV12_NoDiscontNoAE_ChangableDuration(argModelParam, argPat);
}


// ---------------- END for acute HCV phase -------------------------------------------