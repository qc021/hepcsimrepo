#include "data_type.h"


// Convert annual transition probability according to cycle length probabilities

double funcConvertProb(double transProb)
{
	transProb = 1.0 - pow((1.0 - transProb), (CYCLE_LENGTH / CYCLES_PER_YEAR));  // p2 = 1-(1-p1)^(1/52) for weekly cycle 52.1775
	return transProb;
}

// Convert annual rates to cycle length probabilities

double funcConvertRate(double transRate)
{
	double transProb; //convert annual rate to cycle-lenth probability
	transProb = 1.0 - exp(-transRate*CYCLE_LENGTH / CYCLES_PER_YEAR);  // prob = 1 - exp(-rate*t)
	return transProb;
}


double funcConvertProbForCycle(double transProb, double numCycle)
{
	double v= 1.0 - pow((1.0 - transProb), (1.0/ numCycle));
	return v;
}

// Decide progression rate based on age and gender

double transitionType::funcFibrosisProgression(double FibAge, char FibGender)
{
	double fibProgression;

	if (FibGender == 'M') {
		if (FibAge < 50) {
			fibProgression = pr_FibProg_M50under;
		}
		else if (FibAge >= 50 && FibAge < 60) {
			fibProgression = pr_FibProg_M5059;
		}
		else if (FibAge >= 60 && FibAge < 69) {
			fibProgression = pr_FibProg_M6069;
		}
		else {
			fibProgression = pr_FibProg_M70plus;
		}
	}
	else if (FibGender == 'F') {
		if (FibAge < 50) {
			fibProgression = pr_FibProg_F50under;
		}
		else if (FibAge >= 50 && FibAge < 60) {
			fibProgression = pr_FibProg_F5059;
		}
		else if (FibAge >= 60 && FibAge < 69) {
			fibProgression = pr_FibProg_F6069;
		}
		else if (FibAge >= 70 && FibAge < 79) {
			fibProgression = pr_FibProg_F7079;
		}
		else {
			fibProgression = pr_FibProg_F80plus;
		}
	} //end of else-if loop
	return fibProgression;
} //end of function


double transitionType::funcFibrosisF0F1(double FibAge, char FibGender, int FibGenotype)
{
	double fibProgression;
	int maleProp, genotypeProp;

	if (FibGender == 'M')
		maleProp = 1;
	else
		maleProp = 0;

	if (FibGenotype == 1 || FibGenotype == 11 || FibGenotype == 12) // updated 12/12/2015
		genotypeProp = 1;
	else
		genotypeProp = 0;

	// modified by Qiushi @ 3/31/2017
	double dur = FibAge - _ageHCVacquisition;
	if (_flagUseUpdatedMetaRegressionForFibProgr) {
		fibProgression = exp(-2.0124 - 0.07589 * dur + 0.3247 * 0.9009 + 0.5063*maleProp + 0.4839*genotypeProp);
	}
	else {
		fibProgression = exp(-2.0124 - 0.07589*4.2735 + 0.3247*0.9009 + 0.5063*maleProp + 0.4839*genotypeProp);
	}
	return fibProgression;
}

double transitionType::funcFibrosisF1F2(double FibAge, char FibGender, int FibGenotype)
{
	double fibProgression;
	int maleProp, genotypeProp;

	if (FibGender == 'M')
		maleProp = 1;
	else
		maleProp = 0;

	if (FibGenotype == 1 || FibGenotype == 11 || FibGenotype == 12) // updated 12/12/2015
		genotypeProp = 1;
	else
		genotypeProp = 0;

	// modified by Qiushi @ 3/31/2017
	double dur = FibAge - _ageHCVacquisition;
	if (_flagUseUpdatedMetaRegressionForFibProgr) {
		fibProgression = exp(-1.5387 - 0.06146 * dur + 0.8001*0.19);
	}
	else {
		fibProgression = exp(-1.5387 - 0.06146*14.4294 + 0.8001*0.19);
	}
	return fibProgression;
}


double transitionType::funcFibrosisF2F3(double FibAge, char FibGender, int FibGenotype)
{
	double fibProgression;
	double ageHCV = max(0.0, FibAge - 24.4784); //current age minus the mean duration of HCV infection (which depends on fibrosis score)
	double durationHCV = min(FibAge - 0.0, 24.4784); // = 0 if current age is less than the mean duration of infection

	// modified by Qiushi @ 3/31/2017	
	if (_flagUseUpdatedMetaRegressionForFibProgr) {
		double dur = FibAge - _ageHCVacquisition;
		fibProgression = exp(-1.6038 + 0.0172*_ageHCVacquisition - 0.05939*dur + 0.4539*0.19);
	}
	else {
		fibProgression = exp(-1.6038 + 0.0172*ageHCV - 0.05939*durationHCV + 0.4539*0.23);
	}
	return fibProgression;
}

double transitionType::funcFibrosisF3F4(double FibAge, char FibGender, int FibGenotype)
{
	double fibProgression;
	int genotypeProp;

	double bloodProp = 0.31;


	if (FibGenotype == 1 || FibGenotype == 11 || FibGenotype == 12) // updated 12/12/2015
		genotypeProp = 1;
	else
		genotypeProp = 0;


	// modified by Qiushi @ 3/31/2017	
	if (_flagUseUpdatedMetaRegressionForFibProgr) {
		double dur = FibAge - _ageHCVacquisition;
		fibProgression = exp(-2.2898 + 0.01689* _ageHCVacquisition - 0.03694*dur + 0.5963* _iduProp + 1.1682*bloodProp - 0.4652*genotypeProp);
	}
	else {
		double iduProp = 0.41;
		double ageHCV = max(0.0, FibAge - 32.9554); //current age minus the mean duration of HCV infection (which depends on fibrosis score)
		double durationHCV = min(FibAge - 0.0, 32.9554); // = 0 if current age is less than the mean duration of infection

		fibProgression = exp(-2.2898 + 0.01689*ageHCV - 0.03694*durationHCV + 0.5963*iduProp + 1.1682*bloodProp - 0.4652*genotypeProp);
	}
	return fibProgression;
}



patientType::patientType()
{
	// assign default values

	_cohort = cohort_defaultSingleCohort;
	_bcScr_group = bc_group_general;
	_tx_group = tx_group_general;

	state = 0; //default state of patient is F0
	currentAge = 40; //default age of patient (if nothing is specified in the main function)
	gender = 'M'; //default gender of patient (if nothing is specified in the main function)
	flagDeCirr = 0;
	flagHCC = 0;
	flagLivTr = 0;
	flagTxEx = 0;
	flagETR = 0;
	flagCured = 0;
	flagNR = 0;
	flagReL = 0;
	flagTxComplete = 0;
	flagDeathLiv = 0;
	flagAEAnemia = 0;
	numCycleLivTr = 0;
	numCycleDeCirr = 0;

	flagEverTreated = false; // added 2016/12/1

	// ---- updated for disease burden model

	_yrBirth = 1900;


	_flagReceivable = false;
	_flagTrR = false;
	_flagBirthCohortForScr = false;

	//_insurStatus = insr_uninsured;
	_flagTrInsur = false;
	_flagAware = false;
	_flagCured = false;

	_flagInterferonTol = true;

	_cycleOnDrug_PEGRBV = 0;
	_cycleOnDrug_PI = 0;
	_cycleOnDrug_DAA = 0;
	_cycleSVRF4 = 0;

	_flagTxEligible = false;
	_txCat = txCat_NoTx;
	_flagDelayTx = false;

	//typeStateDisBurdnModel _state;
	//typeTxResponseStatus _curTxResps;

	_flagAwareByScreening = false;
	_flagAwareByUsualCare = false;

	_flagEverFailedNS5A = false; // added 7/7/2017 for US Salvage analysi
	_flagEverReceivedNS5A = false;

}

counterType::counterType()
{
	/********************* initialize counters *******************************/
	countDeCirr = 0;
	countHCC = 0;
	countLivTr = 0;
	countTxEx = 0;
	countNR = 0;
	countReL = 0;
	countDeathLiv = 0;

	for (int i = 0; i < MAXCYCLE; i++) {
		//aliveCount[i] = 0;
		//incidentDeathCount[i] = 0;
		//incidentDeathLiv[i] = 0;
		//txCountDouble[i] = 0;
		//txCountTriple[i] = 0;
		//txCountNone[i] = 0;

		//incidentDisTxDouble[i] = 0;
		//incidentDisTxTriple[i] = 0;
		//incidentStateCountF0[i] = 0;
		//incidentStateCountF1[i] = 0;
		//incidentStateCountF2[i] = 0;
		//incidentStateCountF3[i] = 0;
		//incidentStateCountCoCirr[i] = 0;
		//incidentStateCountDeCirr[i] = 0;
		//incidentStateCountDeCirr1yrPlus[i] = 0;
		//incidentStateCountHCC[i] = 0;
		//incidentStateCountLivTr[i] = 0;
		//incidentStateCountLivTr1yrPlus[i] = 0;
		//incidentStateCountCured[i] = 0;
		//incidentStateCountETR[i] = 0;

		//prevStateCountF0[i] = 0;
		//prevStateCountF1[i] = 0;
		//prevStateCountF2[i] = 0;
		//prevStateCountF3[i] = 0;
		//prevStateCountCoCirr[i] = 0;
		//prevStateCountDeCirr[i] = 0;
		//prevStateCountDeCirr1yrPlus[i] = 0;
		//prevStateCountHCC[i] = 0;
		//prevStateCountLivTr[i] = 0;
		//prevStateCountLivTr1yrPlus[i] = 0;
		//prevStateCountCured[i] = 0;
		//prevStateCountETR[i] = 0;

		//prevAEAnemia[i] = 0;
		//incidentAEAnemia[i] = 0;

		aliveCount.push_back(0);
		incidentDeathCount.push_back(0);
		incidentDeathLiv.push_back(0);
		txCountDouble.push_back(0);
		txCountTriple.push_back(0);
		txCountNone.push_back(0);

		incidentDisTxDouble.push_back(0);
		incidentDisTxTriple.push_back(0);
		incidentStateCountF0.push_back(0);
		incidentStateCountF1.push_back(0);
		incidentStateCountF2.push_back(0);
		incidentStateCountF3.push_back(0);
		incidentStateCountCoCirr.push_back(0);
		incidentStateCountDeCirr.push_back(0);
		incidentStateCountDeCirr1yrPlus.push_back(0);
		incidentStateCountHCC.push_back(0);
		incidentStateCountLivTr.push_back(0);
		incidentStateCountLivTr1yrPlus.push_back(0);
		incidentStateCountCured.push_back(0);
		incidentStateCountETR.push_back(0);

		prevStateCountF0.push_back(0);
		prevStateCountF1.push_back(0);
		prevStateCountF2.push_back(0);
		prevStateCountF3.push_back(0);
		prevStateCountCoCirr.push_back(0);
		prevStateCountDeCirr.push_back(0);
		prevStateCountDeCirr1yrPlus.push_back(0);
		prevStateCountHCC.push_back(0);
		prevStateCountLivTr.push_back(0);
		prevStateCountLivTr1yrPlus.push_back(0);
		prevStateCountCured.push_back(0);
		prevStateCountETR.push_back(0);

		prevAEAnemia.push_back(0);
		incidentAEAnemia.push_back(0);
	}
}

costType::costType()
{
	// convert cost data across years using CPI:
	//2014 annual	435.292
	//2016 half1	459.07
	//2016 half2	468.29
	//2016 annual	463.68


	// ----- 2014 dollar value -----------
	//HEALTH_INFLATION = 386.567/272.8 = 1.41703446; using medical care component of CPI from 2001 to 2010								
	//c_acute = 0; //
	//c_F0 = 728.0; //								
	//c_F1 = 728.0; //								
	//c_F2 = 737.0; //								
	//c_F3 = 1496.0; //								
	//c_CoCirr = 1745.0; //								
	//c_DeCirr = 19389.0; //								
	//c_DeCirr1yrPlus = 19389.0;  //								
	//c_HCC = 35655.0; //								
	//c_LivTr = 103102.0; //								
	//c_PostLivTr = 27057.0; //	
	//					   // added 11/23/2016
	//c_testing_RNA = 92;		// He2016
	//c_testing_genotype = 408; //He2016
	//c_testing_antiHCV = 0;
	//c_testing_preTx = 0;
	//c_testing_postTx = 0;

	//// ----- 2016 dollar value - CPI 2016 HALF1 -----------
	//c_acute = 0; //
	//c_F0 = 768; //								
	//c_F1 = 768; //								
	//c_F2 = 777; //								
	//c_F3 = 1578.0; //								
	//c_CoCirr = 1840.0; //								
	//c_DeCirr = 20448.0; //								
	//c_DeCirr1yrPlus = 20448.0;  //								
	//c_HCC = 37602.0; //								
	//c_LivTr = 108733.0; //								
	//c_PostLivTr = 28535.0; //	
	//
	//// added 11/23/2016
	//c_testing_RNA = 97;		// He2016
	//c_testing_genotype = 430; //He2016
	//c_testing_antiHCV = 0;
	//c_testing_preTx = 0;
	//c_testing_postTx = 0;
	//c_background = 0; // default assumption
	
	// ----- 2016 dollar value - CPI 2016 annual  -----------
	c_acute = 0; //
	c_F0 = 775; //								
	c_F1 = 775; //								
	c_F2 = 785; //								
	c_F3 = 1594.0; //								
	c_CoCirr = 1859.0; //								
	c_DeCirr = 20653.0; //								
	c_DeCirr1yrPlus = 20653.0;  //								
	c_HCC = 37980.0; //								
	c_LivTr = 109826.0; //								
	c_PostLivTr = 28822.0; //	

	// added 11/23/2016
	c_testing_RNA = 98;		
	c_testing_genotype = 435; 

	c_testing_antiHCV = 35.0; // = 33.0 in 2014$, He 2016

	c_testing_preTx = 0;
	c_testing_postTx = 0;

	// added 03/08/2016
	c_background = 0; // default assumption


	// ---------- unchanged ----------------
	c_ETR = 0;
	c_SE_Boc = 0;
	c_SVR = 0;
	c_PEG = 588.00; //weekly cost								
	c_RBV = 309.00; //weekly cost							
	c_BOC = 1100.00; //weekly cost								
	c_TEL = 4100.00;	// weekly cost								
	c_Epo = 0; //499.00;	// weekly cost of treating anemia with Epo								
	c_SOF = 7000.00;	// weekly cost								
	c_SMV = 5530.00;    // weekly cost								
	c_LDV = 875;	// weekly cost of LDP, the SOF+LDP combo = $7875/week [updated on Oct 10, 2014]

	c_DCV = 5250;	// $63,000 for 12 weeks (Jag's email) http://www.hepmag.com/articles/Daklinza_Sovaldi_daclatasvir_2501_27556.shtml
	c_PrOD = 6943;	// $83,320 for 12 weeks (Jag's email)


	c_ASV = 1200;
	// source [1]
	// http://www.apexbt.com/asunaprevir.html?gclid=Cj0KEQiAqK-zBRC2zaXc8MOiwfIBEiQAXPHrXn2Eb8rqDvOeunfXx_DesxKGBof0WPQzDdqMQjPsxo8aAuLa8P8HAQ
	//5mg	$287.00	In stock	
	//10mg	$484.00	In stock	
	//50mg	$871.00	In stock	
	// regimen: 100mg twice daily
	// weekly cost = 484.00 * 10 * 7 = $33880

	// source [2]: Cost Effectiveness of Daclatasvir/Asunaprevir Versus Peginterferon/Ribavirin and Protease Inhibitors for the Treatment of Hepatitis c Genotype 1b Naïve Patients in Chile
	// $77419 for DCV+ASV for 12 wks
	// ASV per week = $77419/12 - 5250 ~ $1200 per week


	// screening cost
	c_screening = 2874;
	// treatment cost - see _table_tx_cost initialized by reading from external table

	// map for drug price discount. 
	// <upperbound of year, discount>
	map_drug_cost_discount[2013] = 0.23;
	map_drug_cost_discount[2014] = 0.23;
	map_drug_cost_discount[2015] = 0.46;

	override_DAA_cost_value = 900;
	override_DAA_cost_since_when = 2030;

	scalingCoef = 1e6;



	// ------ for HCC screening analysis, added 6/20/2017 ---------
	hccscr_c_tx_transpl_1y = 0; 
	hccscr_c_tx_transpl_1yPlus = 0;
	hccscr_c_tx_res = 0;
	hccscr_c_tx_abl = 0;
	hccscr_c_tx_pal = 0;
	hccscr_c_test_surveillance = 0;
	hccscr_c_test_diagnostic = 0;


}

const double costType::GetDrugCost(drugType theDrug)
{
	switch (theDrug)
	{
	case RBV:
		return c_RBV;
		break;

	case PEG:
		return c_PEG;
		break;

	case BOC:
		return c_BOC;
		break;

	case SOF:
		return c_SOF;
		break;

	case TEL:
		return c_TEL;
		break;
	case SMV:
		return c_SMV;
		break;

	case LDP:
		return c_LDV;
		break;

	case DCV:
		return c_DCV;
		break;

	case PrOD:
		return c_PrOD;
		break;

	case ASV:
		return c_ASV;
		break;

	default:
		ExitWithMsg("[Error]GetDrugCost(): unknown drug type = " + basicToStr(theDrug));
		return 0;

	}
}

transitionType::transitionType()
{
	// Convert rates to probabilities using prob = 1 - exp(-rate*t)										
	pr_FibProg_M50under = 0.054; //Fibrosis progression rate for male age under 50 (reference	 Salomon JAMA 2003)								
	pr_FibProg_M5059 = 0.125; //Fibrosis progression rate for male age 50-59 (reference	 Salomon JAMA 2003)								
	pr_FibProg_M6069 = 0.221; //Fibrosis progression rate for male age 60-69 (reference	 Salomon JAMA 2003)								
	pr_FibProg_M70plus = 0.301; //Fibrosis progression rate for male age 70 plus (reference	 Salomon JAMA 2003)								
	pr_FibProg_F50under = 0.028; //Fibrosis progression rate for female age under 50 (reference	 Salomon JAMA 2003)								
	pr_FibProg_F5059 = 0.065; //Fibrosis progression rate for female age 50-59 (reference	 Salomon JAMA 2003)								
	pr_FibProg_F6069 = 0.114; //Fibrosis progression rate for female age 60-69 (reference	 Salomon JAMA 2003)								
	pr_FibProg_F7079 = 0.154; //Fibrosis progression rate for female age 70-79 (reference	 Salomon JAMA 2003)								
	pr_FibProg_F80plus = 0.210; //Fibrosis progression rate for female age 80 plus (reference	 Salomon JAMA 2003)								

	pr_F0_F1 = 0.117; // Thein; range (0.104-0.130)									
	pr_F1_F2 = 0.085; // Thein; range (0.075-0.096)									
	pr_F2_F3 = 0.120; // Thein; range (0.109-0.133)									
	pr_F3_CoCirr = 0.116; // Thein; range  (0.104-0.129)									
	pr_F3_HCC = 0.0;
	pr_F3SVR_HCC = 0;
	pr_regression = 0;
	pr_CoCirr_DeCirr = 0.039; // Pooled analysis; range from other studies	(0.030-0.083)								
	pr_CoCirr_HCC = 0.014; // Fattovich et al.									
	pr_SVR_CoCirr_DeCirr = 0.008; //Cardosa et al.									
	pr_SVR_CoCirr_HCC = 0.005; //Cardosa et al.									
	pr_DeCirr_HCC = 0.068; // Planas									
	pr_DeCirr_LivTr = 0.023; //estimated									
	pr_DeCirr_DeathLiv = 0.182; // Planas									
	pr_DeCirr1yrPlus_DeathLiv = 0.112; // Planas 									
	pr_HCC_LivTr = 0.040; // Lang									
	pr_HCC_DeathLiv = 0.427; // Fattovich									
	pr_LivTr_DeathLiv = 0.116; //Wolfe									
	pr_LivTr1yrPlus_DeathLiv = 0.044; //Wolfe								
									  




	
	pr_SVR_Delta_DAA = 0.0; // decrement in SVR rates used for sensitivity analysis									
	pr_SVR_Delta_oSOC = 0.0; // decrement in SVR rates used for sensitivity analysis									
	

	pr_SVR_SOC = 125.0 / 311.0; // probability of achieving SVR (using data from Phase III trial)									
	pr_ETR_SOC = 176.0 / 311.0;// probability of achieving ETR (using data from Phase III trial) 									
	pr_SVR_BOC = 213.0 / 311.0; // probability of achieving SVR (using data from Phase III trial) 									
	pr_ETR_BOC = 241.0 / 311.0; // probability of achieving ETR (using data from Phase III trial) 									

	pr_SVR_TxNaive_RGT28 = 143.0 / 147.0;
	pr_ETR_TxNaive_RGT28 = 147.0 / 147.0;
	pr_SVR_TxNaive_RGT48 = 68.0 / 127.0;
	pr_ETR_TxNaive_RGT48 = 88.0 / 127.0;

	pr_TxFailTW24_SOC = 92.0 / 311.0; // probability of detectable HCV-RNA at TW 24 (treatment failure) in SOC									
	pr_TxFailTW24_PRB48 = 33.0 / 311.0; // probability of detectable HCV-RNA at TW 24 (treatment failure) in PRB48									
	pr_TxFailTW24_RGT = 42.0 / 316.0;
	pr_ResponseTW24_TxNaive = 147.0 / (147.0 + 127.0); //probability of response between TW8-TW24 									

	// Treatment-experienced study probabilities from Phase III trial (P05101)									
	pr_TxFailTW12_TxExp_SOC = 49.0 / 80.0; // probability of detectable HCV-RNA at TW 12 (treatment failure) in SOC									
	pr_TxFailTW12_TxExp_PRB48 = 29.0 / 161.0; // probability of detectable HCV-RNA at TW 12 (treatment failure) in BOC									
	pr_ResponseTW8_TxExp = 74.0 / (74.0 + 72.0); // probability of undetectable HCV-RNA at TW 8 in RGT									
	pr_ResponseNoneTW8_TxExp = 72.0 / (74.0 + 72.0); // probability of detectable HCV-RNA at TW 8 in RGT 									
	pr_TxFailTW12_TxExp_RGT36 = (82.1 - 74.0) / 82.1; // 82.1 patients go in arm 2a and 74 have response at TW12// probability of detectable HCV-RNA at TW 12 (treatment failure) in RGT 36 week therapy									
	pr_TxFailTW12_TxExp_RGT48 = (79.9 - 52.0) / 79.9; // 79.9 patients go in arm 2b and 40 have response at TW12//probability of detectable HCV-RNA at TW 12 (treatment failure) in RGT 48 week therapy									

	pr_SVR_TxExp_PRB48 = 107.0 / 161.0;// (Phase III)									
	pr_ETR_TxExp_PRB48 = 124.0 / 161.0;// (Phase III)									
	pr_SVR_TxExp_RGT36 = 64.0 / 82.1;// (Phase III)									
	pr_ETR_TxExp_RGT36 = 72.0 / 82.1;// (Phase III)									
	pr_SVR_TxExp_RGT48 = 29.0 / 79.9;// (Phase III)									
	pr_ETR_TxExp_RGT48 = 40.0 / 79.9;// (Phase III)									
	pr_SVR_TxExp_SOC = 17.0 / 80.0;// (Phase III)									
	pr_ETR_TxExp_SOC = 25.0 / 80.0;// (Phase III)									

	//AE related probabilities									
	pr_AEAnemia_SOC = 87.0 / 363.0; // probability of getting anemia on PEG+RIB 									
	pr_AEAnemia_RGT = 159.0 / 368.0;
	pr_AEAnemia_BOC = 159.0 / 366.0; // probability of getting anemia on Boceprevir 									
	pr_AERash_SOC = 0.37; // probability of getting rash on PEG+RIB 									
	pr_AERash_BOC = 0.40; // probability of getting rash on Bocepravir 									
	pr_AEAnemia_TxExp_SOC = 17.0 / 80.0; // #epo=17--Table 39; #anemia = 20--Table 22, whereas in Table 29&30 #anemia=16 									
	pr_AEAnemia_TxExp_RGT = 66.0 / 162.0; // #epo=66--Table 39; #anemia = 78--Table 22, whereas in Table 29&30 #anemia=70 									
	pr_AEAnemia_TxExp_PRB48 = 74.0 / 161.; // #epo=74--Table 39; #anemia = 79--Table 22, whereas in Table 29&30 #anemia=75&74 									

	pr_ContinueTx_TxNaive_SOC = 0.993212504;// 0.996065250; // weekly probability of continuing treatment									
	pr_ContinueTx_TxNaive_RGT = 0.994482859; //0.996589633; // weekly probability of continuing treatment									
	pr_ContinueTx_TxNaive_BOC = 0.992369264; //0.996048611; // weekly probability of continuing treatment									
	pr_ContinueTx_TxExp_SOC = 0.9960;//0.999046483; // weekly probability of continuing treatment									
	pr_ContinueTx_TxExp_RGT = 0.9957;//0.994591991; // (48 week arm only) weekly probability of continuing treatment									
	pr_ContinueTx_TxExp_BOC = 0.9954;//0.996782276; // weekly probability of continuing treatment									

	pr_EpoUse_2W = 0.24;
	pr_EpoUse_8W = 0.47;
	pr_EpoUse_18W = 0.23;
	pr_EpoUse_30W = 0.06;
	pr_EpoUse_42W = 0;





	_fibDistr.push_back(0.0975);
	_fibDistr.push_back(0.2543);
	_fibDistr.push_back(0.2066);
	_fibDistr.push_back(0.1960);
	_fibDistr.push_back(0.2457);

	_genoDistr.push_back(0.796);
	_genoDistr.push_back(0.13);
	_genoDistr.push_back(0.063);
	_genoDistr.push_back(0.011);





	_p8Week_LDV = 0.57;
	_pMale = 0.64;
	_pIFN_intol = 0.23;
	_pTE = 0.39;
	_pBOC = 0.5;

	_mfDistr.push_back(_pMale);
	_mfDistr.push_back(1 - _pMale);

	_wks_acute_svr12 = 12; // 12 weeks
	_wks_acute_wait_and_see = 26; // 6 months = 26 weeks	
	


	_pr_acute_spont_clearance = 0.25;
	_pr_svr_for_SA = 0.981;			// use SOF+LDV for non-cirrhotic patients, 12 wk, TARGET trial
	_pr_svr_for_SA_acute = 0.981;	// use SOF+LDV for non-cirrhotic patients, 12 wk, TARGET trial
	_pr_svr_for_SA_F0 = 0.981;		// use SOF+LDV for non-cirrhotic patients, 12 wk, TARGET trial
	_pr_acute_prob_complete_tx = 1.0;
	_pr_chronic_prob_follow_up = 1.0;


	_flag_override_pefct_tx = false;

	// modified by Qiushi # 03/31/2017: update metaregression model for fibrosis progression
	_ageHCVacquisition = 25.5; // thein(2008), average age of HCV acquisition in US model
	_iduProp = 0.41; // thein(2008), average % of IDU
	_flagUseUpdatedMetaRegressionForFibProgr = false; // default: use the previuos implementation of meta-regression formula for US model


	// added by Qiushi @ 6/19/2017 for HCC screening analysis
	_map_inci_HCC[scr_fib_F3] = 0.008;
	_map_inci_HCC[scr_fib_F4] = 0.016;
	_map_inci_HCC[scr_fib_DC] = 0.078;
	_map_inci_HCC[scr_fib_F3_SVR] = 0.002;
	_map_inci_HCC[scr_fib_F4_SVR] = 0.005;
	_pr_monthly_hcc_small_med = 0.056;
	_pr_monthly_hcc_med_large = 0.036;

	_mort_hcc_large = 0.75;
	_mort_F4 = 0.051; 
	_mort_DC = 0.265;

	_mort_tx_transpl_1y = 0;  
	_mort_tx_transpl_1yPlus = 0;
	_mort_tx_res_small = 0;
	_mort_tx_res_med = 0;
	_mort_tx_res_large = 0;
	_mort_tx_abl_cc_small = 0;
	_mort_tx_abl_cc_med = 0;
	_mort_tx_abl_dc_small = 0;
	_mort_tx_abl_dc_med = 0;
	_mort_tx_pal = 0;


	_sens_surveillance = 0.75;
	_spec_surveillance = 0.94;
	_sens_diagnostic = 0.9;

	_pr_transpl_hcv = funcConvertProbForCycle(0.364, 3.0);
	_pr_transpl_hcc = funcConvertProbForCycle(0.442, 3.0);

}

qolType::qolType()
{
	q_acute = 0.93;
	q_F0 = 0.93;//0.73/0.93; //0.95;					
	q_F1 = 0.93;//0.73/0.93;//0.95;					
	q_F2 = 0.93;//0.73/0.93;//0.92;					
	q_F3 = 0.93;//0.73/0.93;//0.92;					
	q_CoCirr = 0.90;//0.73/0.93;//0.89;					
	q_DeCirr = 0.80;//0.69/0.93;// = 0.65*0.62+0.55*0.28+0.53*0.1; //0.81;					
	q_HCC = 0.79;//0.51/0.93;//0.81;					
	q_LivTr = 0.84;//0.70/0.93;//0.86;					
	q_PostLivTr = 0.84;//0.70/0.93;//0.86;					
	q_SVR = 1.0;//0.77/0.93;						
	q_ETR = 1.0;//0.77/0.93;//0.90					
	q_TX_oSOC = 0.90; //multiplying factor					
	q_TX_DAA = 0.95; //multiplying factor					
	q_Dec_Anemia = 0.83; //multiplying factor					

	//q_PEG_Rib_Boc =	0.90;				
	//q_PEG_Rib = 0.90; //multiplying factor					
	//q_Boc_noSE = 0.90; //multiplying factor					
	//q_Boc_SE = 0.79; //not used					



	// ------ for HCC screening analysis, added 6/20/2017 ---------
	q_tx_transpl_1y = 0;
	q_tx_transpl_1yPlus = 0;
	qChange_tx_res_1y = 0;
	qChange_tx_res_1yPlus = 0;
	qChange_tx_abl_1y = 0;
	qChange_tx_abl_1yPlus = 0;
	qChange_tx_pal_1y = 0;
	qChange_tx_pal_1yPlus = 0;
}

txProfileType::txProfileType()
{
	// assign the default values
	_assignedEpoUse = 0; // total assigned time for Epo use for anemia
	_beginEpoUse = 2; //week to begin treatment for anemia with Epo
	_adjustEpoUse = 0; //adjustment for Epo use to incorporate AE discontinuations

	_pr_TxFailTW12 = 0.0;
}

double txProfileType::GetTxCostWeek(costType & argCost, int txWeek)
{
	double c = 0;
	for (map<drugType, pair<int, int> >::const_iterator it = _tx_timeline.begin();
		it != _tx_timeline.end(); it++) {
		if (txWeek >= it->second.first && txWeek <= it->second.second) {
			c = c + argCost.GetDrugCost(it->first);
		}
	}
	return c;
}

int psaDistrType::ReadDistrForPSA(string argFileStr)
{
	cout << "Reading parameter distributions for PSA ..." << endl;
	_distrPSA_varName.clear();
	_distrPSA_distrName.clear();
	_distrPSA_param1.clear();
	_distrPSA_param2.clear();

	ifstream inf;
	inf.open(argFileStr);

	string varName, line;
	while (inf >> varName) {
		if ("//" == varName) {
			getline(inf, line);
			continue;
		}
		string distrName;
		inf >> distrName;
		double p1, p2;
		inf >> p1;
		inf >> p2;
		_distrPSA_varName.push_back(varName);
		_distrPSA_distrName.push_back(distrName);
		_distrPSA_param1.push_back(p1);
		_distrPSA_param2.push_back(p2);
	}
	assert(_distrPSA_distrName.size() == _distrPSA_param1.size() && _distrPSA_param1.size() == _distrPSA_param2.size());

	inf.close();

	return 0;
}

int modelParamType::ReduceDrugCost(double reductPerct)
{
	_costData.c_TEL = _costData.c_TEL * (1 + reductPerct);
	_costData.c_BOC = _costData.c_BOC * (1 + reductPerct);
	_costData.c_RBV = _costData.c_RBV * (1 + reductPerct);
	_costData.c_PEG = _costData.c_PEG * (1 + reductPerct);
	_costData.c_SMV = _costData.c_SMV * (1 + reductPerct);
	_costData.c_SOF = _costData.c_SOF * (1 + reductPerct);
	_costData.c_LDV = _costData.c_LDV * (1 + reductPerct);
	_costData.c_DCV = _costData.c_DCV * (1 + reductPerct);
	_costData.c_PrOD = _costData.c_PrOD * (1 + reductPerct);

	_costData.c_ASV = _costData.c_ASV * (1 + reductPerct);

	//_costData.c_TEL=_costData.c_TEL * (1.0 -0.11);
	//_costData.c_BOC=_costData.c_BOC * (1.0 -0.11);
	//_costData.c_RBV=_costData.c_RBV * (1.0 -0.11);
	//_costData.c_PEG=_costData.c_PEG * (1.0 -0.11);
	//_costData.c_SMV=_costData.c_SMV * (1.0 -0.11);
	//_costData.c_SOF=_costData.c_SOF * (1.0 -0.46);
	//_costData.c_LDV=_costData.c_LDV * (1.0 -0.46);

	//_costData.c_DCV=_costData.c_DCV * (1.0 -0.46);
	//_costData.c_PrOD=_costData.c_PrOD * (1.0 -0.46);
	return 0;


}

DALYType::DALYType() {
	_discountRate = 0.03;
	_beta = 0.04;
	_C = 0.1658;
	_K = 0;

	ReadLifeExpectancyData("./Input_life_expectancy_WHO.txt");

	_dw_f0 = 0.1;
	_dw_f1 = 0.1;
	_dw_f2 = 0.1;
	_dw_f3 = 0.1;
	_dw_CoCirr = 0.3;
	_dw_DeCirr = 0.5;
	_dw_DeCirr1yrPlus = 0.5;
	_dw_HCC = 0.6;
	_dw_LivTr = 0.2;
	_dw_LivTr1yrPlus = 0.2;
	_dw_SVR = 0;

	// for testing:
	//cout << GetYLL(0.1, 'F') << endl;
	//cout << GetYLL(2.6, 'F') << endl;
	//cout << GetYLL(18.1, 'F') << endl;
	//cout << GetYLL(37.5, 'F') << endl;
	//cout << GetYLL(67.7, 'F') << endl;
	//cout << GetYLL(82.4, 'F') << endl;
	//cout << GetYLD(2.5, s_DeCirr, 2.0) << endl;
	//cout << GetYLD(10, s_DeCirr, 3) << endl;
	//cout << GetYLD(22.5, s_DeCirr, 4) << endl;
	//cout << GetYLD(37.5, s_DeCirr, 5) << endl;
	//cout << GetYLD(52.5, s_DeCirr, 6) << endl;
	//cout << GetYLD(65.0, s_DeCirr, 15) << endl;


	//_discountRate = 0;
	//_K = 0;
	//cout << GetYLL(0.1, 'F') << endl;
	//cout << GetYLL(2.6, 'F') << endl;
	//cout << GetYLL(18.1, 'F') << endl;
	//cout << GetYLL(37.5, 'F') << endl;
	//cout << GetYLL(67.7, 'F') << endl;
	//cout << GetYLL(82.4, 'F') << endl;
	//cout << GetYLD(2.5, s_DeCirr, 2.0) << endl;
	//cout << GetYLD(10, s_DeCirr, 3) << endl;
	//cout << GetYLD(22.5, s_DeCirr, 4) << endl;
	//cout << GetYLD(37.5, s_DeCirr, 5) << endl;
	//cout << GetYLD(52.5, s_DeCirr, 6) << endl;
	//cout << GetYLD(65.0, s_DeCirr, 15) << endl;
}

double DALYType::GetYLD(double argAge, int argState, double argL) {

	double dw = GetDisutilityWeight(argState);
	if (_discountRate < EPSILON) {//no time discounting
		double partWithWeight = _C * exp(-_beta * argAge) / pow(_beta, 2.0) * (exp(-_beta * argL) * (-_beta * (argL + argAge) - 1)
			- (-_beta * argAge - 1));
		double partWithoutWeight = argL;
		return dw * (_K * partWithWeight + (1 - _K) * partWithoutWeight);

	}
	else {


		double rPlusBeta = _discountRate + _beta;
		double partWithWeight = _C * exp(_discountRate * argAge) / pow(rPlusBeta, 2.0) * (exp(-rPlusBeta * (argL + argAge)) * (-rPlusBeta * (argL + argAge) - 1)
			- exp(-rPlusBeta * argAge) * (-rPlusBeta * argAge - 1));
		double partWithoutWeight = (1 - exp(-_discountRate * argL)) / _discountRate;
		return dw * (_K * partWithWeight + (1 - _K) * partWithoutWeight);

	}
}

double DALYType::GetYLL(double argAge, char argGender) {

	double L;
	if (argGender == 'M') {
		L = Interpolate(argAge, _lifeExpct_male);
	}
	else {
		L = Interpolate(argAge, _lifeExpct_female);
	}

	if (_discountRate < EPSILON) {//no time discounting

		double partWithWeight = _C * exp(-_beta * argAge) / pow(_beta, 2.0) * (exp(-_beta * L) * (-_beta * (L + argAge) - 1)
			- (-_beta * argAge - 1));
		double partWithoutWeight = L;
		return _K * partWithWeight + (1 - _K) * partWithoutWeight;

	}
	else { // with time discounting
		double rPlusBeta = _discountRate + _beta;
		double partWithWeight = _C * exp(_discountRate * argAge) / pow(rPlusBeta, 2.0) * (exp(-rPlusBeta * (L + argAge)) * (-rPlusBeta * (L + argAge) - 1)
			- exp(-rPlusBeta * argAge) * (-rPlusBeta * argAge - 1));
		double partWithoutWeight = (1 - exp(-_discountRate * L)) / _discountRate;
		return _K * partWithWeight + (1 - _K) * partWithoutWeight;

	}

}

double DALYType::GetDisutilityWeight(int argState) {
	switch (argState) {
	case s_F0:	return _dw_f0;
	case s_F1:	return _dw_f1;
	case s_F2:	return _dw_f2;
	case s_F3:	return _dw_f3;
	case s_CoCirr:	return _dw_CoCirr;
	case s_DeCirr:	return _dw_DeCirr;
	case s_DeCirr1yrPlus:	return _dw_DeCirr1yrPlus;
	case s_HCC:	return _dw_HCC;
	case s_LivTr:	return _dw_LivTr;
	case s_LivTr1yrPlus:	return _dw_LivTr1yrPlus;
	case s_SVR:	return _dw_SVR;
	default: return ExitWithMsg("Incorrect input state = "+basicToStr(argState)+" to get DisutilityWeight @ class DALYType");
	}
}





int DALYType::ReadLifeExpectancyData(string argFile) {
	_lifeExpct_male.clear();
	_lifeExpct_female.clear();
	ifstream _infLE;
	_infLE.open(argFile);
	int val;
	while (_infLE >> val) {
		double m, f;
		_infLE >> f; _infLE >> m;
		_lifeExpct_male[val] = m;
		_lifeExpct_female[val] = f;
	}
	_infLE.close();
	return 0;
};





int modelParamType::ApplyIndiaParameters() {
	_transData._flagUseUpdatedMetaRegressionForFibProgr = true;
	_transData._ageHCVacquisition = 23.5; // Ahmad study, lowest age observed in F0 category
	_transData._iduProp = 0.41; // Ahmad study, lowest age observed in F0 category


	//HEALTH_INFLATION = 386.567/272.8 = 1.41703446; using medical care component of CPI from 2001 to 2010								


	// India Rupees
	_costData.c_testing_antiHCV = 0;
	_costData.c_testing_RNA = 0;
	_costData.c_testing_genotype = 0;
	_costData.c_testing_preTx = 8000; // used in India CEA manuscript
	_costData.c_testing_postTx = 6000; // used in India CEA manuscript

	// annual cost
	_costData.c_F0 = 2000.0; //								
	_costData.c_F1 = 2000.0; //								
	_costData.c_F2 = 2000.0; //								
	_costData.c_F3 = 2000.0; //								
	_costData.c_CoCirr = 10000;				// updated 11/23/2016	// 40000.0; //								
	_costData.c_DeCirr = 40000;				// updated 11/23/2016	// 100000.0; //								
	_costData.c_DeCirr1yrPlus = 40000;		// updated 11/23/2016	// 100000.0;  //								
	_costData.c_HCC = 60000;				// updated 11/23/2016	// 100000.0; //								
	_costData.c_LivTr = 0; //								
	_costData.c_PostLivTr = 0; //								
	_costData.c_ETR = 0;
	_costData.c_SE_Boc = 0;
	_costData.c_SVR = 0;

	// weekly cost
	_costData.c_PEG = 0;//	588.00; //weekly cost								
	_costData.c_RBV = 0;//	309.00; //weekly cost							
	_costData.c_BOC = 0;// 1100.00; //weekly cost								
	_costData.c_TEL = 0;// 4100.00;	// weekly cost								
	//_costData.c_Epo = 0; //499.00;	// weekly cost of treating anemia with Epo								
	_costData.c_SOF = 6711.0 / 4.0; // 6000.00;	// weekly cost								
	_costData.c_SMV = 0;// 5530.00;    // weekly cost								
	_costData.c_LDV = 0;// 875;	// weekly cost of LDP, the SOF+LDP combo = $7875/week [updated on Oct 10, 2014]

	_costData.c_DCV = 0;// 5250;	// $63,000 for 12 weeks (Jag's email) http://www.hepmag.com/articles/Daklinza_Sovaldi_daclatasvir_2501_27556.shtml
	_costData.c_PrOD = 0;// 6943;	// $83,320 for 12 weeks (Jag's email)


	_costData.c_ASV = 0;// 1200;
	// source [1]
	// http://www.apexbt.com/asunaprevir.html?gclid=Cj0KEQiAqK-zBRC2zaXc8MOiwfIBEiQAXPHrXn2Eb8rqDvOeunfXx_DesxKGBof0WPQzDdqMQjPsxo8aAuLa8P8HAQ
	//5mg	$287.00	In stock	
	//10mg	$484.00	In stock	
	//50mg	$871.00	In stock	
	// regimen: 100mg twice daily
	// weekly cost = 484.00 * 10 * 7 = $33880

	// source [2]: Cost Effectiveness of Daclatasvir/Asunaprevir Versus Peginterferon/Ribavirin and Protease Inhibitors for the Treatment of Hepatitis c Genotype 1b Naïve Patients in Chile
	// $77419 for DCV+ASV for 12 wks
	// ASV per week = $77419/12 - 5250 ~ $1200 per week




	_transData.pr_DeCirr_HCC = 0.068; // US value for now
	_transData.pr_DeCirr_LivTr = 0;
	//_transData.pr_DeCirr_DeathLiv = 0.182; // Planas									
	//_transData.pr_DeCirr1yrPlus_DeathLiv = 0.112; // Planas 									
	_transData.pr_HCC_LivTr = 0;
	//_transData.pr_HCC_DeathLiv = 0.427; // Fattovich									
	//_transData.pr_LivTr_DeathLiv = 0.116 ; //Wolfe									
	//_transData.pr_LivTr1yrPlus_DeathLiv = 0.044 ; //Wolfe	


	////// add population data //////////////////////////////////////////
	_pplData._distr_fib[s_F0] = 0.184;
	_pplData._distr_fib[s_F1] = 0.248;
	_pplData._distr_fib[s_F2] = 0.217;
	_pplData._distr_fib[s_F3] = 0.217;
	_pplData._distr_fib[s_CoCirr] = 0.134;
	// G1/2/3/4 : 0.312,0.025,0.618,0.045, ignore G2, renormalize
	_pplData._distr_genotype["G1"] = 0.32;
	_pplData._distr_genotype["G2"] = 0;
	_pplData._distr_genotype["G3"] = 0.634;
	_pplData._distr_genotype["G4"] = 0.046;

	_pplData._distr_gender['M'] = 0.58;
	_pplData._distr_gender['F'] = 0.42;


	////// Overwrite DALY data ////////////////////////////////////////
	_dalyData._discountRate = 0.0;
	_dalyData._beta = 0.04;
	_dalyData._C = 0.1658;
	_dalyData._K = 0;

	//_dalyData.ReadLifeExpectancyData("./Input_life_expectancy_India.txt");
	_dalyData.ReadLifeExpectancyData("./Input_life_expectancy_WHO.txt");
	// Egypt study: Estes et al 205 APT
	_dalyData._dw_f0 = 0;
	_dalyData._dw_f1 = 0;
	_dalyData._dw_f2 = 0;
	_dalyData._dw_f3 = 0;
	_dalyData._dw_CoCirr = 0;
	_dalyData._dw_DeCirr = 0.194;
	_dalyData._dw_DeCirr1yrPlus = 0.194;
	_dalyData._dw_HCC = 0.508;
	_dalyData._dw_LivTr = 0;
	_dalyData._dw_LivTr1yrPlus = 0;
	_dalyData._dw_SVR = 0;

	//FILE_background_mortality = "Input_mortality_male_female_India.in";
	return 0;
}
