#include "sim.h"


HepCSim::HepCSim()
{
	_randomSeed = 0;
	_tx_delayed_weeks = 0;
	_print_results_every_patient = false;
	// initialize
	_tempOutf.open("output_temp.txt", ios::app);
	_tempOutf << fixed << setprecision(9);

	_psa_seed.seed(PSA_SEED);

	_psa_sampler_age.Initialize(_psa_seed(), UNIFREAL);
	_psa_sampler_state.Initialize(_psa_seed(), UNIFREAL);
	_psa_sampler_gender.Initialize(_psa_seed(), UNIFREAL);
	_rnd_extra.Initialize(SIM_SEED+1, UNIFREAL);
}


int HepCSim::SampleRandomNumbers_Tx(double * rand_pr0_CRN, double * rand_prContinueTx_CRN, double & rand_prSVR, double & rand_prETR, double & rand_prAE_Anemia,
	double & rand_prEpoUse, double & rand_prTxFail, double & rand_ETRgivenTxDiscontinue, double & rand_prTxResponse) {
	for (int txWeek = 0; txWeek < MAX_TXWEEK; txWeek++)
	{
		rand_pr0_CRN[txWeek] = lcgrand(1); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); and store random number in array; it is used again in other treatment arms i.e. common random number (CRN)
		rand_prContinueTx_CRN[txWeek] = lcgrand(1); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); //store random number in array; it is used again in other treatment arms i.e. common random number (CRN)
	}

	rand_prSVR = lcgrand(1); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); 
	rand_prETR = lcgrand(1);
	rand_prAE_Anemia = lcgrand(1);
	rand_prEpoUse = lcgrand(1);
	rand_prTxFail = lcgrand(1);
	rand_ETRgivenTxDiscontinue = lcgrand(1);
	rand_prTxResponse = lcgrand(1); // this random number is only used in RGT arm (NOTE: don't use lcgrand(1) beyond this point to preserve same sequence with other arms)

	return 0;
}

int HepCSim::Run(modelParamType argModelParam)
{

	clock_t startclock = clock();
	modelParamType theParam = argModelParam;
	ResetCounters();
	ConvertParameters(theParam);	// convert to probabilities/costs for each cycle

	ifstream inFile1, inFile2;
	ofstream outFile, outFileWebModel, outFileSummary;


	//double a, b, c;
	//bool random; 
	//enum stateType {s_F0, s_F1, s_F2, s_F3, s_CoCirr, s_DeCirr, s_HCC, s_LivTr, s_LivTr1yrPlus, s_SVR, s_Death}; // this corresponds to state = {0,1,2,...}; No need to remember numbers corresponding to each state.

	double rand_pr0; //random number associated with the transition probabilities during treatment
	double rand_pr0_CRN[MAX_TXWEEK]; // arry that stores common random number
	double rand_pr; //random number associated with the transition probabilities during the natural history of disease 
	double rand_prContinueTx; //random number associated with the probability of continuing treatment in a given week 
	double rand_prContinueTx_CRN[MAX_TXWEEK]; // arry that stores common random number
	double rand_prSVR, rand_prETR; // random number associated with achieving SVR/ETR 
	double rand_prAE_Anemia; //random number to flag patients for anemia
	double rand_prEpoUse; // random number associated with the duration of Epo use
	double rand_prProfile_State; //random number to assign state for base patient profile
	double rand_prProfile_Age; //random number to assign age for base patient profile
	double rand_prProfile_Gender; //random number to assign gender for base patient profile
	int rand_patProfile; //random patient profile number (integer value)
	double rand_prTxFail; //random number associated with treatment failure (futility rule)
	double rand_ETRgivenTxDiscontinue;
	double rand_prTxResponse; //random number associated with length of RGT (36 or 48 week arm)


	vector<patientProfile> patProfile;
	int countProfile; //total number of profiles read from input file

	double eQALY, eCost, eLY;
	double c_intermed, q_intermed;  		// intermediate cost/QALY of treatment in the treatment cycle = c_INF_Rib + c_Bocep (depending on the week in treatment)
	double avgQALY, avgCost, avgLY;
	double c_drug; //total treatment cost
	double c_testing; // total testing cost
	double discountFactAdjC, discountFactAdjQ, discountFactAdjDALY;
	int  txWeek, maxCycleNum;

	int cycleNum;

	// [QC] move following declaration to txEffType class
	//double pr_ContinueTx;// = argModelParam._transData.pr_ContinueTx_TxExp_RGT; //weekly probability of continuing treatment
	//double pr_TxComplete;// = pow(pr_ContinueTx,48); //convert weekly probability to overall probability of completing treatment
	//double pr_ETRgivenEOT;//
	//double pr_SVRgivenETR;//
	////pr_TxComplete = pow(pr_ContinueTx,48); //convert weekly probability to overall probability of completing treatment
	//double pr_AEAnemia;// = argModelParam._transData.pr_AEAnemia_TxExp_RGT; //CODEFLAG
	////the following are only used in the two arms (a & b) in the if-loops
	//double pr_TxFailTW12;// = pr_TxFailTW12_RGT36;
	//double pr_ETRgivenEOT;// 
	//double pr_SVRgivenETR;//
	//double pr_ETRgivenTxDiscontinue; // probability of ETR in subjects who discontinued treatment
	//int durationEpoUse;
	//int assignedEpoUse = 0; // total assigned time for Epo use for anemia
	//int beginEpoUse = 2; //week to begin treatment for anemia with Epo
	//double adjustEpoUse = 0; //adjustment for Epo use to incorporate AE discontinuations



	double pr_Death; //probability of death from natural causes in each cycle
	double pr_mDeathAllCause[121]; //male annual mortality rate for ages 0-119;
	double pr_fDeathAllCause[121]; //female annual mortality rate for ages 0-119;
	double q_mNormal[121]; //male age-based normal health utility values. Assume HUI = 0.5 for age >100
	double q_fNormal[121]; //female age-based normal health utility values. Assume HUI = 0.5 for age >100
	double q_AgeNormal; //age-based normal health utility value within each iteration
	double q_baseline = 1;


	// define cumulative probabilities of duration of Epo use below:
	theParam._transData.pr_EpoUse_2W = 0.09;
	theParam._transData.pr_EpoUse_8W = 0.21 + 0.09;
	theParam._transData.pr_EpoUse_18W = 0.29 + 0.21 + 0.09;
	theParam._transData.pr_EpoUse_30W = 0.36 + 0.29 + 0.21 + 0.09;
	theParam._transData.pr_EpoUse_42W = 0.05 + 0.36 + 0.29 + 0.21 + 0.09;

	//*********************Open input/outpur files*******************************//

	if (PROFILE_BOOTSTRAP == 1) inFile1.open(FILE_base_cohort_profile.c_str());	//read patient profiles if running the model for multiple profiles 

	inFile2.open(FILE_background_mortality.c_str());

	if (PRINT_PATLEVEL_OUT == 1)
		outFile.open("output_PatientResults_" + theParam._armName + ".out");
	outFile << fixed << showpoint;

	if (!inFile1 || !inFile2) {
		cout << "Cannot open input file(s) in Arm: " << theParam._armName << endl
			<< "The program terminates." << endl;
		return 1;
	}
	if (!outFile || !outFileWebModel || !outFileSummary) {
		cout << "Cannot open output file(s) in Arm:" << theParam._armName << endl
			<< "The program terminates." << endl;
		return 1;
	}

	//read patient profiles into patientProfile structure
	countProfile = 0;
	while (!inFile1.eof()) {
		patientProfile pf;
		inFile1 >> pf.bAge >> pf.bState >> pf.bGender;
		patProfile.push_back(pf);
		countProfile++;
	}
	//NOTE: j is (number of profiles) + 1

	// read mortality rates (based on 



	while (!inFile2.eof()) {
		int i;
		inFile2 >> i;
		// 		if (i==0) {
		// 			string tem;
		// 			inFile2>>tem>>tem>>tem>>tem;
		// 			continue;
		// 		}// [QC] skip the first line
		// 		i=i-1;
		double dm, df, qm, qf;
		inFile2 >> dm >> df >> qm >> qf;
		pr_mDeathAllCause[i] = dm;
		pr_fDeathAllCause[i] = df;
		q_mNormal[i] = qm;
		q_fNormal[i] = qf;

		//pr_mDeathAllCause[i]=dm*1;
		//pr_fDeathAllCause[i]=dm*1;
		//q_mNormal[i]=qm;
		//q_fNormal[i]=qm;
	}

	//convert annual discount factor according to cycle length
	discountFactAdjC = pow((1 + DISCOUNT_FACT_C), (1 / CYCLES_PER_YEAR)) - 1; //weekly rate = (1 + annual rate)^(1/52) - 1
	discountFactAdjQ = pow((1 + DISCOUNT_FACT_Q), (1 / CYCLES_PER_YEAR)) - 1; //weekly rate = (1 + annual rate)^(1/52) - 1
	discountFactAdjDALY = pow((1 + DISCOUNT_FACT_DALY), (1 / CYCLES_PER_YEAR)) - 1; //weekly rate = (1 + annual rate)^(1/52) - 1

	/********************* initialize patient/average level records *******************************/
	avgQALY = 0;
	avgLY = 0;
	avgCost = 0;
	maxCycleNum = 0;
	c_drug = 0;
	c_testing = 0;
	_sim_DALY_YLD = 0;	// added @ 20160826
	_sim_DALY_YLL = 0;  // added @ 20160826

	// argModelParam._qolData.q_TX = 0.95;
	// counter data has been initialized in constructor.


	// [QC added, 10/05/2014, temporary var]
	int numPatTreated = 0;
	int numPatUncured = 0;

	// QC added 12/6/2016, temporary var
	int tempNumCleared = 0;
	int tempNumUncleared = 0;
	double temp_eQALY_start = 0;
	double temp_avgQALY_cleared = 0;
	double temp_avgQALY_uncleared = 0;

	/***************************************************************************************/
	/************** the following simulates each patient for the lifetime ******************/
	/***************************************************************************************/

	for (long patientNum = 1; patientNum <= NUM_PATIENTS; patientNum++) {

		// [QC] create a new patient object
		patientType pat;

		// #######################################################################################
		// ##### ATTENTION: MUST USE DIFFERENT SEEDS FOR DIFFERENT STREAMS!
		// #####            OTHERWISE, IT INTRODUCES CORRELATION RATHER THAN INDEPENDENCY!!!
		// #####																Qiushi, 3/2/2015
		// #######################################################################################
		// set seed of the random number based on patientNum of each run. In this way we can generate same sequence in each treatment arm; 
		lcgrandst(patientNum + _randomSeed + 1, 1); //stream 1 is reserved for treatment/follow-up component.
		lcgrandst(patientNum + _randomSeed + 2, 2); //stream 2 is reserved for "cured" component. [achieve SVR]
		lcgrandst(patientNum + _randomSeed + 3, 3); //stream 3 is reserved for natural history component.
		lcgrandst(patientNum + _randomSeed + 4, 4); //stream 4 is reserved for patient profile 
		lcgrandst(patientNum + _randomSeed + 5, 5); //stream 5 is reserved for acute phase
		

		/****************** Generate and store random number that will be used as Common Random Number in different arms of trial**********/
		SampleRandomNumbers_Tx(rand_pr0_CRN, rand_prContinueTx_CRN, rand_prSVR, rand_prETR, rand_prAE_Anemia,
			rand_prEpoUse, rand_prTxFail, rand_ETRgivenTxDiscontinue, rand_prTxResponse);

		/********************* Create patient profile *******************************/
		/* There are three ways to select patient profiles:
		1. Use a distribution observed in trial to generate age,gender, fibrosis state
		2. Bootstrap from the cohort used in the trial
		3. Use a single base profile
		*/


		if (OPTION_CREATE_PATIENT_PROFILE == 1) {
			rand_prProfile_State = lcgrand(4); //this allows to randomly select patient's base age
			rand_prProfile_Age = lcgrand(4); //this allows to randomly select patient's base age
			rand_prProfile_Gender = lcgrand(4); //this allows to randomly select patient's base age

			//distribution of fibrosis states in total population in clinical trial (ignoring subjects with missing states)
			if (rand_prProfile_State <= 0.05)
				pat.state = s_F0;
			else if (rand_prProfile_State > 0.05 && rand_prProfile_State <= 0.53 + 0.05)
				pat.state = s_F1;
			else if (rand_prProfile_State > 0.53 + 0.05 && rand_prProfile_State <= 0.21 + 0.53 + 0.05)
				pat.state = s_F2;
			else if (rand_prProfile_State > 0.21 + 0.53 + 0.05 && rand_prProfile_State <= 0.08 + 0.21 + 0.53 + 0.05)
				pat.state = s_F3;
			else
				pat.state = s_CoCirr;

			if (rand_prProfile_Gender <= 0.67)
				pat.gender = 'M';
			else
				pat.gender = 'F';

			pat.currentAge = 26 + rand_prProfile_Age*(74 - 26);
			pat.genotype = theParam._cohortData.baseGenotype;

		}
		else if (OPTION_CREATE_PATIENT_PROFILE == 2) {// Bootstrap from the cohort used in the trial
			if (PROFILE_BOOTSTRAP == 1) {
				rand_patProfile = static_cast<int>(lcgrand(4)*(countProfile - 1)); // randomly select a row number for patient profile					
				pat.state = patProfile[rand_patProfile].bState;
				pat.currentAge = patProfile[rand_patProfile].bAge;
				pat.gender = patProfile[rand_patProfile].bGender;
				pat.genotype = theParam._cohortData.baseGenotype;	// [QC] genotype in line with the baseline genotype
			}
			else {
				ExitWithMsg("[Error: HepCSim::Run()]: inconsistent config param: PROFILE_BOOTSTRP and OPTION_CREATE_PATIENT_PROFILE");
			}
		}
		else if (OPTION_CREATE_PATIENT_PROFILE == 3) {// Use a single base profile
		   // Store base profile of patient
		   // copy from patientType.baseXXXX
			pat.state = theParam._cohortData.baseState;
			pat.currentAge = theParam._cohortData.baseAge;
			pat.gender = theParam._cohortData.baseGender;
			pat.genotype = theParam._cohortData.baseGenotype;

		}
		else {
			ExitWithMsg("[Error: HepCSim::Run()]: Unknown value for OPTION_CREATE_PATIENT_PROFILE = " + basicToStr(OPTION_CREATE_PATIENT_PROFILE));
		}

		// [QC] copy the value of race/IL28B/priorTxRes from the baseline cohort
		pat.race = theParam._cohortData.baseRace;
		pat.IL28B = theParam._cohortData.baseIL28B;
		pat.priorTxRes = theParam._cohortData.basePriorTxRes;
		// [QC] save the initial state and age
		pat.initialAge = pat.currentAge;
		pat.initialState = pat.state;
		// [QC] finish defining patient's profile


		//assign baseline fibrosis qol utility weight
		if (pat.state == s_F0)
			q_baseline = theParam._qolData.q_F0;
		else if (pat.state == s_F1)
			q_baseline = theParam._qolData.q_F1;
		else if (pat.state == s_F2)
			q_baseline = theParam._qolData.q_F2;
		else if (pat.state == s_F3)
			q_baseline = theParam._qolData.q_F3;
		else if (pat.state == s_CoCirr)
			q_baseline = theParam._qolData.q_CoCirr;



		/********************* initialize patient records *******************************/
		eQALY = 0;
		eLY = 0;
		eCost = 0;
		cycleNum = 0;


		int durationEpoUse = 0;

		/********************* UPDATE treatment efficacy data *******************************/


		// [QC updated 9/30/2014]
		double v1 = 0; double v2 = 0;

		double age0 = pat.currentAge;

		//TREATMENT=false;
		//pat.currentAge=pat.initialAge+_tx_delayed_weeks/CYCLES_PER_YEAR;

		int theDelayedWeekForThisPt = _tx_delayed_weeks;

		//if(patientNum<90000){
		//	Run_NaturalHistory(pat, cycleNum, rand_pr, pr_Death, pr_mDeathAllCause, pr_fDeathAllCause,
		//		q_AgeNormal, q_mNormal,  q_fNormal,  theParam,  eQALY,  eCost,  discountFactAdjQ,  discountFactAdjC, _tx_delayed_weeks);

		//	pat.state=s_CoCirr;	
		//	pat.currentAge=pat.initialAge+_tx_delayed_weeks/CYCLES_PER_YEAR;
		//	cycleNum=_tx_delayed_weeks;
		//	
		//}else{
		//	pat.state=s_HCC;
		//	eQALY=0;
		//}


		// moved from before-treatment-loop to here
		TREATMENT = theParam._armName == "NoTx" ? false : true;

		if (OVERRIDE_NOTREATMENT) {
			TREATMENT = false;
		}

		// RNA and genotype testing before treatment
		if ((!pat.flagEverTreated) && TREATMENT) {
			eCost = eCost + (theParam._costData.c_testing_genotype + theParam._costData.c_testing_RNA) / pow((1 + discountFactAdjC), cycleNum);
			c_testing = c_testing + (theParam._costData.c_testing_genotype + theParam._costData.c_testing_RNA) / pow((1 + discountFactAdjC), cycleNum);
		}

		// ================================================================================ // 
		// Acute phase of HCV
		// Added on 11/22/2016 by Qiushi Chen
		// ----------------------------------------

		bool flag_sim_complete = false;
		if (pat.state == s_Acute) {

			bool flag_cleared = Run_AcutePhaseHCV(pat, cycleNum, theParam, eQALY, eLY, eCost,
				pr_mDeathAllCause, pr_fDeathAllCause, q_mNormal, q_fNormal, discountFactAdjQ, discountFactAdjC, c_drug);

	
			if (flag_cleared && pat.state == s_Death) {
				ExitWithMsg("Error: flag_cleared && pat.state == s_Death");
			}

			if (pat.state == s_Death) { flag_sim_complete == true; }

			// ----------  if cured/spontaneously cleared, follow background mortality only --------
			if (flag_cleared && pat.state != s_Death) {
				tempNumCleared++;
				temp_eQALY_start = eQALY;
				while (pat.currentAge < 120 && cycleNum <= static_cast<int>(TIME_HORIZON * CYCLES_PER_YEAR)) //Maximum patient age = 120; if patient is not treated and alive
				{

					_simCounter.txCountNone[cycleNum]++;
					_simCounter.prevStateCountCured[cycleNum]++;
					double pr_Death, q_AgeNormal, rand_pr0;
					if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
						pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
						q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
					}
					else {
						pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
						q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
					}

					if (EXPECTED_LIFE_ONLY == 1) { q_AgeNormal = 1; }

					pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length
					rand_pr0 = lcgrand(2); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); //generates random number between 0 and 1; RAND_MAX = 32767


					if (rand_pr0 <= pr_Death) {
						pat.state = s_Death;
						//			pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
						_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
																	//cycleNum ++; //advance the cycle clock; used for future discounting
						break; // exit the while loop
					}
					else {
						eQALY += argModelParam._qolData.q_SVR * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
						eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
						//_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
						eCost += argModelParam._costData.c_SVR * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
						eCost += argModelParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
						pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
						_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
						cycleNum++; //advance the cycle clock; used for future discounting
					}



				}   //end of while loop
				flag_sim_complete = true;


				temp_avgQALY_cleared += (eQALY - temp_eQALY_start);
			}// end of if(flag_cleared && pat.state != s_Death)


		}




		// ================================================================================ // 
		if (!flag_sim_complete) {



			// ************** CHANGE ****************************************
			// Qiushi removed the module for delaying treatment @ 2016/12/1
			// 
			if (_tx_delayed_weeks > 0) {
				ExitWithMsg("[ERROR] The module of delaying treatment has been removed on 2016/12/1 by Qiushi. Check the source code in Run() @ sim.cpp");
			}
			//Run_NaturalHistory(pat, cycleNum, rand_pr, pr_Death, pr_mDeathAllCause, pr_fDeathAllCause,
			//	q_AgeNormal, q_mNormal, q_fNormal, theParam, eQALY, eCost, discountFactAdjQ, discountFactAdjC, _tx_delayed_weeks);
			// **************************************************************

			//if(s_CoCirr==pat.state){
			//		
			//	pat.currentAge=pat.initialAge+_tx_delayed_weeks/CYCLES_PER_YEAR;
			//	cycleNum=_tx_delayed_weeks;
			//}


			//TREATMENT=false; // skip the treatment

			double temp_q_beforeTx = eQALY;
			bool treated = false;
			bool temp_flag_alive = false;

			if (pat.state != s_Death) temp_flag_alive = true;
			if (pat.state == s_CoCirr) treated = true;

			// --- moved up before AcutePhase() --------------
			//if (theParam._armName == "NoTx") {
			//	TREATMENT = false;
			//}
			//else {
			//	TREATMENT = true;
			//}

			if (TREATMENT == true && pat.state <= s_CoCirr) {		// [QC updated 9/30/2014]
				tempNumUncleared++;



				// --- Qiushi added @ 11/27/2016, testing cost before CHRONIC treatment ----- //
				eCost = eCost + theParam._costData.c_testing_preTx / pow((1 + discountFactAdjC), cycleNum);
				c_testing = c_testing + theParam._costData.c_testing_preTx / pow((1 + discountFactAdjC), cycleNum);
				// -------------------------------------------------------------------//

				// 2016/12/2: added the addition loop for "retreatment"
				int nLoopTx = ALLOW_RETREATMENT ? MAX_LINES_TREATMENT : 1;

				for (int kTx = 0; kTx < nLoopTx; kTx++) {

					txProfileType txEffData;
					txEffData = GetTx(theParam, pat);

					// UPDATE the common random number arrays
					if (kTx > 0) {
						SampleRandomNumbers_Tx(rand_pr0_CRN, rand_prContinueTx_CRN, rand_prSVR, rand_prETR, rand_prAE_Anemia,
							rand_prEpoUse, rand_prTxFail, rand_ETRgivenTxDiscontinue, rand_prTxResponse);

						txEffData = GetTx_Retreatment(theParam, pat);
					}


					// 2017/03/14: added the override with perfect treatment
					if (argModelParam._transData._flag_override_pefct_tx) {
						txEffData._pr_TxComplete = 1;
						txEffData._pr_ContinueTx = 1;
						txEffData._pr_ETRgivenEOT = 1.0;
						txEffData._pr_SVRgivenETR = 1.0 / txEffData._pr_TxComplete * (1.0 - argModelParam._transData.pr_SVR_Delta_DAA);
						txEffData._pr_AEAnemia = 0.0;
						txEffData._assignedEpoUse = 0;
					}
					//-------------------------------------------------------

					int cycleStartTx = cycleNum; // Qiushi Modified @ 12/6/2016



					/***************************************************************************************/
					/************** the following simulates the treatment duration of patients *************/
					/************** treatment based on Arm 2 of Phase III trial ****************************/
					/***************************************************************************************/
					numPatTreated++;
					treated = true;

					//flag patients for adverse events
					if (rand_prAE_Anemia < txEffData._pr_AEAnemia)
						pat.flagAEAnemia = 1;
					else
						pat.flagAEAnemia = 0;

					// [QC] delete assignEpoUse=2 here, moved it to GetTx functions.

					/***************************************************************************************/
					/************** treatment based on Phase III trial ****************************/
					/***************************************************************************************/

					// [QC] remove the below value assignment, move to GetTx functions now.
					//pr_TxFailTW12 = argModelParam._transData.pr_TxFailTW12_TxExp_RGT36;
					//pr_ETRgivenEOT = pr_ETRgivenEOT_SOF12;
					//pr_SVRgivenETR = pr_SVRgivenETR_SOF12;

					// if no. of patients completing tx (model) < no. of ETR (trial), then assign ETR to discontinued patients to match with the trial.
					if (txEffData._pr_ETRgivenEOT > 1.0)
						txEffData._pr_ETRgivenTxDiscontinue = (txEffData._pr_ETRgivenEOT - 1.0) / (1 - txEffData._pr_TxComplete)*calibrationFactor; //fraction of non-completers who achieved ETR; 
					else
						txEffData._pr_ETRgivenTxDiscontinue = 0;

					int actualTxWeeks = 1;	// [QC added, 10/5/2014]

					for (txWeek = 0; txWeek < txEffData._txDuration; txWeek++) //get X week treatment
					{



						//rand_pr0 = rand_pr0_CRN[cycleNum]; //it is used again in other treatment arms i.e. common random number (CRN)  
						//rand_prContinueTx = rand_prContinueTx_CRN[cycleNum]; 

						//int txWeekIdx=cycleNum-_tx_delayed_weeks;
						//if(txWeek != txWeekIdx){
						//	ExitWithMsg("[ERROR] Run(): unmatched txWeek and txWeekIdx:"+basicToStr(txWeek)+" vs. "+basicToStr(txWeekIdx));
						//}

						if (USE_COMMON_RANDOM_NUMBER) {
							rand_pr0 = rand_pr0_CRN[txWeek]; //it is used again in other treatment arms i.e. common random number (CRN)  
							rand_prContinueTx = rand_prContinueTx_CRN[txWeek];
						}
						else {
							rand_pr0 = lcgrand(1);
							rand_prContinueTx = lcgrand(1);
						}
						if (OVERRIDE_SYNC_RND_STRM) { rand_pr0 = lcgrand(2); }

						//if(patientNum<100) _tempOutf<<cycleNum<<"\t"<<rand_pr0<<"\t"<<rand_prContinueTx<<"\t"<<endl;

						if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
							pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
							q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
						}
						else {
							pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
							q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
						}

						pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length	

						if (EXPECTED_LIFE_ONLY == 1)
							q_AgeNormal = 1;


						//_tempOutf<<patientNum<<"\tweek:"<<txWeek<<"\trand_pr0:"<<rand_pr0<<"\trand_prContinueTx:"<<rand_prContinueTx<<"\t"<<"\tpr_Death:"<<pr_Death<<endl;


						//ASSUMPTION: During treatment, no progression in fibrosis states
						c_intermed = txEffData.GetTxCostWeek(theParam._costData, txWeek + 1);	// week starts from 1 in this function
						// [QC] modified the q_oSOC or DAA part
						double tempQ = txEffData._oSOC ? theParam._qolData.q_TX_oSOC : theParam._qolData.q_TX_DAA;
						q_intermed = tempQ * q_AgeNormal * q_baseline;
						_simCounter.txCountTriple[cycleNum]++;
						c_drug += c_intermed; //assign cost of drug therapy


						/*********Adjust cost and HRQoL factor according to adverse events ********/
						//apply anemia related costs and HRQoL

						if (txWeek >= txEffData._beginEpoUse) { //assumption: Anemia and Tx with Epo starts at week X
							if (pat.flagAEAnemia == 1 && durationEpoUse < txEffData._assignedEpoUse) {
								c_intermed = c_intermed + theParam._costData.c_Epo;
								q_intermed = q_intermed * theParam._qolData.q_Dec_Anemia;
								_simCounter.prevAEAnemia[cycleNum]++;
								durationEpoUse++;
							}
						}

						if (rand_pr0 <= pr_Death) {
							pat.state = s_Death;
							//		pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
							_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
																		//cycleNum ++; //advance the cycle clock; used for future discounting
							break;
						}
						else {
							eQALY += q_intermed * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
							eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
							_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
							eCost += c_intermed / pow((1 + discountFactAdjC), cycleNum);
							eCost += argModelParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
							pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
							_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle  

							cycleNum++; //advance the cycle clock; used for future discounting
						}

						// check if the patient will continue treatment in the next treatment week (CODEFLAG)
						//if (rand_prContinueTx > txEffData._pr_ContinueTx && txWeek < txEffData._txDuration -1){ // condition for patient to discontinue treatment (also check if tx week is not week 12 (i.e. tx < 11))


						if (txWeek < txEffData._txDuration - 1) {
							if (rand_prContinueTx > txEffData._pr_ContinueTx) {

								pat.flagNR = 1;
								pat.flagTxComplete = 0;
								pat.flagCured = 0;
								_simCounter.incidentDisTxTriple[cycleNum]++;
								pat.flagETR = 0;


								if (rand_ETRgivenTxDiscontinue < txEffData._pr_ETRgivenTxDiscontinue) { //patients who discontinue treatment but attained ETR
									pat.flagETR = 1;
									_simCounter.incidentStateCountETR[cycleNum]++;
								}

								break;
							}
							else {			//continue treatment
								actualTxWeeks++;	//[QC added, 10/5/2014]
								//txWeek = txWeek + static_cast<int>(CYCLE_LENGTH); //advance for loop according to cycle length
							}
						}
						else { //check if the patient completed her treatment // txWeek starts from 0
							pat.flagTxComplete = 1; // sane as end of treatment (EOT)
							assert(actualTxWeeks == txEffData._txDuration);	//[QC added, 10/5/2014]
							//break;
						}

					} //end of "for" loop of treatment-weeks

					//check condition: if the patient achieves ETR (CODEFLAG)


					if (pat.state != s_Death) {
						if (pat.flagTxComplete == 1 && rand_prETR < txEffData._pr_ETRgivenEOT) { //patients who complete treatment may achieve ETR (undetectable HCV RNA at the end of tx)
							pat.flagETR = 1;
							_simCounter.incidentStateCountETR[cycleNum]++;
						}
					}

					v1 = eQALY;


					/***************************************************************************************/
					/************** the following simulates the follow-up duration of patients *************/
					/************** based on XXX trial ********************************************/
					/***************************************************************************************/
					if (pat.state != s_Death && pat.flagETR == 1)
					{
						for (int kSVR = 0; kSVR < 12; kSVR++)
							//while (cycleNum < cycleStartTx+ theDelayedWeekForThisPt + txEffData._txDuration + 12) // modified @ 12/6/2016
							//while (cycleNum < theDelayedWeekForThisPt + txEffData._txDuration + 12) //
								//while (cycleNum < _tx_delayed_weeks + actualTxWeeks + 12 ) //
								//while (cycleNum < _tx_delayed_weeks + 12 ) //follow-up phase of 12 weeks in used for SOF
						{

							_simCounter.txCountNone[cycleNum]++;
							_simCounter.prevStateCountETR[cycleNum]++;

							int txWeekIdx = cycleNum - theDelayedWeekForThisPt - cycleStartTx;
							if (USE_COMMON_RANDOM_NUMBER) {
								rand_pr0 = rand_pr0_CRN[txWeekIdx];
							}
							else {
								rand_pr0 = lcgrand(1);
							}
							if (OVERRIDE_SYNC_RND_STRM) { rand_pr0 = lcgrand(2); }

							if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
								pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
								q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
							}
							else {
								pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
								q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
							}

							pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length	

							if (EXPECTED_LIFE_ONLY == 1)
								q_AgeNormal = 1;


							if (rand_pr0 <= pr_Death) {
								pat.state = s_Death;
								//		pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
								_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
								//cycleNum ++; //advance the cycle clock; used for future discounting
								break; // exit the for loop
							}
							else {
								eQALY += theParam._qolData.q_ETR * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
								eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
								_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
								eCost += theParam._costData.c_ETR * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								eCost += theParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
								_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
								cycleNum++; //advance the cycle clock; used for future discounting
							}


						}   //end of while loop
					}	// end of if loop for death condition	


					/****************************************************************************************************/
					/*************** check condition: if the patient achieves SVR ***************************************/
					/****************************************************************************************************/

					// --- Qiushi added @ 11/27/2016, testing cost AFTER treatment ----- //
					eCost = eCost + theParam._costData.c_testing_postTx / pow((1 + discountFactAdjC), cycleNum);
					eCost = eCost + theParam._costData.c_testing_RNA / pow((1 + discountFactAdjC), cycleNum);
					c_testing = c_testing + (theParam._costData.c_testing_postTx + theParam._costData.c_testing_RNA) / pow((1 + discountFactAdjC), cycleNum);

					// -------------------------------------------------------------------//

					if (pat.state != s_Death && pat.flagETR == 1 && rand_prSVR < txEffData._pr_SVRgivenETR) {

						pat.flagCured = 1;
						_simCounter.incidentStateCountCured[cycleNum]++;
					}
					else {
						pat.flagReL = 1; // QC: no use (2/26/2015)
					}


					if (pat.flagCured) { break; }	// 2016/12/2: added for retreatment
					if (pat.state == s_Death) { break; } // 2017/4/17: added for retreatment

				}// end of for(int kTx = 0; kTx < nLoopTx; kTx++) for retreatment

				////****************************************************************************************************/
				////*********following defines the "cured loop" i.e. natural progression of cured patients *************/
				////***************************************************************************************************/

				if (pat.state != s_Death && pat.flagCured == 1)
					while (pat.currentAge < 120 && cycleNum <= static_cast<int>(TIME_HORIZON * CYCLES_PER_YEAR)) //Maximum patient age = 120; if patient is not treated and alive
					{

						_simCounter.txCountNone[cycleNum]++;
						_simCounter.prevStateCountCured[cycleNum]++;

						if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
							pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
							q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
						}
						else {
							pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
							q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
						}


						pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length

						rand_pr0 = lcgrand(2); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); //generates random number between 0 and 1; RAND_MAX = 32767

						if (EXPECTED_LIFE_ONLY == 1)
							q_AgeNormal = 1;

						if (pat.state == s_CoCirr) {
							if (rand_pr0 <= pr_Death) {
								pat.state = s_Death;
								//			pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
								_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
								//cycleNum ++; //advance the cycle clock; used for future discounting
								break; // exit the while loop
							}
							else if (rand_pr0 > pr_Death && rand_pr0 <= pr_Death + theParam._transData.pr_SVR_CoCirr_DeCirr) { //move to DeCirr state
								pat.state = s_DeCirr;
								eQALY += theParam._qolData.q_DeCirr * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
								eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
								_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
								eCost += theParam._costData.c_DeCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								eCost += theParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								_simCounter.incidentStateCountDeCirr[cycleNum]++;
								_simCounter.prevStateCountDeCirr[cycleNum]++;
								pat.flagCured = 0; //change the cure status to not cured
								pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
								_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
								cycleNum++; //advance the cycle clock; used for future discounting
								break;// exit the while loop
							}
							else if (rand_pr0 > pr_Death + theParam._transData.pr_SVR_CoCirr_DeCirr && rand_pr0 <= pr_Death + theParam._transData.pr_SVR_CoCirr_DeCirr + theParam._transData.pr_SVR_CoCirr_HCC) { //move to HCC state
								pat.state = s_HCC;
								eQALY += theParam._qolData.q_HCC * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
								eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
								_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
								eCost += theParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								eCost += theParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								_simCounter.incidentStateCountHCC[cycleNum]++;
								_simCounter.prevStateCountHCC[cycleNum]++;
								pat.flagCured = 0; //change the cure status to not cured
								pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
								_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
								cycleNum++; //advance the cycle clock; used for future discounting
								break;// exit the while loop
							}
							else { //continue in SVR-F4 state (implicitly modeled)
								eQALY += theParam._qolData.q_SVR * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
								eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
								_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
								eCost += theParam._costData.c_SVR * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								eCost += theParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
								_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
								cycleNum++; //advance the cycle clock; used for future discounting
							}

						}
						else if (pat.state <= s_F3) { //natural history of cured F0-F3 patients (remain in cured state). Note s_SVR is not used. Patient's state remains as one of the baseline states

							if (rand_pr0 <= pr_Death) {
								pat.state = s_Death;
								//			pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
								_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
								//cycleNum ++; //advance the cycle clock; used for future discounting
								break; // exit the while loop
							}
							else {
								eQALY += theParam._qolData.q_SVR * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
								eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
								_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
								eCost += theParam._costData.c_SVR * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								eCost += theParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
								pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
								_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
								cycleNum++; //advance the cycle clock; used for future discounting
							}

						} //end of if else loop of (pat.state == s_CoCirr)


					}   //end of while loop

			} //end of if condition for TREATMENT


			//if(abs(pat.currentAge-age0 - eQALY)>EPSILON){
			//	cout<<pat.currentAge-age0<<"\t"<<eQALY<<endl;
			//	ExitWithMsg("[ERROR Run(): after TX: unmatched age and eQALY]");
			//}
			/***************************************************************************************/
			/************** the following simulates the natural history of liver disease ***********/
			/***************************************************************************************/
			//Run_NaturalHistory(pat, cycleNum, rand_pr, pr_Death, pr_mDeathAllCause, pr_fDeathAllCause,
			//	q_AgeNormal, q_mNormal,  q_fNormal,  theParam,  eQALY,  eCost,  discountFactAdjQ,  discountFactAdjC, 208);

			Run_NaturalHistory(pat, cycleNum, rand_pr, pr_Death, pr_mDeathAllCause, pr_fDeathAllCause,
				q_AgeNormal, q_mNormal, q_fNormal, theParam, eQALY, eLY, eCost, discountFactAdjQ, discountFactAdjC);

			//if(abs(pat.currentAge-age0 - eQALY)>EPSILON){
			//	cout<<pat.currentAge-age0<<"\t"<<eQALY<<endl;
			//	ExitWithMsg("[ERROR Run(): final: unmatched age and eQALY]");
			//}


			temp_avgQALY_uncleared += (eQALY - temp_q_beforeTx);
		}// end of if(!flag_sim_complete)



		if (PRINT_PATLEVEL_OUT == 1) {
			if (patientNum == 1) //print header text in output file
				outFile << "PatientID" << " \t " << "eQALY" << " \t " << "eCost" << " \t " << "baseState" << " \t " << "baseGender" << " \t " << "baseAge" << " \t " << "TotalAge" << " \t " << "flagTxComplete" << " \t " << "flagTxETR" << " \t"
				<< "flagCured" << " \t" << "flagReL" << " \t " << "flagAEAnemia" << " \t " << "flagDeathLiv" << " \t " << "flagDeCirr" << " \t "
				<< "flagHCC" << " \t " << "flagLivTr" << " \t " << "numCycleLivTr" << endl;

			outFile << patientNum << " \t " << eQALY << " \t " << eCost << " \t " << pat.initialState << " \t " << pat.gender << " \t " << pat.initialAge << " \t " << pat.currentAge << " \t " << pat.flagTxComplete << " \t " << pat.flagETR << " \t"
				<< pat.flagCured << " \t" << pat.flagReL << " \t" << pat.flagAEAnemia << " \t" << pat.flagDeathLiv << " \t " << pat.flagDeCirr << " \t "
				<< pat.flagHCC << " \t " << pat.flagLivTr << " \t " << pat.numCycleLivTr << endl;
		}
		//	cout << "Control Arm: Iteration\t" << patientNum << endl; //<<" \t " <<eQALY << " \t " <<eCost <<  endl;
		avgQALY = avgQALY + eQALY;
		avgLY = avgLY + eLY;
		avgCost = avgCost + eCost;

		if (maxCycleNum < cycleNum)
			maxCycleNum = cycleNum;

		if (_print_results_every_patient && patientNum % 100 == 0) _tempOutf << patientNum << "\t" << avgQALY << "\t" << avgCost << endl;



		// =============================
		// DELETED @ 20160923 by Qiushi, ONLY count YLL for deaths due to liver disease, not for all deaths.
		// Add "_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);" at each ocurrence of "_simCounter.countDeathLiv++;"
		// -----------------------------
		//// add DALY results, @ 20160826
		//if (pat.state == s_Death) {
		//	_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);
		//}
		// =============================

	} //end of patientNum for loop

	avgQALY = avgQALY / NUM_PATIENTS;
	avgLY = avgLY / NUM_PATIENTS;
	avgCost = avgCost / NUM_PATIENTS;



	if (PRINT_RAWOUTPUT == 1) {

		/*********** print "raw" output for WebModel ********************/
		outFileWebModel << "# " << currentDateTime() << " WebModelInputFileTxExpRGT.txt (Output of Tx-Experienced Arm 2 of Boceprevir Health Economic Model for WebModel)" << endl
			<< "Boceprevir's Health Economic Model for Tx-experienced population" << endl
			<< "Strategy" << "\t" << "integer" << "\t" << "2" << endl;

		outFileWebModel << "ScenarioParameters" << endl;
		outFileWebModel << "totalPopulation" << "\t" << "integer" << "\t" << NUM_PATIENTS << endl
			<< "# cycleLength in weeks" << endl << "cycleLength" << "\t" << "decimal" << "\t" << CYCLE_LENGTH << endl
			<< "discountFactorCost" << "\t" << "decimal" << "\t" << DISCOUNT_FACT_C << endl
			<< "discountFactorQALYs" << "\t" << "decimal" << "\t" << DISCOUNT_FACT_Q << endl;

		outFileWebModel << "ModelOutput" << endl;
		outFileWebModel << "Structure(Name)" << "\t" << "aliveCount" << "\t" << "incidentDeathCount" << "\t" << "incidentDeathLiv" << "\t" << "txCountDouble" << "\t"
			<< "txCountTriple" << "\t" << "txCountNone" << "\t" << "incidentDisTxDouble" << "\t" << "incidentDisTxTriple" << "\t"
			<< "incidentStateCountF0" << "\t" << "incidentStateCountF1" << "\t" << "incidentStateCountF2" << "\t"
			<< "incidentStateCountF3" << "\t" << "incidentStateCountCoCirr" << "\t" << "incidentStateCountDeCirr" << "\t" << "incidentStateCountDeCirr1yrPlus" << "\t"
			<< "incidentStateCountHCC" << "\t" << "incidentStateCountLivTr" << "\t" << "incidentStateCountLivTr1yrPlus" << "\t"
			<< "incidentStateCountETR" << "\t" << "incidentStateCountCured" << "\t" << "prevStateCountF0" << "\t"
			<< "prevStateCountF1" << "\t" << "prevStateCountF2" << "\t" << "prevStateCountF3" << "\t" << "prevStateCountCoCirr" << "\t"
			<< "prevStateCountDeCirr" << "\t" << "prevStateCountDeCirr1yrPlus" << "\t" << "prevStateCountHCC" << "\t" << "prevStateCountLivTr" << "\t" << "prevStateCountLivTr1yrPlus" << "\t"
			<< "prevStateCountETR" << "\t" << "prevStateCountCured" << "\t" << "prevAEAnemia" << endl
			<< "# Follow-up time (week)" << endl << "t" << endl;

		for (int i = 0; i < maxCycleNum; i++)
		{
			outFileWebModel << i << "\t" << _simCounter.aliveCount[i] << "\t" << _simCounter.incidentDeathCount[i] << "\t" << _simCounter.incidentDeathLiv[i] << "\t" << _simCounter.txCountDouble[i] << "\t"
				<< _simCounter.txCountTriple[i] << "\t" << _simCounter.txCountNone[i] << "\t" << _simCounter.incidentDisTxDouble[i] << "\t" << _simCounter.incidentDisTxTriple[i] << "\t"
				<< _simCounter.incidentStateCountF0[i] << "\t" << _simCounter.incidentStateCountF1[i] << "\t" << _simCounter.incidentStateCountF2[i] << "\t"
				<< _simCounter.incidentStateCountF3[i] << "\t" << _simCounter.incidentStateCountCoCirr[i] << "\t" << _simCounter.incidentStateCountDeCirr[i] << "\t" << _simCounter.incidentStateCountDeCirr1yrPlus[i] << "\t"
				<< _simCounter.incidentStateCountHCC[i] << "\t" << _simCounter.incidentStateCountLivTr[i] << "\t" << _simCounter.incidentStateCountLivTr1yrPlus[i] << "\t"
				<< _simCounter.incidentStateCountETR[i] << "\t" << _simCounter.incidentStateCountCured[i] << "\t" << _simCounter.prevStateCountF0[i] << "\t" << _simCounter.prevStateCountF1[i] << "\t"
				<< _simCounter.prevStateCountF2[i] << "\t" << _simCounter.prevStateCountF3[i] << "\t" << _simCounter.prevStateCountCoCirr[i] << "\t"
				<< _simCounter.prevStateCountDeCirr[i] << "\t" << _simCounter.prevStateCountDeCirr1yrPlus[i] << "\t" << _simCounter.prevStateCountHCC[i] << "\t" << _simCounter.prevStateCountLivTr[i] << "\t"
				<< _simCounter.prevStateCountLivTr1yrPlus[i] << "\t" << _simCounter.prevStateCountETR[i] << "\t" << _simCounter.prevStateCountCured[i] << "\t" << _simCounter.prevAEAnemia[i] << endl;
		}
	}	 //end of PRINT_RAWOUTPUT loop


	//// [QC] Below part is removed as we only evaluate one arm at one time
	//parameterType::QALY_2 = avgQALY;
	//parameterType::COST_2 = avgCost;
	//parameterType::COUNT_DECIRR_2 = _simCounter.countDeCirr;
	//parameterType::COUNT_HCC_2 = _simCounter.countHCC;
	//parameterType::COUNT_LIVTR_2 = _simCounter.countLivTr;
	//parameterType::COUNT_DEATHLIV_2 = _simCounter.countDeathLiv;

	_sim_QALY = avgQALY;
	_sim_LY = avgLY;
	_sim_cost = avgCost;
	_sim_cost_tx = c_drug / NUM_PATIENTS;
	_sim_cost_test = c_testing / NUM_PATIENTS;
	_sim_DALY_YLD = _sim_DALY_YLD / NUM_PATIENTS;
	_sim_DALY_YLL = _sim_DALY_YLL / NUM_PATIENTS;

	outFileSummary << "Average expected QALYs (" << theParam._armName << ")" << "\t" << avgQALY << endl
		<< "Average expected Cost (" << theParam._armName << ")" << "\t" << avgCost << endl
		<< "countDeCirr (" << theParam._armName << ")" << "\t" << _simCounter.countDeCirr << endl
		<< "countHCC (" << theParam._armName << ")" << "\t" << _simCounter.countHCC << endl
		<< "countLivTr (" << theParam._armName << ")" << "\t" << _simCounter.countLivTr << endl
		<< "countDeathLiv (" << theParam._armName << ")" << "\t" << _simCounter.countDeathLiv << endl
		<< "DALY-YLL (" << theParam._armName << ")" << "\t" << _sim_DALY_YLL << endl
		<< "DALY-YLD (" << theParam._armName << ")" << "\t" << _sim_DALY_YLD << endl
		<< "DALY     (" << theParam._armName << ")" << "\t" << _sim_DALY_YLL + _sim_DALY_YLD << endl;

	if (PRINT_SCREEN == 1) {
		cout << "Average expected QALYs (" << theParam._armName << ")" << "\t" << avgQALY << endl
			<< "Average expected Cost (" << theParam._armName << ")" << "\t" << avgCost << endl
			<< "countDeCirr (" << theParam._armName << ")" << "\t" << _simCounter.countDeCirr << endl
			<< "countHCC (" << theParam._armName << ")" << "\t" << _simCounter.countHCC << endl
			<< "countLivTr (" << theParam._armName << ")" << "\t" << _simCounter.countLivTr << endl
			<< "countDeathLiv (" << theParam._armName << ")" << "\t" << _simCounter.countDeathLiv << endl
			<< "countPatientsTreated (" << theParam._armName << ")" << "\t" << numPatTreated << "\t" << numPatUncured << endl
			<< "DALY-YLL (" << theParam._armName << ")" << "\t" << _sim_DALY_YLL << endl
			<< "DALY-YLD (" << theParam._armName << ")" << "\t" << _sim_DALY_YLD << endl
			<< "DALY     (" << theParam._armName << ")" << "\t" << _sim_DALY_YLL + _sim_DALY_YLD << endl
			<< "Simulation finished in " << Time(startclock) << " sec ..." << endl << endl;
	}

	inFile1.close(); inFile2.close();
	outFile.close(); outFileWebModel.close(); outFileSummary.close();

	//cout << tempNumCleared << " [" << temp_avgQALY_cleared / tempNumCleared << "] - " <<tempNumUncleared<<" ["<< temp_avgQALY_uncleared / tempNumUncleared << "];\t";

	return 0;
}

int HepCSim::Run_HealthyPopulation(modelParamType argModelParam)
{
	/********************* initialize patient records *******************************/
	double eQALY = 0;
	double eLY = 0;
	double eCost = 0;

	double pr_Death; //probability of death from natural causes in each cycle
	double pr_mDeathAllCause[121]; //male annual mortality rate for ages 0-119;
	double pr_fDeathAllCause[121]; //female annual mortality rate for ages 0-119;
	double q_mNormal[121]; //male age-based normal health utility values. Assume HUI = 0.5 for age >100
	double q_fNormal[121]; //female age-based normal health utility values. Assume HUI = 0.5 for age >100
	double q_AgeNormal; //age-based normal health utility value within each iteration
	double q_baseline = 1;
	
	ifstream inFile2;
	inFile2.open(FILE_background_mortality.c_str());
	while (!inFile2.eof()) {
		int i;
		inFile2 >> i;
		double dm, df, qm, qf;
		inFile2 >> dm >> df >> qm >> qf;
		pr_mDeathAllCause[i] = dm;
		pr_fDeathAllCause[i] = df;
		q_mNormal[i] = qm;
		q_fNormal[i] = qf;
	}
	double discountFactAdjC = pow((1 + DISCOUNT_FACT_C), (1 / CYCLES_PER_YEAR)) - 1; //weekly rate = (1 + annual rate)^(1/52) - 1
	double discountFactAdjQ = pow((1 + DISCOUNT_FACT_Q), (1 / CYCLES_PER_YEAR)) - 1; //weekly rate = (1 + annual rate)^(1/52) - 1
	
	_rnd_extra.SetSeed(SIM_SEED);
	for (long patientNum = 1; patientNum <= NUM_PATIENTS; patientNum++) {
		lcgrandst(patientNum + _randomSeed + 1, 2); //stream 2 is reserved for "cured" component. [achieve SVR]

		// [QC] create a new patient object
		patientType pat;

		//convert annual discount factor according to cycle length
		/********************* UPDATE treatment efficacy data *******************************/
		pat.currentAge = argModelParam._cohortData.baseAge;
		pat.gender = argModelParam._cohortData.baseGender;
		int cycleNum = 0;
		while (pat.currentAge < 120 && cycleNum <= static_cast<int>(TIME_HORIZON * CYCLES_PER_YEAR)) //Maximum patient age = 120; if patient is not treated and alive
		{
			double pr_Death, q_AgeNormal, rand_pr0;
			if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
				pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
				q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
			}
			else {
				pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
				q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
			}

			if (EXPECTED_LIFE_ONLY == 1) { q_AgeNormal = 1; }

			pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length

			//double rnd = _rnd_extra.GetU01();
			double rnd = lcgrand(2);
			if (rnd <= pr_Death) {
				pat.state = s_Death;
				break; // exit the while loop
			}
			else {
				eQALY += q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
				eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
				pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
				cycleNum++; //advance the cycle clock; used for future discounting
			}
		}
	}
	//r.push_back(mySim->GetAvgQALY());
	//r.push_back(mySim->GetAvgCost());
	//r.push_back(mySim->GetCounter().countDeCirr);
	//r.push_back(mySim->GetCounter().countHCC);
	//r.push_back(mySim->GetCounter().countLivTr);
	//r.push_back(mySim->GetCounter().countDeathLiv);
	//r.push_back(mySim->GetTxCost());
	//r.push_back(mySim->GetAvgLY());
	
	_sim_QALY = eQALY / NUM_PATIENTS;
	_sim_cost = 0;
	_simCounter.countDeCirr = 0;
	_simCounter.countHCC = 0;
	_simCounter.countLivTr = 0;
	_simCounter.countDeathLiv = 0;
	_sim_cost_tx = 0;
	_sim_LY = eLY / NUM_PATIENTS;
	return 0;
}


// added Qiushi Chen @11/28/2016 for Acute HCV Treatment
bool HepCSim::Run_AcutePhaseHCV(patientType & pat, int & cycleNum, modelParamType & argModelParam, double & eQALY,  double & eLY, double & eCost,
	double * pr_mDeathAllCause, double * pr_fDeathAllCause, double * q_mNormal, double * q_fNormal, double & discountFactAdjQ, double & discountFactAdjC, double & arg_c_drug) {

	bool flag_acute_cleared = false;
	double q_baseline = argModelParam._qolData.q_acute;

	// *** NOTE ***
	// arm names can only be in the following format
	// F0_XXXXX
	// Acute_XXXXX
	vector<string> ret = SplitVec(argModelParam._armName, '_');
	if (ret[0] == "F0") {
		double rnd_clearance = lcgrand(5);
		//*********************************************************************************
		//			Arm 1: Wait until F0 (after 6 month of diagnosis of acute HCV 
		//*********************************************************************************
		// run additional 6 months, wait, no cost
		while (cycleNum < argModelParam._transData._wks_acute_wait_and_see) {

			_simCounter.txCountNone[cycleNum]++;
			_simCounter.prevStateCountETR[cycleNum]++;


			double rand_pr0 = lcgrand(5);
			double pr_Death, q_AgeNormal;

			if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
				pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
				q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
			}
			else {
				pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
				q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
			}
			pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length	

			if (EXPECTED_LIFE_ONLY == 1) { q_AgeNormal = 1; }


			if (rand_pr0 <= pr_Death) {
				pat.state = s_Death;
				//		pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
				_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
															//cycleNum ++; //advance the cycle clock; used for future discounting
				//break; // exit the for loop
				return false;
			}
			else {
				eQALY += argModelParam._qolData.q_acute * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
				eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
				//_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
				eCost += argModelParam._costData.c_acute * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
				eCost += argModelParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
				pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
				_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
				cycleNum++; //advance the cycle clock; used for future discounting
			}


		}   //end of while loop

		
		if(_rnd_extra.GetU01() < argModelParam._transData._pr_chronic_prob_follow_up){// added by Qiushi @ 4/18/2017
			// test confirming clearance
			eCost = eCost + argModelParam._costData.c_testing_RNA / pow((1 + discountFactAdjC), cycleNum);

		
			if (rnd_clearance < argModelParam._transData._pr_acute_spont_clearance) {
				flag_acute_cleared = true;
				TREATMENT = false;
				// if cleared, follow the background mortality only
			}
			else {
				pat.state = s_F0;
				flag_acute_cleared = false;
				TREATMENT = true;
				if (OVERRIDE_NOTREATMENT) { TREATMENT = false; }
			}

		}
		else {
			// ----------- added by Qiushi @ 4/18/2017 --------------------------------
			// if no follow-up, patients won't get treatment and start progression.
			// ------------------------------------------------------------------------
			TREATMENT = false;
			if (rnd_clearance < argModelParam._transData._pr_acute_spont_clearance) {
				// if cleared, follow the background mortality only
				flag_acute_cleared = true;								
			}
			else {
				pat.state = s_F0;
				flag_acute_cleared = false;
				
			}
		}


	}
	else if (ret[0] == "Acute") {
		//*********************************************************************************
		//						Arm 2: Immediate Treatment in Acute Phase
		//*********************************************************************************
		// 2016/12/2: added the addition loop for "retreatment"



		txProfileType txEffData_Acute;
		txEffData_Acute = GetTx_Acute(argModelParam, pat);

		int nLoopTx = ALLOW_RETREATMENT ? MAX_LINES_TREATMENT : 1;
		bool flag_loss_follow_up = false; // added 4/20/2017: it is set to be true if the patient discountinue acute treatment

		for (int kTx = 0; kTx < nLoopTx; kTx++) {
			// UPDATE the common random number arrays
			if (kTx > 0) {
				if (flag_loss_follow_up) {
					//break;
				}
				txEffData_Acute = GetTx_Retreatment(argModelParam, pat);
			}

			
			int durationEpoUse = 0;
			int actualTxWeeks = 1;
			// ==========================================================================
			// ===== during treatment periods (typically 6 wks for acute HCV) ===========
			// ==========================================================================
			for (int txWeek = 0; txWeek < txEffData_Acute._txDuration; txWeek++) {
				double rand_pr0 = lcgrand(5);
				double rand_prContinueTx = lcgrand(5);
				
	
				double pr_Death, q_AgeNormal;
				if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
					pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
					q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];
				}
				else {
					pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
					q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
				}
				pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length	

				if (EXPECTED_LIFE_ONLY == 1) q_AgeNormal = 1;


				// *************************************************************************************
				// Run natural history with background mortality only
				// *************************************************************************************
				//ASSUMPTION: During treatment, no progression in fibrosis states
				double c_intermed = txEffData_Acute.GetTxCostWeek(argModelParam._costData, txWeek + 1);	// week starts from 1 in this function
				double tempQ = txEffData_Acute._oSOC ? argModelParam._qolData.q_TX_oSOC : argModelParam._qolData.q_TX_DAA;
				double q_intermed = tempQ * q_AgeNormal * q_baseline;

				//assign cost of drug therapy
				arg_c_drug += c_intermed;

				if (rand_pr0 <= pr_Death) {
					pat.state = s_Death;
					//		pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
					_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
																//cycleNum ++; //advance the cycle clock; used for future discounting
					return false;
				}
				else {
					eQALY += q_intermed * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
					eCost += c_intermed / pow((1 + discountFactAdjC), cycleNum);
					eCost += argModelParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
					_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle  

					cycleNum++; //advance the cycle clock; used for future discounting
				}



				// *************************************************************************************
				// Discontinue treatment? loss of follow-up?
				// check if the patient will continue treatment in the next treatment week (CODEFLAG)
				// *************************************************************************************
				if (txWeek < txEffData_Acute._txDuration - 1) {
					// ------------- modified @ 4/18/2017 ----------------------
					// if it is the treamtment for acute HCV (i.e., kTx = 0)
					// override the treatment continuation probability using param._transData.					
					double pr_continueTx = txEffData_Acute._pr_ContinueTx;
					if (1 - txEffData_Acute._pr_ContinueTx > EPSILON) { ExitWithMsg("Error @ Run_AcutePhaseHCV: we assume no discontinuation due to treatment side effects, check the paramteres in the model!"); }

					if (kTx == 0) {
						pr_continueTx = pow(argModelParam._transData._pr_acute_prob_complete_tx, 1.0 / (double)(txEffData_Acute._txDuration));
					}

					if (rand_prContinueTx > pr_continueTx) {
						//if (rand_prContinueTx > txEffData_Acute._pr_ContinueTx) {

						pat.flagNR = 1;
						pat.flagTxComplete = 0;
						pat.flagCured = 0;
						//_simCounter.incidentDisTxTriple[cycleNum]++; // removed on 4/21/2017 by Qiushi
						pat.flagETR = 0;
						//if (rand_ETRgivenTxDiscontinue < txEffData_Acute._pr_ETRgivenTxDiscontinue) { //patients who discontinue treatment but attained ETR
						//	pat.flagETR = 1;
						//	_simCounter.incidentStateCountETR[cycleNum]++;
						//}

						flag_acute_cleared = false;
						TREATMENT = false;

						flag_loss_follow_up = true; // added 4/20/2017
						break;
					}
					else {			
						//continue treatment
						//[QC added, 10/5/2014]
						actualTxWeeks++;	
						//txWeek = txWeek + static_cast<int>(CYCLE_LENGTH); //advance for loop according to cycle length
					}
				}
				else { //check if the patient completed her treatment // txWeek starts from 0
					pat.flagTxComplete = 1; // same as end of treatment (EOT)
					pat.flagETR = 1;
					_simCounter.incidentStateCountETR[cycleNum]++;
					assert(actualTxWeeks == txEffData_Acute._txDuration);	//[QC added, 10/5/2014]
																			//break;
				}
								

			} //end of "for" loop of treatment-weeks


			// ==========================================================================
			// Follow-periods after completing treatment
			// [1] extra 12 weeks for those who completed treatment: waiting for confirming SVR12
			// [2] until 6 months (26 weeks from the beginning) for those who dropped out.
			// ==========================================================================
			if (pat.state != s_Death && pat.flagETR == 1) {
				int endingCycle;
				if (flag_loss_follow_up) {
					endingCycle = argModelParam._transData._wks_acute_wait_and_see;
				}
				else {
					endingCycle = txEffData_Acute._txDuration + argModelParam._transData._wks_acute_svr12;
				}
							
				while (cycleNum < endingCycle) {

					_simCounter.txCountNone[cycleNum]++;
					_simCounter.prevStateCountETR[cycleNum]++;

					int txWeekIdx = cycleNum;
					double rand_pr0 = lcgrand(2);
										
					double pr_Death = (pat.gender == 'M' ? pr_mDeathAllCause[static_cast<int>(pat.currentAge)] : pr_fDeathAllCause[static_cast<int>(pat.currentAge)]);
					pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length	
					double q_AgeNormal = (pat.gender == 'M' ? q_mNormal[static_cast<int>(pat.currentAge)] : q_fNormal[static_cast<int>(pat.currentAge)]);
					if (EXPECTED_LIFE_ONLY == 1) { q_AgeNormal = 1; }


					if (rand_pr0 <= pr_Death) {
						pat.state = s_Death;
						//		pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR)/2 ; //half-cycle correction
						_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
																	
						return false;
						
					}
					else {
						eQALY += argModelParam._qolData.q_ETR * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
						eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
						//_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
						eCost += argModelParam._costData.c_ETR * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
						eCost += argModelParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
						pat.currentAge = pat.currentAge + (CYCLE_LENGTH / CYCLES_PER_YEAR); //advance patient's age if the patient is alive in the given cycle
						_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle
						cycleNum++; //advance the cycle clock; used for future discounting
					}

				}   //end of while loop for follow-up process after treatment
				
			}	// end of if (pat.state != s_Death && pat.flagETR == 1)
			

			// patients with death state have been excluded
			if (flag_loss_follow_up) {
				// if patients don't follow up
				// won't get further treatment
				TREATMENT = false;
				// some can get spontaneous clearance
				flag_acute_cleared = _rnd_extra.GetU01() < argModelParam._transData._pr_acute_spont_clearance;
				pat.flagCured = flag_acute_cleared;
				
				if(! flag_acute_cleared){
					pat.state = s_F0;		// and continue to natural history.			
				}

				return flag_acute_cleared;

			}else{
				//****************************************************************************************************/
				//*************** check condition: if the patient achieves SVR ***************************************/
				//****************************************************************************************************/

				// --- Qiushi added @ 11/27/2016, testing cost before treatment ----- //
				eCost = eCost + argModelParam._costData.c_testing_RNA / pow((1 + discountFactAdjC), cycleNum);
				// -------------------------------------------------------------------//
				double rand_prSVR = lcgrand(5);
				if (pat.flagETR == 1 && rand_prSVR < txEffData_Acute._pr_SVRgivenETR) {
					pat.flagCured = 1;
				}

				if (pat.flagCured) { break; }// 2016/12/2: added for retreatment
					
			}
		}// end of for(int kTx = 0; kTx < nLoopTx; kTx++) for retreatment




			if (pat.flagCured) {
				flag_acute_cleared = true;// 2016/11/28: added for Acute HCV treatment
				_simCounter.incidentStateCountCured[cycleNum]++;
			}
			else {
				flag_acute_cleared = false;// 2016/11/28: added for Acute HCV treatment			
				pat.state = s_F0;
				pat.flagReL = 1; // QC: no use (2/26/2015)
			}
			TREATMENT = false; // disable the treatment even if F0
	}
	else {
		ExitWithMsg("[Error] Undefined arm (" + argModelParam._armName + ") in the acute phase @ Run_AcutePhaseHCV() in sim.cpp");
	}




	return flag_acute_cleared;

}

int HepCSim::Run_NaturalHistory(patientType & pat, int & cycleNum, double & rand_pr, double & pr_Death, double * pr_mDeathAllCause, double * pr_fDeathAllCause,
	double & q_AgeNormal, double * q_mNormal, double * q_fNormal, modelParamType & argModelParam, double & eQALY, double & eLY, double & eCost, double & discountFactAdjQ, double & discountFactAdjC,
	int delayedWeek, int horz)
{
	double discountFactAdjDALY = pow((1 + DISCOUNT_FACT_DALY), (1 / CYCLES_PER_YEAR)) - 1; //weekly rate = (1 + annual rate)^(1/52) - 1


	//modelParamType argModelParam=theModelParam;

	// Modified by QC, 9/30/2014
	double totalHorizon;
	if (delayedWeek == -1) {
		totalHorizon = (TIME_HORIZON * CYCLES_PER_YEAR); //static_cast<int>(TIME_HORIZON * CYCLES_PER_YEAR)+1;
	}
	else {
		totalHorizon = delayedWeek;
	}

	if (horz != -1) { totalHorizon = horz; }// added @ 2016/12/1

	if (pat.state != s_Death && pat.flagCured == false)
		while (pat.currentAge <= 100 && cycleNum < totalHorizon) // [QC modified 9/30/2014] Maximum patient age = 120; if patient is not treated and alive
		{

			_simCounter.txCountNone[cycleNum]++;
			rand_pr = lcgrand(3); //static_cast<float>(rand()) / static_cast<float>(RAND_MAX); //generates random number between 0 and 1; RAND_MAX = 32767


			if (pat.gender == 'M') { //select annual mortality rate and qol according to patient's sex
				pr_Death = pr_mDeathAllCause[static_cast<int>(pat.currentAge)];
				q_AgeNormal = q_mNormal[static_cast<int>(pat.currentAge)];

			}
			else {
				pr_Death = pr_fDeathAllCause[static_cast<int>(pat.currentAge)];
				q_AgeNormal = q_fNormal[static_cast<int>(pat.currentAge)];
			}

			pr_Death = funcConvertProb(pr_Death); //convert annual transition probability according to cycle length

			if (EXPECTED_LIFE_ONLY == 1)
				q_AgeNormal = 1;

			/*
			if(PROGRESSION_AGE_GENDER == 1){
			argModelParam._transData.pr_F0_F1 = funcFibrosisProgression (pat.currentAge, pat.gender);
			argModelParam._transData.pr_F1_F2 = funcFibrosisProgression (pat.currentAge, pat.gender);
			argModelParam._transData.pr_F2_F3 = funcFibrosisProgression (pat.currentAge, pat.gender);
			argModelParam._transData.pr_F3_CoCirr = funcFibrosisProgression (pat.currentAge, pat.gender);
			}
			*/

			if (PROGRESSION_AGE_GENDER == 1 && PSA_OPTION == 0) {
				// [QC] Get annual risk
				argModelParam._transData.pr_F0_F1 = argModelParam._transData.funcFibrosisF0F1(pat.currentAge, pat.gender, pat.genotype);
				argModelParam._transData.pr_F1_F2 = argModelParam._transData.funcFibrosisF1F2(pat.currentAge, pat.gender, pat.genotype);
				argModelParam._transData.pr_F2_F3 = argModelParam._transData.funcFibrosisF2F3(pat.currentAge, pat.gender, pat.genotype);
				argModelParam._transData.pr_F3_CoCirr = argModelParam._transData.funcFibrosisF3F4(pat.currentAge, pat.gender, pat.genotype);

				// [QC] Conver to cycle-risk
				argModelParam._transData.pr_F0_F1 = funcConvertProb(argModelParam._transData.pr_F0_F1); //convert annual transition probability according to cycle length
				argModelParam._transData.pr_F1_F2 = funcConvertProb(argModelParam._transData.pr_F1_F2); //convert annual transition probability according to cycle length
				argModelParam._transData.pr_F2_F3 = funcConvertProb(argModelParam._transData.pr_F2_F3); //convert annual transition probability according to cycle length
				argModelParam._transData.pr_F3_CoCirr = funcConvertProb(argModelParam._transData.pr_F3_CoCirr); //convert annual transition probability according to cycle length
			}

			// print out the current state



			switch (pat.state)
			{
			case s_F0:
				if (rand_pr <= argModelParam._transData.pr_F0_F1) {
					pat.state = s_F1;
					eQALY += argModelParam._qolData.q_F1 * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F1 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountF1[cycleNum]++;
					_simCounter.prevStateCountF1[cycleNum]++;


				}
				else if (rand_pr > argModelParam._transData.pr_F0_F1 && rand_pr < argModelParam._transData.pr_F0_F1 + pr_Death) {
					pat.state = s_Death;
					//  eQALY += 0 * CYCLE_LENGTH;
				}
				else {    //remain in s_F0
					eQALY += argModelParam._qolData.q_F0 * q_AgeNormal* (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F0 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountF0[cycleNum]++;

					//_tempOutf<<argModelParam._qolData.q_F0 <<endl<< q_AgeNormal <<endl<<(CYCLE_LENGTH / CYCLES_PER_YEAR)  <<endl<< discountFactAdjQ<<endl<<cycleNum<<endl << pow((1 + discountFactAdjQ),cycleNum)<<endl;


				}

				break;

			case s_F1:
				if (rand_pr <= argModelParam._transData.pr_F1_F2) { //move to s_F2
					pat.state = s_F2;
					eQALY += argModelParam._qolData.q_F2 * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F2 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountF2[cycleNum]++;
					_simCounter.prevStateCountF2[cycleNum]++;
				}
				else if (rand_pr > argModelParam._transData.pr_F1_F2 && rand_pr < argModelParam._transData.pr_F1_F2 + pr_Death) { //move to s_Death
					pat.state = s_Death;
				}
				else {  //remain in s_F1
					eQALY += argModelParam._qolData.q_F1 * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum);
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F1 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountF1[cycleNum]++;
				}
				break;

			case s_F2:
				if (rand_pr <= argModelParam._transData.pr_F2_F3) { //move to s_F3
					pat.state = s_F3;
					eQALY += argModelParam._qolData.q_F3 * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F3 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountF3[cycleNum]++;
					_simCounter.prevStateCountF3[cycleNum]++;
				}
				else if (rand_pr > argModelParam._transData.pr_F2_F3 && rand_pr < argModelParam._transData.pr_F2_F3 + pr_Death) { //move to s_Death
					pat.state = s_Death;
				}
				else {  //remain in s_F2
					eQALY += argModelParam._qolData.q_F2 * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F2 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountF2[cycleNum]++;
				}
				break;

			case s_F3:
				if (rand_pr <= argModelParam._transData.pr_F3_CoCirr) { //move to s_CoCirr
					pat.state = s_CoCirr;
					eQALY += argModelParam._qolData.q_CoCirr * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_CoCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountCoCirr[cycleNum]++;
					_simCounter.prevStateCountCoCirr[cycleNum]++;
				}
				else if (rand_pr > argModelParam._transData.pr_F3_CoCirr && rand_pr <= argModelParam._transData.pr_F3_CoCirr + argModelParam._transData.pr_F3_HCC) { //move to s_HCC
					pat.state = s_HCC;
					eQALY += argModelParam._qolData.q_HCC * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountHCC[cycleNum]++;
					_simCounter.prevStateCountHCC[cycleNum]++;
				}
				else if (rand_pr > argModelParam._transData.pr_F3_CoCirr + argModelParam._transData.pr_F3_HCC && rand_pr < argModelParam._transData.pr_F3_CoCirr + argModelParam._transData.pr_F3_HCC + pr_Death) { //move to s_Death
					pat.state = s_Death;
				}
				else {  //remain in s_F3
					eQALY += argModelParam._qolData.q_F3 * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_F3 * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountF3[cycleNum]++;
				}
				break;

			case s_CoCirr:
				if (rand_pr <= argModelParam._transData.pr_CoCirr_DeCirr) { //move to s_DeCirr
					pat.state = s_DeCirr;
					eQALY += argModelParam._qolData.q_DeCirr * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_DeCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountDeCirr[cycleNum]++;
					_simCounter.prevStateCountDeCirr[cycleNum]++;
				}
				else if (rand_pr > argModelParam._transData.pr_CoCirr_DeCirr && rand_pr <= argModelParam._transData.pr_CoCirr_DeCirr + argModelParam._transData.pr_CoCirr_HCC) { //move to s_HCC
					pat.state = s_HCC;
					eQALY += argModelParam._qolData.q_HCC * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountHCC[cycleNum]++;
					_simCounter.prevStateCountHCC[cycleNum]++;
				}
				else if (rand_pr > argModelParam._transData.pr_CoCirr_DeCirr + argModelParam._transData.pr_CoCirr_HCC
					&& rand_pr <= argModelParam._transData.pr_CoCirr_DeCirr + argModelParam._transData.pr_CoCirr_HCC + pr_Death) { //move to s_Death
					pat.state = s_Death;
				}
				else {  //remain in s_CoCirr
					eQALY += argModelParam._qolData.q_CoCirr * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_CoCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountCoCirr[cycleNum]++;
				}

				//// [QC updated 10/5/2014]
				//if(rand_pr <= pr_Death){ //move to s_Death
				//	pat.state = s_Death;  
				//}
				//if (rand_pr > pr_Death && rand_pr <= pr_Death+ argModelParam._transData.pr_CoCirr_DeCirr){ //move to s_DeCirr
				//	pat.state = s_DeCirr;
				//	eQALY += argModelParam._qolData.q_DeCirr * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ),cycleNum);
				//	eCost += argModelParam._costData.c_DeCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC),cycleNum);
				//	_simCounter.incidentStateCountDeCirr[cycleNum]++;
				//	_simCounter.prevStateCountDeCirr[cycleNum]++;
				//}
				//else if (rand_pr >  pr_Death+ argModelParam._transData.pr_CoCirr_DeCirr && rand_pr <=  pr_Death+ argModelParam._transData.pr_CoCirr_DeCirr + argModelParam._transData.pr_CoCirr_HCC){ //move to s_HCC
				//	pat.state = s_HCC;
				//	eQALY += argModelParam._qolData.q_HCC * q_AgeNormal * (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ),cycleNum);
				//	eCost += argModelParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC),cycleNum);
				//	_simCounter.incidentStateCountHCC[cycleNum]++;
				//	_simCounter.prevStateCountHCC[cycleNum]++;
				//}

				//else {  //remain in s_CoCirr
				//	eQALY += argModelParam._qolData.q_CoCirr * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ),cycleNum);
				//	eCost += argModelParam._costData.c_CoCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC),cycleNum);
				//	_simCounter.prevStateCountCoCirr[cycleNum]++;
				//}
				break;

			case s_DeCirr:
				pr_Death = pr_Death*(1 - argModelParam._transData.pr_DeCirr_DeathLiv); //adjust natural probability of death to factor out death due to DeCirr
				if (pat.flagDeCirr == 0) { //during first visit to this state, flag the patient and increment the counter
					pat.flagDeCirr = 1;
					_simCounter.countDeCirr++;
				}
				if (pat.numCycleDeCirr < 52) { //first year in DeCirr 
					if (rand_pr <= argModelParam._transData.pr_DeCirr_HCC) { //move to s_HCC
						pat.state = s_HCC;
						eQALY += argModelParam._qolData.q_HCC * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
						eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
						_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
						eCost += argModelParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
						_simCounter.incidentStateCountHCC[cycleNum]++;
						_simCounter.prevStateCountHCC[cycleNum]++;

						if (pat.flagHCC == 0) { //during first visit to this state, flag the patient and increment the counter
							pat.flagHCC = 1;
							_simCounter.countHCC++;
						}
					}
					else if (rand_pr > argModelParam._transData.pr_DeCirr_HCC && rand_pr <= argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr && pat.currentAge <= maxLivTrAge) { //move to s_LivTr
						pat.state = s_LivTr;
						eQALY += argModelParam._qolData.q_LivTr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
						eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
						_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
						_simCounter.incidentStateCountLivTr[cycleNum]++;
						_simCounter.prevStateCountLivTr[cycleNum]++;
						pat.numCycleLivTr++; //increment the number of cycles spent in liver transplant state
						if (pat.flagLivTr == 0) { //during first visit to this state, flag the patient and increment the counter
							pat.flagLivTr = 1;
							_simCounter.countLivTr++;
							eCost += argModelParam._costData.c_LivTr / pow((1 + discountFactAdjC), cycleNum); //one time cost associated with liver transplant (not dependent on time spent in this state)
						}
					}
					else if (rand_pr > argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr
						&& rand_pr <= argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr + argModelParam._transData.pr_DeCirr_DeathLiv) { //move to s_Death due to liver complications
						pat.state = s_Death;
						if (pat.flagDeathLiv == 0) { //during first visit to this state, flag the patient and increment the counter
							pat.flagDeathLiv = 1;
							_simCounter.countDeathLiv++;
							_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);
							_simCounter.incidentDeathLiv[cycleNum]++;
						}
					}

					else if (rand_pr > argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr + argModelParam._transData.pr_DeCirr_DeathLiv
						&& rand_pr <= argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr + argModelParam._transData.pr_DeCirr_DeathLiv + pr_Death) { //move to s_Death
						pat.state = s_Death;
					}

					else {  //remain in s_DeCirr
						eQALY += argModelParam._qolData.q_DeCirr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
						eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
						_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
						eCost += argModelParam._costData.c_DeCirr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
						_simCounter.prevStateCountDeCirr[cycleNum]++;
					}
				}

				if (pat.numCycleDeCirr == 52) { // if alive after 1st year in DeCirr, move to s_DeCirr1yrPlus
					_simCounter.incidentStateCountDeCirr1yrPlus[cycleNum]++;
					pat.state = s_DeCirr1yrPlus;
					eQALY += argModelParam._qolData.q_DeCirr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_DeCirr1yrPlus * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountDeCirr1yrPlus[cycleNum]++;
				}

				pat.numCycleDeCirr++;

				//// [QC changed on 2/25/2015]------------
				//pat.numCycleDeCirr++;				
				//if (pat.numCycleDeCirr == 52){ // if alive after 1st year in DeCirr, move to s_DeCirr1yrPlus
				//	pat.state = s_DeCirr1yrPlus;
				//}
				//// --- END OF CHANGE -------------------



				break;

			case s_DeCirr1yrPlus:
				pr_Death = pr_Death*(1 - argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv); //adjust natural probability of death to factor out death due to DeCirr
				if (rand_pr <= argModelParam._transData.pr_DeCirr_HCC) { //move to s_HCC
					pat.state = s_HCC;
					eQALY += argModelParam._qolData.q_HCC * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.incidentStateCountHCC[cycleNum]++;
					_simCounter.prevStateCountHCC[cycleNum]++;

					if (pat.flagHCC == 0) { //during first visit to this state, flag the patient and increment the counter
						pat.flagHCC = 1;
						_simCounter.countHCC++;
					}
				}
				else if (rand_pr > argModelParam._transData.pr_DeCirr_HCC && rand_pr <= argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr && pat.currentAge <= maxLivTrAge) { //move to s_LivTr
					pat.state = s_LivTr;
					eQALY += argModelParam._qolData.q_LivTr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					_simCounter.incidentStateCountLivTr[cycleNum]++;
					_simCounter.prevStateCountLivTr[cycleNum]++;
					pat.numCycleLivTr++; //increment the number of cycles spent in liver transplant state
					if (pat.flagLivTr == 0) { //during first visit to this state, flag the patient and increment the counter
						pat.flagLivTr = 1;
						_simCounter.countLivTr++;
						eCost += argModelParam._costData.c_LivTr / pow((1 + discountFactAdjC), cycleNum); //one time cost associated with liver transplant (not dependent on time spent in this state)
					}
				}
				else if (rand_pr > argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr && rand_pr <= argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr + argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv) { //move to s_Death due to liver complications
					pat.state = s_Death;
					if (pat.flagDeathLiv == 0) { //during first visit to this state, flag the patient and increment the counter
						pat.flagDeathLiv = 1;
						_simCounter.countDeathLiv++;
						_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);
						_simCounter.incidentDeathLiv[cycleNum]++;
					}
				}
				else if (rand_pr > argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr + argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv && rand_pr <= argModelParam._transData.pr_DeCirr_HCC + argModelParam._transData.pr_DeCirr_LivTr + argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv + pr_Death) { //move to s_Death
					pat.state = s_Death;
				}

				else {  //remain in s_DeCirr1yrPlus
					eQALY += argModelParam._qolData.q_DeCirr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_DeCirr1yrPlus * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountDeCirr1yrPlus[cycleNum]++;
				}

				pat.numCycleDeCirr++;
				break;

			case s_HCC:
				pr_Death = pr_Death*(1 - argModelParam._transData.pr_HCC_DeathLiv); //adjust natural probability of death to factor out death due to HCC
				if (pat.flagHCC == 0) { //during first visit to this state, flag the patient and increment the counter
					pat.flagHCC = 1;
					_simCounter.countHCC++;
				}
				if (rand_pr <= argModelParam._transData.pr_HCC_LivTr && pat.currentAge <= maxLivTrAge) { //move to s_LivTr
					pat.state = s_LivTr;
					eQALY += argModelParam._qolData.q_LivTr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					_simCounter.incidentStateCountLivTr[cycleNum]++;
					_simCounter.prevStateCountLivTr[cycleNum]++;
					pat.numCycleLivTr++; //increment the number of cycles spent in liver transplant state

					if (pat.flagLivTr == 0) { //during first visit to this state, flag the patient and increment the counter
						pat.flagLivTr = 1;
						_simCounter.countLivTr++;
						eCost += argModelParam._costData.c_LivTr / pow((1 + discountFactAdjC), cycleNum); //one time cost associated with liver transplant (not dependent on time spent in this state)
					}
				}
				else if (rand_pr > argModelParam._transData.pr_HCC_LivTr && rand_pr <= argModelParam._transData.pr_HCC_LivTr + argModelParam._transData.pr_HCC_DeathLiv) { //move to liver related s_Death
					pat.state = s_Death;
					if (pat.flagDeathLiv == 0) { //during first visit to this state, flag the patient and increment the counter
						pat.flagDeathLiv = 1;
						_simCounter.countDeathLiv++;
						_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);
						_simCounter.incidentDeathLiv[cycleNum]++;
					}
				}
				else if (rand_pr > argModelParam._transData.pr_HCC_LivTr + argModelParam._transData.pr_HCC_DeathLiv && rand_pr <= argModelParam._transData.pr_HCC_LivTr + argModelParam._transData.pr_HCC_DeathLiv + pr_Death) { //move to s_Death
					pat.state = s_Death;
				}
				else {  //remain in s_HCC
					eQALY += argModelParam._qolData.q_HCC * q_AgeNormal *(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_HCC * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountHCC[cycleNum]++;
				}
				break;

			case s_LivTr:
				pr_Death = pr_Death*(1 - argModelParam._transData.pr_LivTr_DeathLiv); //adjust natural probability of death to factor out death due to LivTr
				if (pat.flagLivTr == 0) { //during first visit to this state, flag the patient and increment the counter
					pat.flagLivTr = 1;
					_simCounter.countLivTr++;
					_simCounter.incidentStateCountLivTr[cycleNum]++;
					_simCounter.prevStateCountLivTr[cycleNum]++;
				}
				if (pat.numCycleLivTr < 52) { //first year remain in LivTr or move to s_Death (from liver complication or natural causes)
					if (rand_pr <= argModelParam._transData.pr_LivTr_DeathLiv) { //move to liver related s_Death 
						pat.state = s_Death;
						if (pat.flagDeathLiv == 0) { //during first visit to this state, flag the patient and increment the counter
							pat.flagDeathLiv = 1;
							_simCounter.countDeathLiv++;
							_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);
							_simCounter.incidentDeathLiv[cycleNum]++;
						}
					}
					else if (rand_pr > argModelParam._transData.pr_LivTr_DeathLiv && rand_pr <= argModelParam._transData.pr_LivTr_DeathLiv + pr_Death) { //move to s_Death due to natual causes
						pat.state = s_Death;
					}
					else { // remain in s_LivTr
						eQALY += argModelParam._qolData.q_LivTr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
						eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
						_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
						_simCounter.prevStateCountLivTr[cycleNum]++;
						if (pat.flagLivTr == 0) { //during first visit to this state, flag the patient and increment the counter
							pat.flagLivTr = 1;
							_simCounter.countLivTr++;
							eCost += argModelParam._costData.c_LivTr / pow((1 + discountFactAdjC), cycleNum); //one time cost associated with liver transplant (not dependent on time spent in this state)
						}
					}
				}

				if (pat.numCycleLivTr == 52) { // if alive after 1st year of liver transplant, move to s_LivTr1yrPlus
					_simCounter.incidentStateCountLivTr1yrPlus[cycleNum]++;
					pat.state = s_LivTr1yrPlus;
					eQALY += argModelParam._qolData.q_PostLivTr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_PostLivTr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountLivTr1yrPlus[cycleNum]++;

				}
				pat.numCycleLivTr++; //increment the number of cycles spent in liver transplant state

				//// [QC changed on 2/25/2015]------------
				//pat.numCycleLivTr++;				
				//if (pat.numCycleLivTr == 52){ // if alive after 1st year in DeCirr, move to s_DeCirr1yrPlus
				//	pat.state = s_LivTr1yrPlus;
				//}
				//// --- END OF CHANGE -------------------
				break;


			case s_LivTr1yrPlus:
				pr_Death = pr_Death*(1 - argModelParam._transData.pr_LivTr1yrPlus_DeathLiv); //adjust natural probability of death to factor out death due to LivTr1yrPlus
				if (rand_pr <= argModelParam._transData.pr_LivTr1yrPlus_DeathLiv) { //move to liver related s_Death 
					pat.state = s_Death;
					if (pat.flagDeathLiv == 0) { //during first visit to this state, flag the patient and increment the counter
						pat.flagDeathLiv = 1;
						_simCounter.countDeathLiv++;
						_sim_DALY_YLL += argModelParam._dalyData.GetYLL(pat.currentAge, pat.gender);
						_simCounter.incidentDeathLiv[cycleNum]++;
					}
				}
				else if (rand_pr > argModelParam._transData.pr_LivTr1yrPlus_DeathLiv && rand_pr <= argModelParam._transData.pr_LivTr1yrPlus_DeathLiv + pr_Death) { //move to s_Death due to natual causes
					pat.state = s_Death;
				}
				else { // remain in s_LivTr1yrPlus
					eQALY += argModelParam._qolData.q_PostLivTr * q_AgeNormal*(CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjQ), cycleNum);
					eLY += (CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + DISCOUNT_FACT_LY), cycleNum); 
					_sim_DALY_YLD += argModelParam._dalyData.GetYLD(pat.currentAge, pat.state, CYCLE_LENGTH / CYCLES_PER_YEAR) / pow((1 + discountFactAdjDALY), cycleNum); // added @ 20160826
					eCost += argModelParam._costData.c_PostLivTr * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);
					_simCounter.prevStateCountLivTr1yrPlus[cycleNum]++;
				}
				pat.numCycleLivTr++;
				break;

			} //end of switch statement
			//         outFile << patientNum << "\t" << pat.currentAge <<" \t " <<"State: " << pat.state << " \t " <<eQALY << endl;


			if (pat.state == s_Death) { //exit the while loop if patient dies
				_simCounter.incidentDeathCount[cycleNum]++; //increment the counter for death incident
				break;
			}

			eCost += argModelParam._costData.c_background * CYCLE_LENGTH / pow((1 + discountFactAdjC), cycleNum);

			pat.currentAge = pat.currentAge + CYCLE_LENGTH / CYCLES_PER_YEAR; //advance patient's age in the model if patient is alive in the given cycle
			cycleNum++;
			_simCounter.aliveCount[cycleNum]++; //increment counter for patient alive in the current cycle


		} //end of while loop
	return 0;
}


int HepCSim::ConvertParameters(modelParamType & argModelParam)
{
	//***************************************************************************************//
	//if estimating life years (instead of QALYs)

	if (EXPECTED_LIFE_ONLY == 1) {
		argModelParam._qolData.q_acute = 1.0;//0.98;
		argModelParam._qolData.q_F0 = 1.0;//0.98;
		argModelParam._qolData.q_F1 = 1.0;//0.94;
		argModelParam._qolData.q_F2 = 1.0;//0.94;
		argModelParam._qolData.q_F3 = 1.0;//0.94;
		argModelParam._qolData.q_CoCirr = 1.0;//0.89;
		argModelParam._qolData.q_DeCirr = 1.0;//0.81;
		argModelParam._qolData.q_HCC = 1.0;//0.81;
		argModelParam._qolData.q_LivTr = 1.0;//0.86;
		argModelParam._qolData.q_PostLivTr = 1.0;//
		argModelParam._qolData.q_SVR = 1.0;//1;	
		argModelParam._qolData.q_ETR = 1.0;//0.9;
		argModelParam._qolData.q_TX_oSOC = 1.0;//0.9;
		argModelParam._qolData.q_TX_DAA = 1.0;//0.9;
		argModelParam._qolData.q_Dec_Anemia = 1.0;//0.90;
	}

	//***************************************************************************************//
	//to exclude adverse event, the following loop is executed

	if (ADVERSE_EVENTS == 0) {
		argModelParam._transData.pr_AEAnemia_SOC = 0;//0.34;
		argModelParam._transData.pr_AEAnemia_RGT = 0;//0.34;
		argModelParam._transData.pr_AEAnemia_BOC = 0;//0.56;
		argModelParam._transData.pr_AERash_SOC = 0;//0.37; 
		argModelParam._transData.pr_AERash_BOC = 0;//0.40; 
		argModelParam._transData.pr_AEAnemia_TxExp_SOC = 0;// = 16.0/60.0; 
		argModelParam._transData.pr_AEAnemia_TxExp_RGT = 0;// = 70.0/162.; 
		argModelParam._transData.pr_AEAnemia_TxExp_PRB48 = 0;// = 74.0/161.;  
		argModelParam._costData.c_Epo = 0;//2150
		argModelParam._qolData.q_Dec_Anemia = 1.0;//0.90;
	}

	//***************************************************************************************//	
	//CHOOSE either one of the following conversions********************
	//***************************************************************************************//	

	//convert annual rate according to cycle length probabilities
	argModelParam._transData.pr_FibProg_M50under = funcConvertRate(argModelParam._transData.pr_FibProg_M50under);
	argModelParam._transData.pr_FibProg_M5059 = funcConvertRate(argModelParam._transData.pr_FibProg_M5059);
	argModelParam._transData.pr_FibProg_M6069 = funcConvertRate(argModelParam._transData.pr_FibProg_M6069);
	argModelParam._transData.pr_FibProg_M70plus = funcConvertRate(argModelParam._transData.pr_FibProg_M70plus);
	argModelParam._transData.pr_FibProg_F50under = funcConvertRate(argModelParam._transData.pr_FibProg_F50under);
	argModelParam._transData.pr_FibProg_F5059 = funcConvertRate(argModelParam._transData.pr_FibProg_F5059);
	argModelParam._transData.pr_FibProg_F6069 = funcConvertRate(argModelParam._transData.pr_FibProg_F6069);
	argModelParam._transData.pr_FibProg_F7079 = funcConvertRate(argModelParam._transData.pr_FibProg_F7079);
	argModelParam._transData.pr_FibProg_F80plus = funcConvertRate(argModelParam._transData.pr_FibProg_F80plus);
	//    argModelParam._transData.pr_CoCirr_DeCirr = funcConvertRate(argModelParam._transData.pr_CoCirr_DeCirr);
	//    argModelParam._transData.pr_CoCirr_HCC = funcConvertRate(argModelParam._transData.pr_CoCirr_HCC);
	//    argModelParam._transData.pr_DeCirr_HCC = funcConvertRate(argModelParam._transData.pr_DeCirr_HCC);

	//convert annual transition probability according to cycle length probabilities	
	argModelParam._transData.pr_F0_F1 = funcConvertProb(argModelParam._transData.pr_F0_F1);
	argModelParam._transData.pr_F1_F2 = funcConvertProb(argModelParam._transData.pr_F1_F2);
	argModelParam._transData.pr_F2_F3 = funcConvertProb(argModelParam._transData.pr_F2_F3);
	argModelParam._transData.pr_F3_CoCirr = funcConvertProb(argModelParam._transData.pr_F3_CoCirr);
	argModelParam._transData.pr_F3_HCC = funcConvertProb(argModelParam._transData.pr_F3_HCC);
	argModelParam._transData.pr_CoCirr_DeCirr = funcConvertProb(argModelParam._transData.pr_CoCirr_DeCirr);
	argModelParam._transData.pr_CoCirr_HCC = funcConvertProb(argModelParam._transData.pr_CoCirr_HCC);
	argModelParam._transData.pr_DeCirr_HCC = funcConvertProb(argModelParam._transData.pr_DeCirr_HCC);
	argModelParam._transData.pr_DeCirr_LivTr = funcConvertProb(argModelParam._transData.pr_DeCirr_LivTr);
	argModelParam._transData.pr_DeCirr_DeathLiv = funcConvertProb(argModelParam._transData.pr_DeCirr_DeathLiv);
	argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = funcConvertProb(argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv);
	argModelParam._transData.pr_HCC_LivTr = funcConvertProb(argModelParam._transData.pr_HCC_LivTr);
	argModelParam._transData.pr_HCC_DeathLiv = funcConvertProb(argModelParam._transData.pr_HCC_DeathLiv);
	argModelParam._transData.pr_LivTr_DeathLiv = funcConvertProb(argModelParam._transData.pr_LivTr_DeathLiv);
	argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = funcConvertProb(argModelParam._transData.pr_LivTr1yrPlus_DeathLiv);
	argModelParam._transData.pr_SVR_CoCirr_DeCirr = funcConvertProb(argModelParam._transData.pr_SVR_CoCirr_DeCirr);
	argModelParam._transData.pr_SVR_CoCirr_HCC = funcConvertProb(argModelParam._transData.pr_SVR_CoCirr_HCC);


	//convert costs according to cycle length  	

	argModelParam._costData.c_F0 = argModelParam._costData.c_F0 / CYCLES_PER_YEAR;
	argModelParam._costData.c_F1 = argModelParam._costData.c_F1 / CYCLES_PER_YEAR;
	argModelParam._costData.c_F2 = argModelParam._costData.c_F2 / CYCLES_PER_YEAR;
	argModelParam._costData.c_F3 = argModelParam._costData.c_F3 / CYCLES_PER_YEAR;
	argModelParam._costData.c_CoCirr = argModelParam._costData.c_CoCirr / CYCLES_PER_YEAR;
	argModelParam._costData.c_DeCirr = argModelParam._costData.c_DeCirr / CYCLES_PER_YEAR;
	argModelParam._costData.c_DeCirr1yrPlus = argModelParam._costData.c_DeCirr1yrPlus / CYCLES_PER_YEAR;
	argModelParam._costData.c_HCC = argModelParam._costData.c_HCC / CYCLES_PER_YEAR;
	argModelParam._costData.c_LivTr = argModelParam._costData.c_LivTr; // this is one-time lump-sum cost of liver transpalnt
	argModelParam._costData.c_PostLivTr = argModelParam._costData.c_PostLivTr / CYCLES_PER_YEAR;
	argModelParam._costData.c_ETR = argModelParam._costData.c_ETR / CYCLES_PER_YEAR;
	argModelParam._costData.c_SE_Boc = argModelParam._costData.c_SE_Boc / CYCLES_PER_YEAR;
	argModelParam._costData.c_SVR = argModelParam._costData.c_SVR / CYCLES_PER_YEAR;
	argModelParam._costData.c_PEG = argModelParam._costData.c_PEG; //weekly cost of peginterferon+Ribovarin; based on MAY 2008 WAC prices
	argModelParam._costData.c_BOC = argModelParam._costData.c_BOC; //weekly cost of Bocepravir

	// added 2017/03/09
	argModelParam._costData.c_background = argModelParam._costData.c_background / CYCLES_PER_YEAR;
	return 0;

}

int HepCSim::ReadSmpParamForPSA(ifstream & inf, modelParamType & argModelParam)
{

	inf >> argModelParam._cohortData.baseState >> argModelParam._cohortData.baseGender >> argModelParam._cohortData.baseAge
		>> argModelParam._qolData.q_F0 >> argModelParam._qolData.q_F1 >> argModelParam._qolData.q_F2 >> argModelParam._qolData.q_F3 >> argModelParam._qolData.q_CoCirr
		>> argModelParam._qolData.q_DeCirr >> argModelParam._qolData.q_HCC >> argModelParam._qolData.q_LivTr >> argModelParam._qolData.q_PostLivTr
		>> argModelParam._qolData.q_SVR >> argModelParam._qolData.q_Dec_Anemia >> argModelParam._qolData.q_TX_oSOC >> argModelParam._qolData.q_TX_DAA
		>> argModelParam._costData.c_F0 >> argModelParam._costData.c_F1 >> argModelParam._costData.c_F2 >> argModelParam._costData.c_F3 >> argModelParam._costData.c_CoCirr
		>> argModelParam._costData.c_DeCirr >> argModelParam._costData.c_DeCirr1yrPlus >> argModelParam._costData.c_HCC >> argModelParam._costData.c_LivTr >> argModelParam._costData.c_PostLivTr
		>> argModelParam._transData.pr_F0_F1 >> argModelParam._transData.pr_F1_F2 >> argModelParam._transData.pr_F2_F3 >> argModelParam._transData.pr_F3_CoCirr
		>> argModelParam._transData.pr_CoCirr_DeCirr >> argModelParam._transData.pr_CoCirr_HCC >> argModelParam._transData.pr_DeCirr_HCC
		>> argModelParam._transData.pr_DeCirr_LivTr >> argModelParam._transData.pr_DeCirr_DeathLiv >> argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv
		>> argModelParam._transData.pr_HCC_LivTr >> argModelParam._transData.pr_HCC_DeathLiv >> argModelParam._transData.pr_LivTr_DeathLiv
		>> argModelParam._transData.pr_LivTr1yrPlus_DeathLiv >> argModelParam._transData.pr_SVR_CoCirr_DeCirr >> argModelParam._transData.pr_SVR_CoCirr_HCC
		>> argModelParam._transData.pr_SVR_Delta_oSOC >> argModelParam._transData.pr_SVR_Delta_DAA;
	return 0;
}

int HepCSim::OutputBaseCase(ofstream & outf)
{
	/*"Scenario" << "\t" << "Counter" << "\t"
	<< "QALY" << "\t" << "COST" << "\t" << "COUNT_DECIRR" << "\t"
	<< "COUNT_HCC" << "\t" << "COUNT_LIVTR" << "\t" << "COUNT_DEATHLIV"
	<< "\t" << "COST_TX" << endl;*/
	outf << _sim_QALY << "\t" << _sim_cost << "\t"
		<< _simCounter.countDeCirr << "\t" << _simCounter.countHCC << "\t" << _simCounter.countLivTr << "\t" << _simCounter.countDeathLiv << "\t" << _sim_cost_tx <<"\t"
		<< _sim_LY << endl;
	return 0;
}

int HepCSim::PSA_initializeSampler(const vector<string> & vecVarName, const vector<string> & vecDistrName, const vector<double> & vecP1, const vector<double> & vecP2)
{
	seed_type s = _psa_seed();
	for (int k = 0; k < vecVarName.size(); k++) {
		CRndVarGen rvGen;
		seed_type theSeed = PSA_SAME_SEED_FOR_ALL_PARAM ? s : (_psa_seed());
		TYPE_RANDOM_DISTR  theType;
		if ("Beta" == vecDistrName[k])
			theType = BETA;
		else if ("Gamma" == vecDistrName[k])
			theType = GAMMA;
		else if ("Uniform" == vecDistrName[k])
			theType = UNIFREAL;
		else
			ExitWithMsg("[Error]: PSA_initializeSampler: unknown distr type in the input file: " + vecDistrName[k]);
		rvGen.Initialize(theSeed, theType, vecP1[k], vecP2[k]);
		_psa_sampler[vecVarName[k]] = rvGen;



	}
	return 0;


}

int HepCSim::PSA_initializeSampler(const psaDistrType & argPSADistr)
{
	return PSA_initializeSampler(argPSADistr._distrPSA_varName, argPSADistr._distrPSA_distrName, argPSADistr._distrPSA_param1, argPSADistr._distrPSA_param2);
}

int HepCSim::PSA_sampleModelParamValue(modelParamType & argModelParam)
{
	for (map<string, CRndVarGen>::const_iterator it = _psa_sampler.begin(); it != _psa_sampler.end(); it++) {
		string varName = it->first;
		if ("pF0_F1" == varName) {
			argModelParam._transData.pr_F0_F1 = _psa_sampler[varName].GetRV();
		}
		else if ("pF1_F2_SA" == varName) {
			argModelParam._transData.pr_F1_F2 = _psa_sampler[varName].GetRV();
		}
		else if ("pF2_F3_SA" == varName) {
			argModelParam._transData.pr_F2_F3 = _psa_sampler[varName].GetRV();
		}
		else if ("pF3_F4_SA" == varName) {
			argModelParam._transData.pr_F3_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("pF4_DC_SA" == varName) {
			argModelParam._transData.pr_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("pF4_HCC_SA" == varName) {
			argModelParam._transData.pr_CoCirr_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("pDC_HCC_SA" == varName) {
			argModelParam._transData.pr_DeCirr_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("pDC_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_DeCirr_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_dc_cyc_1_SA" == varName) {
			argModelParam._transData.pr_DeCirr_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_dc_cyc_2_SA" == varName) {
			argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pHCC_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_HCC_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_hcc_cyc_SA" == varName) {
			argModelParam._transData.pr_HCC_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_LivTr_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pMort_Post_Liv_Transpl_SA" == varName) {
			argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
		}
		else if ("pr_SVR_CoCirr_DeCirr" == varName) {
			argModelParam._transData.pr_SVR_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("pr_SVR_CoCirr_HCC" == varName) {
			argModelParam._transData.pr_SVR_CoCirr_HCC = _psa_sampler[varName].GetRV();

		}
		else if ("qolType::q_Acute" == varName) {
			argModelParam._qolData.q_acute = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F0" == varName) {
			argModelParam._qolData.q_F0 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F1" == varName) {
			argModelParam._qolData.q_F1 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F2" == varName) {
			argModelParam._qolData.q_F2 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_F3" == varName) {
			argModelParam._qolData.q_F3 = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_CoCirr" == varName) {
			argModelParam._qolData.q_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_DeCirr" == varName) {
			argModelParam._qolData.q_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_HCC" == varName) {
			argModelParam._qolData.q_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_LivTr" == varName) {
			argModelParam._qolData.q_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_PostLivTr" == varName) {
			argModelParam._qolData.q_PostLivTr = _psa_sampler[varName].GetRV();
		}
		else if ("qolType::q_SVR" == varName) {
			argModelParam._qolData.q_SVR = _psa_sampler[varName].GetRV();

		}
		else if ("Anemia_multiplier" == varName) {
			argModelParam._qolData.q_Dec_Anemia = _psa_sampler[varName].GetRV();
		}
		else if ("Therapy-related_multiplier-oSOC" == varName) {
			argModelParam._qolData.q_TX_oSOC = _psa_sampler[varName].GetRV();
		}
		else if ("Therapy-related_multiplier-DAA" == varName) {
			argModelParam._qolData.q_TX_DAA = _psa_sampler[varName].GetRV();

		}
		else if ("costType::c_Acute" == varName) {
			argModelParam._costData.c_acute = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F0" == varName) {
			argModelParam._costData.c_F0 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F1" == varName) {
			argModelParam._costData.c_F1 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F2" == varName) {
			argModelParam._costData.c_F2 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_F3" == varName) {
			argModelParam._costData.c_F3 = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_CoCirr" == varName) {
			argModelParam._costData.c_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_DeCirr" == varName) {
			argModelParam._costData.c_DeCirr = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_DeCirr1yrPlus" == varName) {
			argModelParam._costData.c_DeCirr1yrPlus = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_HCC" == varName) {
			argModelParam._costData.c_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_LivTr" == varName) {
			argModelParam._costData.c_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_PostLivTr" == varName) {
			argModelParam._costData.c_PostLivTr = _psa_sampler[varName].GetRV();

		}
		else if ("pSVR_Delta_DAA" == varName) {
			argModelParam._transData.pr_SVR_Delta_DAA = _psa_sampler[varName].GetRV();
		}
		else if ("pSVR_Delta_oSOC" == varName) {
			argModelParam._transData.pr_SVR_Delta_oSOC = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F0" == varName) {
			argModelParam._dalyData._dw_f0 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F1" == varName) {
			argModelParam._dalyData._dw_f1 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F2" == varName) {
			argModelParam._dalyData._dw_f2 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_F3" == varName) {
			argModelParam._dalyData._dw_f3 = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_CoCirr" == varName) {
			argModelParam._dalyData._dw_CoCirr = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_DeCirr" == varName) {
			double v = _psa_sampler[varName].GetRV();
			argModelParam._dalyData._dw_DeCirr = v;
			argModelParam._dalyData._dw_DeCirr1yrPlus = v;
		}
		else if ("dalyType::dw_HCC" == varName) {
			argModelParam._dalyData._dw_HCC = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_LivTr" == varName) {
			argModelParam._dalyData._dw_LivTr = _psa_sampler[varName].GetRV();
		}
		else if ("dalyType::dw_PostLivTr" == varName) {
			argModelParam._dalyData._dw_LivTr1yrPlus = _psa_sampler[varName].GetRV();
		}

		else if ("costType::c_testing_preTx" == varName) {
			argModelParam._costData.c_testing_preTx = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_testing_postTx" == varName) {
			argModelParam._costData.c_testing_postTx = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_testing_RNA" == varName) {
			argModelParam._costData.c_testing_RNA = _psa_sampler[varName].GetRV();
		}
		else if ("costType::c_testing_genotype" == varName) {
			argModelParam._costData.c_testing_genotype = _psa_sampler[varName].GetRV();
		}
		else if ("pAcuteSelfClearance" == varName) {
			argModelParam._transData._pr_acute_spont_clearance = _psa_sampler[varName].GetRV();
		}
		else if ("pr_SVR" == varName) {
			argModelParam._transData._pr_svr_for_SA = _psa_sampler[varName].GetRV();
		}
		else if ("pr_lost_to_followup" == varName) {
			double rv= 1.0 - _psa_sampler[varName].GetRV();
			argModelParam._transData._pr_chronic_prob_follow_up = rv;
			argModelParam._transData._pr_acute_prob_complete_tx = rv;
		}


		else {
			ExitWithMsg("[Error] SampleModelParamValue(): Unknown parameter " + varName);
		}
	}

	if (PSA_SAME_QOL_F0_F3) {
		argModelParam._qolData.q_F1 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F2 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F3 = argModelParam._qolData.q_F0;


		argModelParam._costData.c_F1 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F2 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F3 = argModelParam._costData.c_F0;

		argModelParam._dalyData._dw_f1 = argModelParam._dalyData._dw_f0;
		argModelParam._dalyData._dw_f2 = argModelParam._dalyData._dw_f0;
		argModelParam._dalyData._dw_f3 = argModelParam._dalyData._dw_f0;
	}

	return 0;
}


int HepCSim::PSA_sampleModelParamValue_Univariate(modelParamType & argModelParam, string argVarName)
{
	map<string, CRndVarGen>::const_iterator it = _psa_sampler.find(argVarName);
	string varName = it->first;
	if ("pF0_F1" == varName) {
		argModelParam._transData.pr_F0_F1 = _psa_sampler[varName].GetRV();
	}
	else if ("pF1_F2_SA" == varName) {
		argModelParam._transData.pr_F1_F2 = _psa_sampler[varName].GetRV();
	}
	else if ("pF2_F3_SA" == varName) {
		argModelParam._transData.pr_F2_F3 = _psa_sampler[varName].GetRV();
	}
	else if ("pF3_F4_SA" == varName) {
		argModelParam._transData.pr_F3_CoCirr = _psa_sampler[varName].GetRV();
	}
	else if ("pF4_DC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("pF4_HCC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("pDC_HCC_SA" == varName) {
		argModelParam._transData.pr_DeCirr_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("pDC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_DeCirr_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_dc_cyc_1_SA" == varName) {
		argModelParam._transData.pr_DeCirr_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_dc_cyc_2_SA" == varName) {
		argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pHCC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_HCC_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_hcc_cyc_SA" == varName) {
		argModelParam._transData.pr_HCC_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pMort_Post_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = _psa_sampler[varName].GetRV();
	}
	else if ("pr_SVR_CoCirr_DeCirr" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("pr_SVR_CoCirr_HCC" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_HCC = _psa_sampler[varName].GetRV();

	}
	else if ("qolType::q_F0" == varName) {
		argModelParam._qolData.q_F0 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_F1" == varName) {
		argModelParam._qolData.q_F1 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_F2" == varName) {
		argModelParam._qolData.q_F2 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_F3" == varName) {
		argModelParam._qolData.q_F3 = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_CoCirr" == varName) {
		argModelParam._qolData.q_CoCirr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_DeCirr" == varName) {
		argModelParam._qolData.q_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_HCC" == varName) {
		argModelParam._qolData.q_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_LivTr" == varName) {
		argModelParam._qolData.q_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_PostLivTr" == varName) {
		argModelParam._qolData.q_PostLivTr = _psa_sampler[varName].GetRV();
	}
	else if ("qolType::q_SVR" == varName) {
		argModelParam._qolData.q_SVR = _psa_sampler[varName].GetRV();

	}
	else if ("Anemia_multiplier" == varName) {
		argModelParam._qolData.q_Dec_Anemia = _psa_sampler[varName].GetRV();
	}
	else if ("Therapy-related_multiplier-oSOC" == varName) {
		argModelParam._qolData.q_TX_oSOC = _psa_sampler[varName].GetRV();
	}
	else if ("Therapy-related_multiplier-DAA" == varName) {
		argModelParam._qolData.q_TX_DAA = _psa_sampler[varName].GetRV();

	}
	else if ("costType::c_F0" == varName) {
		argModelParam._costData.c_F0 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_F1" == varName) {
		argModelParam._costData.c_F1 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_F2" == varName) {
		argModelParam._costData.c_F2 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_F3" == varName) {
		argModelParam._costData.c_F3 = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_CoCirr" == varName) {
		argModelParam._costData.c_CoCirr = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_DeCirr" == varName) {
		argModelParam._costData.c_DeCirr = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_DeCirr1yrPlus" == varName) {
		argModelParam._costData.c_DeCirr1yrPlus = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_HCC" == varName) {
		argModelParam._costData.c_HCC = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_LivTr" == varName) {
		argModelParam._costData.c_LivTr = _psa_sampler[varName].GetRV();
	}
	else if ("costType::c_PostLivTr" == varName) {
		argModelParam._costData.c_PostLivTr = _psa_sampler[varName].GetRV();

	}
	else if ("pSVR_Delta_DAA" == varName) {
		argModelParam._transData.pr_SVR_Delta_DAA = _psa_sampler[varName].GetRV();
	}
	else if ("pSVR_Delta_oSOC" == varName) {
		argModelParam._transData.pr_SVR_Delta_oSOC = _psa_sampler[varName].GetRV();


	}
	else {
		ExitWithMsg("[Error] SampleModelParamValue(): Unknown parameter " + varName);
	}


	if (PSA_SAME_QOL_F0_F3) {
		argModelParam._qolData.q_F1 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F2 = argModelParam._qolData.q_F0;
		argModelParam._qolData.q_F3 = argModelParam._qolData.q_F0;


		argModelParam._costData.c_F1 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F2 = argModelParam._costData.c_F0;
		argModelParam._costData.c_F3 = argModelParam._costData.c_F0;

	}

	return 0;
}

int HepCSim::PSA_sampleCohort(baseCohortType & argCohort)
{
	// sample gender
	if (_psa_sampler_gender.GetU01() < 0.64) {
		argCohort.baseGender = 'M';
	}
	else {
		argCohort.baseGender = 'F';
	}

	// sample initial state
	vector<stateType> label_state;
	vector<double> count_state;
	label_state.push_back(s_F0);		count_state.push_back(0.1);	//F0
	label_state.push_back(s_F1);		count_state.push_back(0.25);	//F1
	label_state.push_back(s_F2);		count_state.push_back(0.2);	//F2
	label_state.push_back(s_F3);		count_state.push_back(0.2);		//F3
	label_state.push_back(s_CoCirr);		count_state.push_back(0.25);	//F4
	argCohort.baseState = DiscreteDistrSampler(label_state, count_state, _psa_sampler_state);

	// sample initial age

	return 0;
}

int HepCSim::ResetCounters()
{
	/********************* initialize counters *******************************/
	_simCounter.countDeCirr = 0;
	_simCounter.countHCC = 0;
	_simCounter.countLivTr = 0;
	_simCounter.countTxEx = 0;
	_simCounter.countNR = 0;
	_simCounter.countReL = 0;
	_simCounter.countDeathLiv = 0;

	for (int i = 0; i < MAXCYCLE; i++) {
		_simCounter.aliveCount[i] = 0;
		_simCounter.incidentDeathCount[i] = 0;
		_simCounter.incidentDeathLiv[i] = 0;
		_simCounter.txCountDouble[i] = 0;
		_simCounter.txCountTriple[i] = 0;
		_simCounter.txCountNone[i] = 0;

		_simCounter.incidentDisTxDouble[i] = 0;
		_simCounter.incidentDisTxTriple[i] = 0;
		_simCounter.incidentStateCountF0[i] = 0;
		_simCounter.incidentStateCountF1[i] = 0;
		_simCounter.incidentStateCountF2[i] = 0;
		_simCounter.incidentStateCountF3[i] = 0;
		_simCounter.incidentStateCountCoCirr[i] = 0;
		_simCounter.incidentStateCountDeCirr[i] = 0;
		_simCounter.incidentStateCountDeCirr1yrPlus[i] = 0;
		_simCounter.incidentStateCountHCC[i] = 0;
		_simCounter.incidentStateCountLivTr[i] = 0;
		_simCounter.incidentStateCountLivTr1yrPlus[i] = 0;
		_simCounter.incidentStateCountCured[i] = 0;
		_simCounter.incidentStateCountETR[i] = 0;

		_simCounter.prevStateCountF0[i] = 0;
		_simCounter.prevStateCountF1[i] = 0;
		_simCounter.prevStateCountF2[i] = 0;
		_simCounter.prevStateCountF3[i] = 0;
		_simCounter.prevStateCountCoCirr[i] = 0;
		_simCounter.prevStateCountDeCirr[i] = 0;
		_simCounter.prevStateCountDeCirr1yrPlus[i] = 0;
		_simCounter.prevStateCountHCC[i] = 0;
		_simCounter.prevStateCountLivTr[i] = 0;
		_simCounter.prevStateCountLivTr1yrPlus[i] = 0;
		_simCounter.prevStateCountCured[i] = 0;
		_simCounter.prevStateCountETR[i] = 0;

		_simCounter.prevAEAnemia[i] = 0;
		_simCounter.incidentAEAnemia[i] = 0;
	}

	return 0;
}


map<string, SParamTriplet> ReadOneWaySAParamValueRange(string argFile)
{
	map<string, SParamTriplet> ranges;
	//read range of values:
	ifstream inf;
	inf.open(argFile);

	string varName, line;
	while (inf >> varName) {
		if ("//" == varName) {
			getline(inf, line);
			continue;
		}

		SParamTriplet p;
		inf >> p._base;
		inf >> p._lb;
		inf >> p._ub;

		ranges.insert(pair<string, SParamTriplet>(varName, p));

	}

	inf.close();

	return ranges;

}

int ReadOneWaySAParamValueRange(string argFile, vector<string> & listVarName, vector<SParamTriplet> & listVal)
{
	//read range of values:
	ifstream inf;
	inf.open(argFile);

	string varName, line;
	while (inf >> varName) {
		if ("//" == varName) {
			getline(inf, line);
			continue;
		}

		SParamTriplet p;
		inf >> p._base;
		inf >> p._lb;
		inf >> p._ub;

		listVarName.push_back(varName);
		listVal.push_back(p);

	}

	inf.close();

	return 0;
}

