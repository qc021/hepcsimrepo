#include "sim_screening.h"


HCCScreenSim::HCCScreenSim()
{
	_randomSeed = 0;
	_print_results_every_patient = false;
	// initialize
	_tempOutf.open("output_temp.txt", ios::app);
	_tempOutf << fixed << setprecision(9);

	_psa_seed.seed(PSA_SEED);

	_psa_sampler_age.Initialize(_psa_seed(), UNIFREAL);
	_psa_sampler_state.Initialize(_psa_seed(), UNIFREAL);
	_psa_sampler_gender.Initialize(_psa_seed(), UNIFREAL);
	_rnd_pt.Initialize(SIM_SEED, UNIFREAL);
	_rnd_extra.Initialize(SIM_SEED+1, UNIFREAL);	

	_scr_interval_months = 3;
	_cycles_of_one_year = 12.0;
}



int HCCScreenSim::Run(modelParamType argModelParam)
{

	clock_t startclock = clock();	
	ResetCounters();
	_rnd_pt.SetSeed(_randomSeed);

	// convert to probabilities/costs for each cycle
	ConvertParameters(argModelParam);


	//convert annual discount factor according to cycle length

	_cumQALY = 0;
	_cumLY = 0;
	_cumCost = 0;
	_cumCost_screening = 0;
	// start simulation 

	for (long ptIdx = 0; ptIdx < NUM_PATIENTS; ptIdx++) {
		// initiation of simulated patients
		HCCScrPatientType pat;
		InitializePatient(pat, argModelParam);


		for (_sim_cycle = 0; _sim_cycle < (int)(TIME_HORIZON * _cycles_of_one_year); _sim_cycle++) {
			if (pat._stateAbs == scr_abs_death_bg || pat._stateAbs == scr_abs_death_liv) {
				break;
			}

			if (pat._stateAbs == scr_abs_none) {
				// --- check ---
				if (pat._stateFib == scr_fib_inAbs || pat._stateHCC == scr_hcc_inAbs) {
					ExitWithMsg("ERROR @ Run_nat_history: patient's state should not be in absorbing states yet...");
				}


				if (_sim_cycle % _scr_interval_months == 0) {
					// ----------------------------------------------- screening period ------------------------------------------------
					double cTest = argModelParam._costData.hccscr_c_test_surveillance;
					
					// *********************************************************************
					// ****************  For true cancer patients **************************
					// *********************************************************************
					if (pat._stateHCC != scr_hcc_none) {
						// ---------------------------------------------------
						// ---- cancer picked up by the surveillance test ----
						// ---------------------------------------------------
						if (_rnd_pt.GetU01() < argModelParam._transData._sens_surveillance) {
							// --- run diagnostic test
							cTest += argModelParam._costData.hccscr_c_test_diagnostic;
							if (_rnd_pt.GetU01() < argModelParam._transData._sens_diagnostic) {
								// - cancer confirmed by diagnostic test -
								// - initiate treatment

								_cumCost += cTest / pow(1 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);;
								_cumCost_screening += cTest;

								// sample treatment 		
								vector<double> distrTx = argModelParam._hccScreeningData._distr_hcc_tx[GetFibHCCStateStr(pat)];
								if (pat._currentAge >= LIV_TRANSP_MAX_AGE) {
									distrTx[1] = 0.0; // set transplant probability =0 for patient above age 75
								}
								pat._stateAbs = DiscreteDistrSampler(argModelParam._hccScreeningData._list_hcc_tx,
									distrTx, _rnd_pt);
								

								if (pat._stateAbs != scr_abs_trspl) {
									pat._stateFibWhenTreated = pat._stateFib;
									pat._stateHCCWhenTreated = pat._stateHCC;
									pat._stateFib = scr_fib_inAbs;
									pat._stateHCC = scr_hcc_inAbs;
								}
								Run_HCCTreatment(pat, argModelParam);
								continue;
							}
							else {
								// - if not confirmed by the diagnostic test
								// do nothing, continue natural history

							} // end if - diagnostic test
						} 
						// ------------------------------------------------
						// ---- cancer MISSED by the surveillance test ----
						// ------------------------------------------------
						else {
							// do nothing, continue natural history
						}
						_cumCost += cTest / pow(1 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);;
						_cumCost_screening += cTest;
					}
					// *********************************************************************
					// ****************  For non-cancer patients ***************************
					// *********************************************************************
					else {
						if (_rnd_pt.GetU01() < argModelParam._transData._spec_surveillance) {
							// --------- mis-diagnosed as cancer patients, i.e., false positive -----
							// --------- run diagnostic test to confirm (assume 100% specificity of diagnostic test)
							cTest += argModelParam._costData.hccscr_c_test_diagnostic;
						}
						else {
							// do nothing, continue natural history
						}
						_cumCost += cTest / pow(1 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
						_cumCost_screening += cTest;
					}

					Run_NaturalHistory(pat, argModelParam);

				}
				else {
					// ----------------------------------------------- non-screening period -----------------------------------------------------------
					Run_NaturalHistory(pat, argModelParam);

				}
			}
			else {
				// treatment process
				Run_HCCTreatment(pat, argModelParam);
			}
			
			_cumLY += 1.0 / _cycles_of_one_year / pow(1+DISCOUNT_FACT_L, (double)_sim_cycle / _cycles_of_one_year);
			pat._currentAge += 1.0 / _cycles_of_one_year;
		}// end for cycles
		

	}// end for patients

	_sim_QALY = _cumQALY / NUM_PATIENTS;
	_sim_LY = _cumLY / NUM_PATIENTS;
	_sim_cost = _cumCost / NUM_PATIENTS;
	_sim_cost_test = _cumCost_screening / NUM_PATIENTS;
	



	return 0;
}

int HCCScreenSim::Run_NaturalHistory(HCCScrPatientType & argPat, const modelParamType & argModelParam)
{

	// ----- check absorbing state ---------
	if (argPat._stateAbs == scr_abs_death_bg || argPat._stateAbs == scr_abs_death_bg) {
		ExitWithMsg("Error @ Run_nat_history: This patient is already dead....");
	}
	if (argPat._stateAbs == scr_abs_res 
		|| argPat._stateAbs == scr_abs_abl 
		|| argPat._stateAbs == scr_abs_pal
		|| (argPat._stateAbs == scr_abs_trspl && argPat._flag_started_transpl)) {
		if (argPat._stateFib != scr_fib_inAbs || argPat._stateHCC != scr_hcc_inAbs) {
			ExitWithMsg("Error @ Run_nat_history(): Unmatched absorbing states");
		}
	}
	
	// -------------------------------------


	double bgMort = (argPat._gender == 'M') ? argModelParam._hccScreeningData._bgMort_male.lower_bound((int)argPat._currentAge)->second
		: argModelParam._hccScreeningData._bgMort_female.lower_bound((int)argPat._currentAge)->second;
	bgMort = funcConvertProbForCycle(bgMort, _cycles_of_one_year);

	double qolAge = (argPat._gender == 'M') ? argModelParam._hccScreeningData._qol_male.lower_bound((int)argPat._currentAge)->second
		: argModelParam._hccScreeningData._qol_female.lower_bound((int)argPat._currentAge)->second;

	double pProgHCC = GetProbHCCProg(argPat, argModelParam);
	double livMort = GetMortLiv(argPat, argModelParam);

	type_scr_fib_state state_next_fib;
	type_scr_hcc_state state_next_hcc;
	type_scr_absb_state state_next_abs;


	state_next_hcc = StateTransHCC(argPat, argModelParam);
	state_next_abs = scr_abs_none;
	
	// ======================================================================================================================================
	// ==================================================== F3 fibrosis state ===============================================================
	// ======================================================================================================================================
	if (argPat._stateFib == scr_fib_F3) {
		//--- state transitions ---
		double pProgFib = argModelParam._transData.pr_F3_CoCirr;
		state_next_fib = _rnd_pt.GetU01() < pProgFib ? scr_fib_F4 : scr_fib_F3;
		if (argPat._stateFib == state_next_fib && argPat._stateHCC == state_next_hcc) {
			StateTransDeath(pProgHCC, pProgFib, bgMort, livMort, state_next_fib, state_next_hcc, state_next_abs);
		}
		
		// --- Update simulation outcomes ---
		_cumQALY += argModelParam._qolData.q_F3 * qolAge / _cycles_of_one_year / pow(1.0+DISCOUNT_FACT_Q, (double) _sim_cycle/_cycles_of_one_year);
		_cumCost += argModelParam._costData.c_F3 / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
	}

	// ======================================================================================================================================
	// ==================================================== F4 fibrosis state ===============================================================
	// ======================================================================================================================================
	else if (argPat._stateFib == scr_fib_F4) {
		//--- state transitions ---
		double pProgFib = argModelParam._transData.pr_CoCirr_DeCirr;
		state_next_fib = _rnd_pt.GetU01() < pProgFib ? scr_fib_DC : scr_fib_F4;
		if (argPat._stateFib == state_next_fib && argPat._stateHCC == state_next_hcc) {
			StateTransDeath(pProgHCC, pProgFib, bgMort, livMort, state_next_fib, state_next_hcc, state_next_abs);
		}

		// --- Update simulation outcomes ---
		_cumQALY += argModelParam._qolData.q_CoCirr * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);
		_cumCost += argModelParam._costData.c_CoCirr / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);

	}

	// ======================================================================================================================================
	// ====================================================     DC state      ===============================================================
	// ======================================================================================================================================
	else if (argPat._stateFib == scr_fib_DC) {
		//--- state transitions ---
		double pProgFib = 0;
		state_next_fib = scr_fib_DC;
		if (argPat._stateFib == state_next_fib && argPat._stateHCC == state_next_hcc) {
			StateTransDeath(pProgHCC, pProgFib, bgMort, livMort, state_next_fib, state_next_hcc, state_next_abs);
		}

		// --- Update simulation outcomes ---
		_cumQALY += argModelParam._qolData.q_DeCirr * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);
		_cumCost += argModelParam._costData.c_DeCirr / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);

	}

	// ======================================================================================================================================
	// ==================================================== F4-SVR fibrosis state ===========================================================
	// ======================================================================================================================================
	else if (argPat._stateFib == scr_fib_F4_SVR) {
		//--- state transitions ---
		double pProgFib = argModelParam._transData.pr_SVR_CoCirr_DeCirr;
		state_next_fib = _rnd_pt.GetU01() < pProgFib ? scr_fib_DC : scr_fib_F4_SVR;
		if (argPat._stateFib == state_next_fib && argPat._stateHCC == state_next_hcc) {
			StateTransDeath(pProgHCC, pProgFib, bgMort, livMort, state_next_fib, state_next_hcc, state_next_abs);
		}

		// --- Update simulation outcomes ---
		_cumQALY += argModelParam._qolData.q_SVR * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);
		_cumCost += 0;

	}
	// ======================================================================================================================================
	// ==================================================== F3-SVR fibrosis state ===========================================================
	// ======================================================================================================================================

	else if (argPat._stateFib == scr_fib_F3_SVR) {
		//--- state transitions ---
		double pProgFib = 0;
		state_next_fib = scr_fib_F3_SVR;
		if (argPat._stateFib == state_next_fib && argPat._stateHCC == state_next_hcc) {
			StateTransDeath(pProgHCC, pProgFib, bgMort, livMort, state_next_fib, state_next_hcc, state_next_abs);
		}

		// --- Update simulation outcomes ---
		_cumQALY += argModelParam._qolData.q_SVR * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);
		_cumCost += 0;

	}
	else {
		ExitWithMsg("Error @ Run_nat_history: Unmknown argPat._stateFib" + basicToStr((int)argPat._stateFib));
	}

	argPat._stateFib = state_next_fib;
	argPat._stateHCC = state_next_hcc;
	argPat._stateAbs = state_next_abs;


	return 0;
}




int HCCScreenSim::Run_HCCTreatment(HCCScrPatientType & argPat, const modelParamType & argModelParam)
{

	

	double pMortBg = (argPat._gender == 'M') ? argModelParam._hccScreeningData._bgMort_male.lower_bound((int)argPat._currentAge)->second
		: argModelParam._hccScreeningData._bgMort_female.lower_bound((int)argPat._currentAge)->second;
	double qolAge = (argPat._gender == 'M') ? argModelParam._hccScreeningData._qol_male.lower_bound((int)argPat._currentAge)->second
		: argModelParam._hccScreeningData._qol_female.lower_bound((int)argPat._currentAge)->second;

	pMortBg = funcConvertProbForCycle(pMortBg, _cycles_of_one_year);
	double pMortLiv = 0;

	// -------------------------------------------------------------
	// ------------ liver transplant treatment ---------------------
	// -------------------------------------------------------------
	if (argPat._stateAbs == scr_abs_trspl) {
		
		/////////////////////////////////////
		/// IMPORTANT NOTE HERE: 
		/// NEED TO EXPLICIT MODEL THE TRANSPLANT WAITING PROCESS
		/// 3-MONTH transplant probability = 0.364 for DC, and 0.442 for HCC patients
		/// 
		if (argPat._flag_started_transpl){
			if (argPat._cycles_on_HCC_treatment != 0) {
				ExitWithMsg("Error @ Run_Treatment for transplant: _cycles_on_HCC_treatment should not be >0 before transplant has started...");
			}

			if (_rnd_pt.GetU01() >= argModelParam._transData._pr_transpl_hcc) {
				// if still waiting, the patient goes through natural history
				Run_NaturalHistory(argPat, argModelParam);

				// if he doesn't die and progressed to large HCC, then switch to palliative treatment next period
				if (argPat._stateAbs != scr_abs_death_bg && argPat._stateAbs != scr_abs_death_liv) {
					if (argPat._stateHCC == scr_hcc_large || argPat._currentAge >= LIV_TRANSP_MAX_AGE) {
						argPat._stateFibWhenTreated = argPat._stateFib;
						argPat._stateHCCWhenTreated = argPat._stateHCC;
						argPat._stateFib = scr_fib_inAbs;
						argPat._stateHCC = scr_hcc_inAbs;
						argPat._stateAbs = scr_abs_pal;
					}
				}
				return 0;
			}
			else {
				argPat._stateFibWhenTreated = argPat._stateFib;
				argPat._stateHCCWhenTreated = argPat._stateHCC;
				argPat._stateFib = scr_fib_inAbs;
				argPat._stateHCC = scr_hcc_inAbs;

				argPat._flag_started_transpl = true;
			}
		}

		// **** for the first year *********************************
		if (argPat._cycles_on_HCC_treatment < _cycles_of_one_year) {
			_cumQALY += argModelParam._qolData.q_tx_transpl_1y * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);
			_cumCost += argModelParam._costData.hccscr_c_tx_transpl_1y / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
			pMortLiv = argModelParam._transData._mort_tx_transpl_1y;
		}
		// **** for the 1+ years ***********************************
		else {
			_cumQALY += argModelParam._qolData.q_tx_transpl_1yPlus * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);
			_cumCost += argModelParam._costData.hccscr_c_tx_transpl_1yPlus / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
			pMortLiv = argModelParam._transData._mort_tx_transpl_1yPlus;
		}


	}
	// -------------------------------------------------------------
	// ------------ resection treatment ----------------------------
	// -------------------------------------------------------------
	else if (argPat._stateAbs == scr_abs_res) {
		double qol = argPat._cycles_on_HCC_treatment < _cycles_of_one_year ? argModelParam._qolData.qChange_tx_res_1y : argModelParam._qolData.qChange_tx_res_1yPlus;
		qol = qol * GetTxBaselineQoL(argPat, argModelParam) * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);

		_cumQALY += qol;
		if (argPat._cycles_on_HCC_treatment == 0) {
			_cumCost += argModelParam._costData.hccscr_c_tx_res / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
		}

		if (argPat._stateHCCWhenTreated == scr_hcc_small) {
			pMortLiv = argModelParam._transData._mort_tx_res_small;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_med) {
			pMortLiv = argModelParam._transData._mort_tx_res_med;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_large) {
			pMortLiv = argModelParam._transData._mort_tx_res_large;
		}
		else {
			ExitWithMsg("ERROR @ RUN_HCCTreatment(): Incorrect HCC state for Resection treatment = "+basicToStr((int)argPat._stateHCCWhenTreated));
		}	

	}
	// -------------------------------------------------------------
	// ------------ ablation treatment -----------------------------
	// -------------------------------------------------------------
	else if (argPat._stateAbs == scr_abs_abl) {
		double qol = argPat._cycles_on_HCC_treatment < _cycles_of_one_year ? argModelParam._qolData.qChange_tx_abl_1y : argModelParam._qolData.qChange_tx_abl_1yPlus;
		qol = qol * GetTxBaselineQoL(argPat, argModelParam) * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);

		_cumQALY += qol;
		if (argPat._cycles_on_HCC_treatment == 0) {
			_cumCost += argModelParam._costData.hccscr_c_tx_abl / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
		}

		//-- check fib state ---
		if (argPat._stateFibWhenTreated == scr_fib_inAbs) { ExitWithMsg("ERROR @ Run_HCCTReatment(): ABL treatment, _stateFibWhenTreated should not be scr_fib_inAbs..."); }

		if (argPat._stateHCCWhenTreated == scr_hcc_small && argPat._stateFibWhenTreated != scr_fib_DC) {
			pMortLiv = argModelParam._transData._mort_tx_abl_cc_small;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_med && argPat._stateFibWhenTreated != scr_fib_DC) {
			pMortLiv = argModelParam._transData._mort_tx_abl_cc_med;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_large && argPat._stateFibWhenTreated != scr_fib_DC) {
			pMortLiv = argModelParam._transData._mort_tx_abl_cc_large;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_small && argPat._stateFibWhenTreated == scr_fib_DC) {
			pMortLiv = argModelParam._transData._mort_tx_abl_dc_small;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_med && argPat._stateFibWhenTreated == scr_fib_DC) {
			pMortLiv = argModelParam._transData._mort_tx_abl_dc_med;
		}
		else if (argPat._stateHCCWhenTreated == scr_hcc_large && argPat._stateFibWhenTreated == scr_fib_DC) {
			pMortLiv = argModelParam._transData._mort_tx_abl_dc_large;
		}
		else {
			ExitWithMsg("ERROR @ RUN_HCCTreatment(): Incorrect state for Ablation treatment");
		}
	}
	// -------------------------------------------------------------
	// ------------ palliative treatment ---------------------------
	// -------------------------------------------------------------
	else if (argPat._stateAbs == scr_abs_pal) {		
		double qol = argPat._cycles_on_HCC_treatment < _cycles_of_one_year ? argModelParam._qolData.qChange_tx_pal_1y : argModelParam._qolData.qChange_tx_pal_1yPlus;
		qol = qol * GetTxBaselineQoL(argPat, argModelParam) * qolAge / _cycles_of_one_year / pow(1.0 + DISCOUNT_FACT_Q, (double)_sim_cycle / _cycles_of_one_year);

		_cumQALY += qol;
		if (argPat._cycles_on_HCC_treatment == 0) {
			_cumCost += argModelParam._costData.hccscr_c_tx_pal / pow(1.0 + DISCOUNT_FACT_C, (double)_sim_cycle / _cycles_of_one_year);
		}		
		pMortLiv = argModelParam._transData._mort_tx_pal;

	}
	else {
		ExitWithMsg("Error @ RUn_HCC_Treatment(): Incorrect treatment absorbing state = " + basicToStr((int)argPat._stateAbs));
	}

	// ================================  state transitions ======================================================================
	double rTrans = _rnd_pt.GetU01();
	if (rTrans < pMortBg) {
		argPat._stateAbs = scr_abs_death_bg;
	}
	else if (rTrans >= pMortBg && rTrans < pMortBg + pMortLiv) {
		argPat._stateAbs = scr_abs_death_liv;
	}
	else {
		// do nothing, no change
	}
	
	argPat._cycles_on_HCC_treatment++;

	return 0;


}


int HCCScreenSim::InitializePatient(HCCScrPatientType & pat, const modelParamType& argModelParam )
{
	pat._gender = _rnd_pt.GetU01() < argModelParam._hccScreeningData._distr_male ? 'M' : 'F';
	pat._genotype = DiscreteDistrSampler_DistrTable(argModelParam._hccScreeningData._distr_genotype, _rnd_pt);

	pat._initialAge = argModelParam._hccScreeningData._initial_age;
	pat._currentAge = pat._initialAge;
	pat._stateFib = DiscreteDistrSampler_DistrTable(argModelParam._hccScreeningData._distr_initial_fib, _rnd_pt);
	pat._stateHCC = scr_hcc_none;
	pat._stateAbs = scr_abs_none;
	pat._cycles_on_HCC_treatment = 0;
	return 0;
}

int HCCScreenSim::StateTransDeath(double pProgHCC, double pProgFib, double bgMort, double livMort, 
	type_scr_fib_state &  state_next_fib, type_scr_hcc_state & state_next_hcc, type_scr_absb_state & state_next_abs)
{
	double pNoProb = (1 - pProgHCC) * (1 - pProgFib);
	// death probability
	double pTrans = _rnd_pt.GetU01();
	if (pTrans < bgMort / pNoProb) {
		state_next_fib = scr_fib_inAbs;
		state_next_hcc = scr_hcc_inAbs;
		state_next_abs = scr_abs_death_bg;
	}
	else if ((pTrans >= bgMort / pNoProb) && (pTrans < (bgMort + livMort) / pNoProb)) {
		state_next_fib = scr_fib_inAbs;
		state_next_hcc = scr_hcc_inAbs;
		state_next_abs = scr_abs_death_liv;
	}
	else {
		// do nothing
	}
	return 0;
}


type_scr_hcc_state HCCScreenSim::StateTransHCC(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram)
{
	type_scr_hcc_state nextHCC;
	double pTrans = _rnd_pt.GetU01();
	if (argPat._stateHCC == scr_hcc_none) {
		double pInci = argModelPaaram._transData._map_inci_HCC.find(argPat._stateFib)->second;
		nextHCC = pTrans < pInci ? scr_hcc_small : scr_hcc_none;
	}
	else if (argPat._stateHCC == scr_hcc_small) {
		nextHCC = pTrans < argModelPaaram._transData._pr_monthly_hcc_small_med ? scr_hcc_med : scr_hcc_small;
	}
	else if (argPat._stateHCC == scr_hcc_med) {
		nextHCC = pTrans < argModelPaaram._transData._pr_monthly_hcc_med_large ? scr_hcc_large : scr_hcc_med;
	}
	else if (argPat._stateHCC == scr_hcc_large) {
		// do nothing, just stay in large
		nextHCC = scr_hcc_large;
	}
	else if (argPat._stateHCC == scr_hcc_inAbs) {
		ExitWithMsg("ERROR @ StateTransHCC(): _stateHCC should not be scr_hcc_inAbs in natural history simulation module...");
	}
	else {
		ExitWithMsg("Error # StateTransHCC: Unknown HCC state = " + basicToStr((int)argPat._stateHCC));
	}
	return nextHCC;

}

double HCCScreenSim::GetMortLiv(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram)
{
	if (argPat._stateHCC == scr_hcc_large) {
		return argModelPaaram._transData._mort_hcc_large;
	}
	else {
		if (argPat._stateFib == scr_fib_F3 || argPat._stateFib == scr_fib_F3_SVR || argPat._stateFib == scr_fib_F4_SVR) {
			return 0.0;
		}
		else if (argPat._stateFib == scr_fib_F4) {
			return argModelPaaram._transData._mort_F4;
		}
		else if (argPat._stateFib == scr_fib_DC) {
			return argModelPaaram._transData._mort_DC;
		}
		else {
			ExitWithMsg("Error @ GetMortLiv(): Unknown patient health state");
		}
	}
	return 0.0;
}

double HCCScreenSim::GetProbHCCProg(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram)
{
	if (argPat._stateHCC == scr_hcc_none) {
		return argModelPaaram._transData._map_inci_HCC.find(argPat._stateFib)->second;
	}
	else if (argPat._stateHCC == scr_hcc_small) {
		return argModelPaaram._transData._pr_monthly_hcc_small_med;
	}
	else if (argPat._stateHCC == scr_hcc_med) {
		return argModelPaaram._transData._pr_monthly_hcc_med_large;
	}
	else if (argPat._stateHCC == scr_hcc_large || argPat._stateHCC == scr_hcc_inAbs) {
		return 0.0;
	}
	else {
		ExitWithMsg("Error @ GetProbHCCProg: unknown patient state = " + basicToStr((int)argPat._stateHCC));
	}


}

double HCCScreenSim::GetTxBaselineQoL(const HCCScrPatientType & argPat, const modelParamType & argModelPaaram)
{
	if (argPat._stateFibWhenTreated == scr_fib_F3_SVR || argPat._stateFibWhenTreated == scr_fib_F4_SVR) {
		return argModelPaaram._qolData.q_SVR;
	}
	else if (argPat._stateFibWhenTreated == scr_fib_F3) {
		return argModelPaaram._qolData.q_F3;
	}
	else if (argPat._stateFibWhenTreated == scr_fib_F4) {
		return argModelPaaram._qolData.q_CoCirr;
	}
	else if (argPat._stateFibWhenTreated == scr_fib_DC) {
		return argModelPaaram._qolData.q_DeCirr;
	}
	else {
		ExitWithMsg("ERROR @ GetTXBaselineQoL(): incorrect baseline state!");
		return 0;
	}


}

string HCCScreenSim::GetFibHCCStateStr(const HCCScrPatientType & argPat)
{
	if ((argPat._stateFib == scr_fib_F3 || argPat._stateFib == scr_fib_F3_SVR) && argPat._stateHCC == scr_hcc_small) {
		return "f3_small";
	}
	else if ((argPat._stateFib == scr_fib_F3 || argPat._stateFib == scr_fib_F3_SVR) && argPat._stateHCC == scr_hcc_med) {
		return "f3_med";
	}
	else if ((argPat._stateFib == scr_fib_F3 || argPat._stateFib == scr_fib_F3_SVR) && argPat._stateHCC == scr_hcc_large) {
		return "f3_large";
	}
	else if ((argPat._stateFib == scr_fib_F4 || argPat._stateFib == scr_fib_F4_SVR) && argPat._stateHCC == scr_hcc_small) {
		return "cc_small";
	}
	else if ((argPat._stateFib == scr_fib_F4 || argPat._stateFib == scr_fib_F4_SVR) && argPat._stateHCC == scr_hcc_med) {
		return "cc_med";
	}
	else if ((argPat._stateFib == scr_fib_F4 || argPat._stateFib == scr_fib_F4_SVR) && argPat._stateHCC == scr_hcc_large) {
		return "cc_large";
	}
	else if (argPat._stateFib == scr_fib_DC && argPat._stateHCC == scr_hcc_small) {
		return "dc_small";
	}
	else if (argPat._stateFib == scr_fib_DC && argPat._stateHCC == scr_hcc_med) {
		return "dc_med";
	}
	else if (argPat._stateFib == scr_fib_DC && argPat._stateHCC == scr_hcc_large) {
		return "dc_large";
	}
	else {
		ExitWithMsg("ERROR @ GetFibHCCStateStr(): Incorrect fib/hcc state!");
		return 0;
	}
}

int HCCScreenSim::ConvertParameters(modelParamType & argModelParam)
{

	//convert annual transition probability according to cycle length probabilities	
	argModelParam._transData.pr_F0_F1 = funcConvertProbForCycle(argModelParam._transData.pr_F0_F1, _cycles_of_one_year);
	argModelParam._transData.pr_F1_F2 = funcConvertProbForCycle(argModelParam._transData.pr_F1_F2, _cycles_of_one_year);
	argModelParam._transData.pr_F2_F3 = funcConvertProbForCycle(argModelParam._transData.pr_F2_F3, _cycles_of_one_year);
	argModelParam._transData.pr_F3_CoCirr = funcConvertProbForCycle(argModelParam._transData.pr_F3_CoCirr, _cycles_of_one_year);
	argModelParam._transData.pr_F3_HCC = funcConvertProbForCycle(argModelParam._transData.pr_F3_HCC, _cycles_of_one_year);
	argModelParam._transData.pr_CoCirr_DeCirr = funcConvertProbForCycle(argModelParam._transData.pr_CoCirr_DeCirr, _cycles_of_one_year);
	argModelParam._transData.pr_CoCirr_HCC = funcConvertProbForCycle(argModelParam._transData.pr_CoCirr_HCC, _cycles_of_one_year);
	argModelParam._transData.pr_DeCirr_HCC = funcConvertProbForCycle(argModelParam._transData.pr_DeCirr_HCC, _cycles_of_one_year);
	argModelParam._transData.pr_DeCirr_LivTr = funcConvertProbForCycle(argModelParam._transData.pr_DeCirr_LivTr, _cycles_of_one_year);
	argModelParam._transData.pr_DeCirr_DeathLiv = funcConvertProbForCycle(argModelParam._transData.pr_DeCirr_DeathLiv, _cycles_of_one_year);
	argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = funcConvertProbForCycle(argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv, _cycles_of_one_year);
	argModelParam._transData.pr_HCC_LivTr = funcConvertProbForCycle(argModelParam._transData.pr_HCC_LivTr, _cycles_of_one_year);
	argModelParam._transData.pr_HCC_DeathLiv = funcConvertProbForCycle(argModelParam._transData.pr_HCC_DeathLiv, _cycles_of_one_year);
	argModelParam._transData.pr_LivTr_DeathLiv = funcConvertProbForCycle(argModelParam._transData.pr_LivTr_DeathLiv, _cycles_of_one_year);
	argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = funcConvertProbForCycle(argModelParam._transData.pr_LivTr1yrPlus_DeathLiv, _cycles_of_one_year);
	argModelParam._transData.pr_SVR_CoCirr_DeCirr = funcConvertProbForCycle(argModelParam._transData.pr_SVR_CoCirr_DeCirr, _cycles_of_one_year);
	argModelParam._transData.pr_SVR_CoCirr_HCC = funcConvertProbForCycle(argModelParam._transData.pr_SVR_CoCirr_HCC, _cycles_of_one_year);

	argModelParam._transData._map_inci_HCC[scr_fib_F3] = funcConvertProbForCycle(argModelParam._transData._map_inci_HCC[scr_fib_F3], _cycles_of_one_year);
	argModelParam._transData._map_inci_HCC[scr_fib_F4] = funcConvertProbForCycle(argModelParam._transData._map_inci_HCC[scr_fib_F4], _cycles_of_one_year);
	argModelParam._transData._map_inci_HCC[scr_fib_DC] = funcConvertProbForCycle(argModelParam._transData._map_inci_HCC[scr_fib_DC], _cycles_of_one_year);
	argModelParam._transData._map_inci_HCC[scr_fib_F3_SVR] = funcConvertProbForCycle(argModelParam._transData._map_inci_HCC[scr_fib_F3_SVR], _cycles_of_one_year);
	argModelParam._transData._map_inci_HCC[scr_fib_F4_SVR] = funcConvertProbForCycle(argModelParam._transData._map_inci_HCC[scr_fib_F4_SVR], _cycles_of_one_year);

	argModelParam._transData._mort_hcc_large = funcConvertProbForCycle(argModelParam._transData._mort_hcc_large, _cycles_of_one_year);
	argModelParam._transData._mort_F4 = funcConvertProbForCycle(argModelParam._transData._mort_F4, _cycles_of_one_year);
	argModelParam._transData._mort_DC = funcConvertProbForCycle(argModelParam._transData._mort_DC, _cycles_of_one_year);

	argModelParam._transData._mort_tx_transpl_1y = funcConvertProbForCycle(argModelParam._transData._mort_tx_transpl_1y, _cycles_of_one_year);
	argModelParam._transData._mort_tx_transpl_1yPlus = funcConvertProbForCycle(argModelParam._transData._mort_tx_transpl_1yPlus, _cycles_of_one_year);
	argModelParam._transData._mort_tx_res_small = funcConvertProbForCycle(argModelParam._transData._mort_tx_res_small, _cycles_of_one_year);
	argModelParam._transData._mort_tx_res_med = funcConvertProbForCycle(argModelParam._transData._mort_tx_res_med, _cycles_of_one_year);
	argModelParam._transData._mort_tx_res_large = funcConvertProbForCycle(argModelParam._transData._mort_tx_res_large, _cycles_of_one_year);
	argModelParam._transData._mort_tx_abl_cc_small = funcConvertProbForCycle(argModelParam._transData._mort_tx_abl_cc_small, _cycles_of_one_year);
	argModelParam._transData._mort_tx_abl_cc_med = funcConvertProbForCycle(argModelParam._transData._mort_tx_abl_cc_med, _cycles_of_one_year);
	argModelParam._transData._mort_tx_abl_cc_large = funcConvertProbForCycle(argModelParam._transData._mort_tx_abl_cc_large, _cycles_of_one_year);
	argModelParam._transData._mort_tx_abl_dc_small = funcConvertProbForCycle(argModelParam._transData._mort_tx_abl_dc_small, _cycles_of_one_year);
	argModelParam._transData._mort_tx_abl_dc_med = funcConvertProbForCycle(argModelParam._transData._mort_tx_abl_dc_med, _cycles_of_one_year);
	argModelParam._transData._mort_tx_abl_dc_large = funcConvertProbForCycle(argModelParam._transData._mort_tx_abl_dc_large, _cycles_of_one_year);
	argModelParam._transData._mort_tx_pal = funcConvertProbForCycle(argModelParam._transData._mort_tx_pal, _cycles_of_one_year);
	



	//convert costs according to cycle length  	

	argModelParam._costData.c_F0 = argModelParam._costData.c_F0 / _cycles_of_one_year;
	argModelParam._costData.c_F1 = argModelParam._costData.c_F1 / _cycles_of_one_year;
	argModelParam._costData.c_F2 = argModelParam._costData.c_F2 / _cycles_of_one_year;
	argModelParam._costData.c_F3 = argModelParam._costData.c_F3 / _cycles_of_one_year;
	argModelParam._costData.c_CoCirr = argModelParam._costData.c_CoCirr / _cycles_of_one_year;
	argModelParam._costData.c_DeCirr = argModelParam._costData.c_DeCirr / _cycles_of_one_year;
	argModelParam._costData.c_DeCirr1yrPlus = argModelParam._costData.c_DeCirr1yrPlus / _cycles_of_one_year;
	argModelParam._costData.c_HCC = argModelParam._costData.c_HCC / _cycles_of_one_year;
	argModelParam._costData.c_LivTr = argModelParam._costData.c_LivTr; // this is one-time lump-sum cost of liver transpalnt
	argModelParam._costData.c_PostLivTr = argModelParam._costData.c_PostLivTr / _cycles_of_one_year;

	// scale the transplant cost; other treatment cost is one-time cost
	argModelParam._costData.hccscr_c_tx_transpl_1y = argModelParam._costData.hccscr_c_tx_transpl_1y / _cycles_of_one_year;
	argModelParam._costData.hccscr_c_tx_transpl_1yPlus = argModelParam._costData.hccscr_c_tx_transpl_1yPlus / _cycles_of_one_year;




	// added 2017/03/09
	argModelParam._costData.c_background = argModelParam._costData.c_background / _cycles_of_one_year;

	return 0;

}

int HCCScreenSim::ReadSmpParamForPSA(ifstream & inf, modelParamType & argModelParam)
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

int HCCScreenSim::OutputBaseCase(ofstream & outf)
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

int HCCScreenSim::PSA_initializeSampler(const vector<string> & vecVarName, const vector<string> & vecDistrName, const vector<double> & vecP1, const vector<double> & vecP2)
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

int HCCScreenSim::PSA_initializeSampler(const psaDistrType & argPSADistr)
{
	return PSA_initializeSampler(argPSADistr._distrPSA_varName, argPSADistr._distrPSA_distrName, argPSADistr._distrPSA_param1, argPSADistr._distrPSA_param2);
}

int HCCScreenSim::PSA_sampleModelParamValue(modelParamType & argModelParam)
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


int HCCScreenSim::PSA_sampleModelParamValue_Univariate(modelParamType & argModelParam, string argVarName)
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

int HCCScreenSim::PSA_sampleCohort(baseCohortType & argCohort)
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

int HCCScreenSim::ResetCounters()
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

