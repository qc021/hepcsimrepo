#ifndef __DATA_TYPE_H__
#define __DATA_TYPE_H__
#include <iostream>
#include "sim_config.h"
#include "disburdn_config.h"
#include "SysUtil.h"

using namespace std;





const int MAXCYCLE = 5000;
double funcConvertProb(double transProb); //convert annual transition probability according to cycle length
double funcConvertRate(double transRate); //convert annual rate according to cycle length probability
double funcConvertProbForCycle(double transProb, double numCycle);

// =======================================
// enum types for Markov simulation
// --------------------------------------
// this corresponds to state = {0,1,2,...}; No need to remember numbers corresponding to each state.
enum txHistryType {TN, TE};
enum txHistryTolType {TN_TOL, TN_INT, TE_NA};
//enum stateType {s_F0, s_F1, s_F2, s_F3, s_CoCirr, s_DeCirr, s_DeCirr1yrPlus, s_HCC, s_LivTr, s_LivTr1yrPlus, s_SVR, s_Death}; 
enum stateType { s_F0, s_F1, s_F2, s_F3, s_CoCirr, s_DeCirr, s_DeCirr1yrPlus, s_HCC, s_LivTr, s_LivTr1yrPlus, s_SVR, s_Death, s_Acute }; // modified 2016/11/28 by Qiushi: incorporate accute state
enum drugType {PEG, RBV, BOC, SOF, TEL, SMV, LDP, DCV, PrOD, ASV};

// =======================================
// enum types for disease burden model
// ---------------------------------------
enum typeScreeningScenario {screen_no_extra, screen_birthCohort, screen_universal_rate, screen_universal_capacity};

const string typeStateDisBurdnModel_label[] = { "F0","F1","F2","F3","F4","DC","DC+","HCC","LivTr","LivTr+",
											"Death","DLR","curedF0F2","curedF3","curedF4","curedAcute","NA"};
enum typeStateDisBurdnModel {state_F0, state_F1, state_F2, state_F3, state_CoCirr, state_DeCirr, state_DeCirr1yrPlus, state_HCC, state_LivTr, state_LivTr1yrPlus,	//states in the model
	state_Death, state_DeathLR, state_curedF0F2, state_curedF3, state_curedCoCirr, state_curedAcute, state_NA};
	//state_curedFib_G1, state_curedCoCirr_G1, state_curedAcute, state_cured_Contra, state_curedFib_G2, state_curedFib_G3, 
	//state_curedFib_G456, state_curedCoCirr_G2, state_curedCoCirr_G3, state_curedCoCirr_G456, state_NA}; 

enum typeCohortDisBurdnModel {cohort_defaultSingleCohort, cohort_incarcerated, cohort_homeless, cohort_hospitalized, cohort_nursingHome, cohort_indianReservation};
const string typeCohort_label[] = { "cohort_default","cohort_incarcerated","cohort_homeless", "cohort_hospitalized", "cohort_nursingHome","cohort_indianRsv" };


enum typeTxResponseStatus {
	trStatus_naive, trStatus_svr,
	trStatus_relap, trStatus_partial, trStatus_null, //treatment response/history status
	trStatus_contraInd_mod, trStatus_contraInd_nonmod,
	trStatus_failed_PI1,
	trStatus_failed_DAA1_nonNS5A, 
	trStatus_failed_DAA2_nonNS5A, trStatus_failed_DAA2_NS5A,
	trStatus_failed_DAA3_nonNS5A, trStatus_failed_DAA3_NS5A,
	trStatus_unknown};

enum typeTxCategory {txCat_PEGRBV , txCat_PI1, txCat_DAA1_nonNS5A, txCat_DAA2_nonNS5A, txCat_DAA2_NS5A,  txCat_DAA3_nonNS5A, txCat_DAA3_NS5A, txCat_NoTx};	//treatment types 
//PI1: triple therapy
//wave1_G1: future DAAs for genptype 1 (2014)
//wave2_G1: future DAAs for genptype 1 higher efficacy for cirrhotics (2017)
//wave3_G1: future therapies higher efficacy for cirrhotic partial and null responders (2020)
//wave1_G23: future DAAs for genptype 2/3 all response states
//prTr_NA: no treatment considered 


//enum typeTxNS5A {tx_nonNS5A, tx_NS5A, tx_NA_NS5A};

enum typeTxPriority {priority_no, priority_PrioritizeF3F4, priority_OnlyF3F4Until2030, pirority_PrioritizeF3F4Until2030}; // default=pirority_PrioritizeF3F4Until2030


const string typeInsurance_label[] = { "insr_uninsured", "insr_private", "insr_medicare", "insr_medicaid", "insr_military", "insr_incarcerated", "insr_indianRsv" };
enum typeInsurance {insr_uninsured, insr_private, insr_medicare, insr_medicaid, insr_military, insr_incarcerated, insr_indianRsv};

enum typeBirthCohortGroup {bc_group_general, bc_group_indianRsv, bc_group_noBirthCohortScr};
enum typeTxCohortGroup { tx_group_general, tx_group_incarcerated, tx_group_indianRsv, tx_group_noTx };

// [genotype, treatment wave, response status of previous treatment, fibrosis state ] -> SVR value
typedef map<int, map<typeTxCategory,map<typeTxResponseStatus,map<typeStateDisBurdnModel, double>>>> map_svr_table;



// added by Qiushi @ 6/19/2017 for HCC screening analysis
enum type_scr_fib_state { scr_fib_F3, scr_fib_F3_SVR, scr_fib_F4, scr_fib_F4_SVR, scr_fib_DC, scr_fib_inAbs };
enum type_scr_hcc_state { scr_hcc_none, scr_hcc_small, scr_hcc_med, scr_hcc_large, scr_hcc_inAbs };
enum type_scr_absb_state { scr_abs_res, scr_abs_abl, scr_abs_trspl, scr_abs_pal, scr_abs_death_bg, scr_abs_death_liv, scr_abs_none };



// -----------------------------------------

class baseCohortType{	// used as cohort profile
public: 
	baseCohortType(){};	// [QC] only used to hold the space in modelParamType constructor()
	baseCohortType(int argState, double argAge, char argGender, char argRace, 
		int argGeno, char argIL28B, char argPriorTxRes ): baseState(argState), 
		baseAge(argAge),
		baseGender(argGender),
		baseRace(argRace),
		baseGenotype(argGeno),
		baseIL28B(argIL28B),
		basePriorTxRes(argPriorTxRes) {};

	~baseCohortType(){};
	//[QC]: 7 attributes
	int baseState;
	double baseAge;
	char baseGender;
	char baseRace;
	int baseGenotype;
	char baseIL28B; //C (CC), X (CT), T (TT), N(none)
	char basePriorTxRes;
};

class disBurdnCounterType {
public:
	disBurdnCounterType() {
		_year_start_lrd_ageDistr = 2015;
		_record_ageDistr_year_lrd = 2020;
	};
	~disBurdnCounterType() {};

	map<double, int> _counter_ageDistr_xsectional_given_year;
	vector<int> _counterPpl, _counterHCV, _counterHCV_complement;

	map<typeCohortDisBurdnModel, vector<int>> _counter_byCohort_HCV;
	map<typeCohortDisBurdnModel, vector<int>> _counter_byCohort_awareness;
	map<typeCohortDisBurdnModel, vector<int>> _counter_byCohort_screening;
	map<typeCohortDisBurdnModel, vector<int>> _counter_byCohort_treatment;

	vector<int> _counter_new_hcv_incidence;
	vector<int> _counter_liverRelatedDeath;
	vector<int> _counter_SVR, _counter_aware, _counter_unaware, _counter_aware_F0F4;
	map<typeStateDisBurdnModel, vector<int>> _counter_aware_byFib;
	vector<int> _counter_failedSVR;
	vector<int> _counter_failedSVR_retreatable_everFailedNS5A, _counter_failedSVR_retreatable_neverFailedNS5A; // added 04/06/2017
	vector<int> _counter_failedSVR_unretreatable;
	vector<int> _counter_failedSVR_everReceivedNS5A, _counter_failedSVR_neverReceivedNS5A; // added 04/06/2017
	vector<int> _counter_treated_alive;

	vector<int> _counter_incidence_DC, _counter_incidence_HCC;

	vector<int> _counter_transplant_DC, _counter_transplant_HCC;

	map<typeStateDisBurdnModel, vector<int>> _counterHealthState;
	map<typeInsurance, vector<int>> _counterInsurance;

	vector<int> _counterTxFailure_numTx;
	map<string, vector<int>> _counterTxFailure_numTx_byNS5A;
	map<string, vector<int>> _counterTxFailure_numTx_byCirr;
	map<int, vector<int>> _counterTxFailure_numTx_byGenotype;

	vector<int> _counterTxFailure_numFailure;
	map<int, vector<int>> _counterTxFailure_numFailure_byGenotype;
	map<string, vector<int>> _counterTxFailure_numFailure_byNS5A;
	map<string, vector<int>> _counterTxFailure_numFailure_byCirr;
	map<string, vector<int>> _counterTxFailure_numFailure_byCirr_NS5A;

	vector<int> _counterTxFailure_numRetxCandidates;

	vector<int> _counter_screening;
	vector<int> _counter_txElig_all, _counter_txElig_F3F4;


	// baby boomers
	vector<int> _counter_birthcohort_infected, _counter_birthcohort_aware, _counter_birthcohort_unaware, 
		_counter_birthcohort_becomingAware, _counter_birthcohort_cured, _counter_birthcohort_gettingTx,
		_counter_birthcohort_infected_insured, _counter_birthcohort_unaware_insured;


	// DALYs
	vector<double> _counter_DALY_YLL, _counter_DALY_YLD, _counter_DALY_Total;
	

	// ------------------------
	int _year_start_lrd_ageDistr, _record_ageDistr_year_lrd;
	map<double,int> _counter_LRD_ageDistr;

	map<double, int> _counter_ageDistr;

	// ------- cost counters -----------
	vector<double> _counterCost;
	map<string, vector<double>> _counterCost_byCategory;
	map<typeInsurance, vector<double>> _counterCost_byInsr;
	map<string, vector<double>> _counterCost_byLivDis;


};
class disBurdnDataType
{
public:
	disBurdnDataType() {
		_flag_use_constant_transplant_capacity = true;
		_flag_multipleCohort = false;		
		_flag_include_immigrants_lpr = false; // lpr = lawful permanent resident = green card holders
		_flag_extended_use_PR_PI = false;
		_num_initial_population = DISBDEN_INITIAL_NUM_PATIENTS;
		_option_screenScenario = screen_birthCohort;

		// initialize values
		_prob_male = 0.6422;
		_prob_receivable = 0.915;
		_prob_TrR = 1.0; // ???
		_prob_chronic = 0.78;
		_prob_tx_coverage_by_medicarePartD =  0.9;
		_prob_naive_initialization = 1.0; // 0.65, or 0.17???
		_prob_chronic_contra = 0.346;
		_prob_mod_given_contra = 0.231 / 0.346;
		_prob_screen_accept = 0.819;		//.91 * .90, Rein 2012. (73.71-90.09)
											//percentage of people who accept to be screened for HepC
											//Rein 2012 assumed that 28%of the chronic cases are aware of their infection, so screening would not benefit them.
		_prob_birthCohortScr_insured = 0.9;
		_prob_birthCohortScr_uninsured = 0.1;
		_prob_delay_tx_until_wave1_f0f2 = 0.75;
		_prob_delay_tx_until_wave1_f3 = 0.25;

		_prob_IFNintol = 0.17;


		// ------ modified on 7/16/2017 -------------------
		// Notes: We now spearate the awareness rate between [initial population] and [new incidences]. 
		// Assign the default values (used in previous manuscripts) for [initial population]
		// Update the awareness rate for [new incidences] based on CDC's estimates = 7% (1 in 13)

		_table_aware_prob_initial_population_uninsured.clear();
		_table_aware_prob_initial_population_uninsured[40] = 0.0581;// 0.42 * 0.1383; // lower_bounds(): < 40
		_table_aware_prob_initial_population_uninsured[50] = 0.1950;// 0.71*0.2747;	// 40-49
		_table_aware_prob_initial_population_uninsured[60] = 0.1629;// 0.61*0.2670;	// 50-59
		_table_aware_prob_initial_population_uninsured[200] = 0.0577;// 0.27*0.2136; // 60-100

		_table_aware_prob_initial_population_insured.clear();
		_table_aware_prob_initial_population_insured[40] = 0.1763;// 0.53*0.3326;
		_table_aware_prob_initial_population_insured[50] = 0.6045;// 0.915*0.6606;
		_table_aware_prob_initial_population_insured[60] = 0.5459;// 0.85*0.6423;
		_table_aware_prob_initial_population_insured[200] = 0.3494;// 0.68*0.5138;
				
		_table_aware_prob_among_new_incidence_uninsured = _table_aware_prob_initial_population_uninsured;
		_table_aware_prob_among_new_incidence_insured = _table_aware_prob_initial_population_insured;

		_flag_override_awareness_for_initial_population = false;
		_awareness_for_initial_population_if_overrided = 0;


		// <genotype-int, probability>
		_table_distr_genotype[1] = 0.73;
		_table_distr_genotype[2] = 0.14;
		_table_distr_genotype[3] = 0.08;
		_table_distr_genotype[456] = 0.05;



		// <fibrosis, probability>
		//  transData.pr_cycle0_CoCirr = atof(argv[28]);		0.0554	0.0554
		//	transData.pr_cycle0_DeCirr = atof(argv[29]);		0.05746	0.00206
		//	transData.pr_cycle0_DeCirr1yrPlus = atof(argv[30]);		0.06674	0.00928
		//	transData.pr_cycle0_HCC = atof(argv[31]);		0.06674	0
		//	transData.pr_cycle0_F0 = atof(argv[32]);		0.51168	0.44494
		//	transData.pr_cycle0_F1 = atof(argv[33]);		0.77642	0.26474
		//	transData.pr_cycle0_F2 = atof(argv[34]);		0.91212	0.1357
		_table_distr_fib[state_CoCirr] = 0.0554;
		_table_distr_fib[state_DeCirr] = 0.00206;
		_table_distr_fib[state_DeCirr1yrPlus] = 0.00928;
		_table_distr_fib[state_HCC] = 0;
		_table_distr_fib[state_F0] = 0.44494;
		_table_distr_fib[state_F1] = 0.26474;
		_table_distr_fib[state_F2] = 0.1357;
		_table_distr_fib[state_F3] = 1 - 0.91212;

		_table_distr_insurance[insr_uninsured] = 0.267;
		_table_distr_insurance[insr_private] = 0.498; // 0.7649 - 0.2666;
		_table_distr_insurance[insr_medicaid] = 0.143; // 0.9078 - 0.7649;
		_table_distr_insurance[insr_military] = 0.092; // 1.0 - 0.9078;

		ACATotalratio[2014] = 0.222;
		ACAPrivateratio[2014] = 0.103;
		ACAMedicaidratio[2014] = 0.120;
		ACAPrivatechange[2014] = 0.103;
		ACAMedicaidchange[2014] = 0.120;

		ACATotalratio[2015] = 0.352;
		ACAPrivateratio[2015] = 0.173;
		ACAMedicaidratio[2015] = 0.179;
		ACAPrivatechange[2015] = 0.070;
		ACAMedicaidchange[2015] = 0.059;

		ACATotalratio[2016] = 0.463;
		ACAPrivateratio[2016] = 0.247;
		ACAMedicaidratio[2016] = 0.216;
		ACAPrivatechange[2016] = 0.074;
		ACAMedicaidchange[2016] = 0.037;

		ACATotalratio[2017] = 0.481;
		ACAPrivateratio[2017] = 0.259;
		ACAMedicaidratio[2017] = 0.222;
		ACAPrivatechange[2017] = 0.013;
		ACAMedicaidchange[2017] = 0.006;


		ACATotalratio[2018] = 0.481;
		ACAPrivateratio[2018] = 0.259;
		ACAMedicaidratio[2018] = 0.222;
		ACAPrivatechange[2018] = 0.0;
		ACAMedicaidchange[2018] = 0.0;

		ACATotalratio[2100] = 0.481;
		ACAPrivateratio[2100] = 0.259;
		ACAMedicaidratio[2100] = 0.222;
		ACAPrivatechange[2100] = 0.0;
		ACAMedicaidchange[2100] = 0.0;

		_table_respStatus_distr_naive[trStatus_contraInd_mod] = _prob_mod_given_contra * _prob_chronic_contra;
		_table_respStatus_distr_naive[trStatus_contraInd_nonmod] = (1-_prob_mod_given_contra) * _prob_chronic_contra;
		_table_respStatus_distr_naive[trStatus_naive] = 1 -  _prob_chronic_contra;

		_table_respStatus_distr_expr_G1[trStatus_relap] = 0.53;
		_table_respStatus_distr_expr_G1[trStatus_partial] = 0.19;
		_table_respStatus_distr_expr_G1[trStatus_null] = 1-0.53-0.19;

		_table_respStatus_distr_expr_G234[trStatus_relap] = 0.47;
		_table_respStatus_distr_expr_G234[trStatus_partial] = 0.16;
		_table_respStatus_distr_expr_G234[trStatus_null] = 1-0.47-0.16;



		_table_universal_scr_rate[-1] = 0;
		_table_universal_scr_rate[0] = 0.142857;
		_table_universal_scr_rate[1] = 0.16667;
		_table_universal_scr_rate[2] = 0.2;
		_table_universal_scr_rate[3] = 0.25;
		_table_universal_scr_rate[4] = 0.333;
		_table_universal_scr_rate[5] = 0.5;
		_table_universal_scr_rate[6] = 1;
		_table_universal_scr_rate[7] = 0;
		_table_universal_scr_rate[100] = 0;


		//_UNINSURED_RATIO = 0.2666;  //insurance distribution among 65-
		//_PRIVATE_RATIO = 0.7649 - 0.2666;	// 0.7649;	 // this is a cumulative probability
		//_MEDI_RATIO = 0.9078 - 0.7649;		// 0.9078;   //	this is a cumulative probability

		// added 4/23/2017
		_flag_inci_proportional_to_prevalence = false;
		_year_start_constant_new_incidence = START_YR;


		_record_ageDistr_year = 2010;
		_record_ageDistr_outputFolder = ".";
	};
	~disBurdnDataType() {};

	bool _flag_multipleCohort; // = true if considering multiple non-NHANES population
	bool _flag_include_immigrants_lpr; // = true if considering immigrants


	long	_num_initial_population;
	typeScreeningScenario _option_screenScenario;	
	double _prob_male;
	double _prob_receivable;
	double _prob_TrR;
	double _prob_chronic;
	double _prob_tx_coverage_by_medicarePartD;
	double _prob_naive_initialization;
	double _prob_chronic_contra, _prob_mod_given_contra;
	double _prob_screen_accept;
	double _prob_birthCohortScr_insured;
	double _prob_birthCohortScr_uninsured;
	double _prob_delay_tx_until_wave1_f0f2;
	double _prob_delay_tx_until_wave1_f3;

	double _prob_IFNintol;
	map_svr_table _table_SVR, _table_SVR_old;

	// <age, value>
	map<int, double> _bgMort_male, _bgMort_female;
	map<int, double> _qol_male, _qol_female;

	// <age, value> added for AmericanIndian population
	map<int, double> _bgMort_AI_male, _bgMort_AI_female;
	map<int, double> _qol_AI_male, _qol_AI_female;


	// <year, screening_ratio_value>
	map<int, double> _table_scr_ratio;
	map<int, int> _table_tx_capacity;

	bool _flag_use_constant_transplant_capacity;
	map<int, int> _table_transplant_capacity;

	// <year, incidence>
	map<int, long> _table_incidence;
	// added for proportional incidence wrt prevalence (for Pakistan analysis) @ 4/27/2017
	bool _flag_inci_proportional_to_prevalence;
	map<int, double> _table_incidence_proportion;
	int _year_start_constant_new_incidence;

	// <year, <age distribution>>
	map<int, vector<double>> _table_inci_age_distr;
	vector<int> _vec_age_category_forInitialization, _vec_age_category_forNewIncidence;
	int _record_ageDistr_year;
	string _record_ageDistr_outputFolder;

	// <response, double>
	map<typeTxResponseStatus, double> _table_respStatus_distr_naive;
	map<typeTxResponseStatus, double> _table_respStatus_distr_expr_G1;
	map<typeTxResponseStatus, double> _table_respStatus_distr_expr_G234;

	
	map<int, map<int, double>> _table_ns5a_marketshare; // <year, <genotype, market share value>>
	map<int, map<int, double>> _table_pr_extended_use_pr_pi; //[year][genotype]->prob
	map<int, map<int, double>> _table_pr_extended_use_pr; //[year][genotype]->prob
	bool _flag_extended_use_PR_PI;

	// <genotype-int, probability>
	map<int, double> _table_distr_genotype;

	// <fibrosis, probability>
	map<typeStateDisBurdnModel, double> _table_distr_fib;

	// --- awareness ----
	map<double, double> _table_aware_prob_initial_population_uninsured;
	map<double, double> _table_aware_prob_initial_population_insured;

	map<double, double> _table_aware_prob_among_new_incidence_uninsured;
	map<double, double> _table_aware_prob_among_new_incidence_insured;
	double _flag_override_awareness_for_initial_population;
	double _awareness_for_initial_population_if_overrided;

	map<typeStateDisBurdnModel, map<bool, map<double, double>>> _table_prob_becoming_aware_by_fib_insur_age;

	// ---- insurance ---- 
	map<typeInsurance, double> _table_distr_insurance;
	map<int, map<typeInsurance, double>> _table_distr_insurance_ACA_by_year;

	map<int,double> ACATotalratio, 	ACAPrivateratio, 	ACAMedicaidratio, ACAPrivatechange, ACAMedicaidchange;

	//double _UNINSURED_RATIO;  //insurance distribution among 65-
	//double _PRIVATE_RATIO;	 // this is a cumulative probability
	//double _MEDI_RATIO;        //this is a cumulative probability

	// --- screening parameters ---
	map<int, double> _table_universal_scr_rate, _table_universal_scr_cap;

};


class patientType	// [QC] used in simulations (simulated patients)
{
public:
	patientType();	// assign default value in the constructor
	~patientType(){};
		
	typeCohortDisBurdnModel _cohort;
	typeBirthCohortGroup _bcScr_group;
	typeTxCohortGroup _tx_group;

	double currentAge, initialAge;
	char gender; // M/F
	char race; // W/B/A/M/P/N
	int genotype; //values = {1,2,3}
	char IL28B; //C (CC), X (CT), T (TT), N(none)
	char priorTxRes; // R (relapser), P (partial responder), N (null-responder), A (all)
	int state, initialState; // using enum, state can take these values {s_F0, s_F1, s_F2, s_F3, s_CoCirr, s_DeCirr, s_HCC, s_LivTr, s_LivTr1yrPlus, s_SVR, s_Death};
	
	bool flagDeCirr; //flag if the person has a history of decompensated cirrhosis
	bool flagHCC; //flag for history of hepatocellular carcinoma
	bool flagLivTr; // flag for liver transplant
	bool flagTxEx; // flag for previous treatment experienced
	bool flagETR; // flag if patient acheives ETR (undetectable HCV RNA at the end of treatment)
	bool flagCured; // flag if patient is cured with the current therapy
	bool flagNR; // flag if non-responder to drug (Bocepravir or MK-7009)
	bool flagReL; // flag for whether or not a person attains ETR and then relapses
	bool flagTxComplete; // flag if person complies to the assigned duration of treatment		
	bool flagDeathLiv; //flag if death was due to liver-related complications
	bool flagAEAnemia; //flag patients who get anemia 
	int numCycleLivTr; //number of cycles spent in liver transplant state
	int numCycleDeCirr; //number of cycles spent in decompensated cirrhosis state

	bool flagEverTreated;


	// =======================================
	// definitions for disease burden model
	// ---------------------------------------
	long long _ptID;
	int _yrBirth;
	int _yrMostRecentTx;
	int _yrMostRecentScr;
	int _yrDiagnosedByUsualCare;


	bool _flagReceivable;
	bool _flagTrR;
	bool _flagBirthCohortForScr;
	
	typeInsurance _insurStatus;
	bool _flagTrInsur;

	bool _flagAware;
	bool _flagAwareByScreening; // true if he is detected by birth-cohort/universal screening
	bool _flagAwareByUsualCare;
	bool _flagCured;


	bool _flagInterferonTol;
	
	bool _flagTxEligible;
	bool _flagDelayTx;
	typeStateDisBurdnModel _state, _stateAtTx, _state_previous;
	typeTxCategory _txCat;
	typeTxResponseStatus _curTxResps;

	int _cycleOnDrug_PEGRBV, _cycleOnDrug_PI, _cycleOnDrug_DAA;
	int _cycleSVRF4;
	
	bool _flagEverFailedNS5A, _flagEverReceivedNS5A; // added @ 04/07/2017 for US salvage analysis


	// =======================================
	// TFA: Added for Treatment Failure Analysis @ 5/23/2016, Qiushi CHen =================
	
	bool _tfaflag_eligible_for_tx;
	int _tfaflag_num_fail_pegrbv;
	int _tfaflag_num_fail_pi;
	typeTxCategory _tfaflag_current_tx_type;
	typeTxCategory _tfaflag_previous_failed_tx_type;
	int _tfaflag_total_tx_times;

	bool _tfaflag_ever_failed_nonns5a_before_2015;
	bool _tfaflag_ever_failed_ns5a;
	bool _tfaflag_ever_failed_nonns5a_after_2015;
	// ================================== END TFA modifications ==============================================================

	

};

// transition probabilities, mainly for natural history
class transitionType{
public:
	transitionType();
	~transitionType(){};
	double funcFibrosisProgression (double FibAge, char FibGender); // Assign progression rate based on age and gender
	double funcFibrosisF0F1 (double FibAge, char FibGender, int FibGenotype); // Assign progression rates based on risk equations from Thein et al.
	double funcFibrosisF1F2 (double FibAge, char FibGender, int FibGenotype); // Assign progression rates based on risk equations from Thein et al.
	double funcFibrosisF2F3 (double FibAge, char FibGender, int FibGenotype); // Assign progression rates based on risk equations from Thein et al.
	double funcFibrosisF3F4 (double FibAge, char FibGender, int FibGenotype); // Assign progression rates based on risk equations from Thein et al.



	double	pr_FibProg_M50under; //Fibrosis progression rate for male age under 50 (reference: Salomon JAMA 2003)
	double	pr_FibProg_M5059; //Fibrosis progression rate for male age 50-59 (reference: Salomon JAMA 2003)
	double	pr_FibProg_M6069; //Fibrosis progression rate for male age 60-69 (reference: Salomon JAMA 2003)
	double	pr_FibProg_M70plus; //Fibrosis progression rate for male age 70 plus (reference: Salomon JAMA 2003)
	double	pr_FibProg_F50under; //Fibrosis progression rate for female age under 50 (reference: Salomon JAMA 2003)
	double	pr_FibProg_F5059; //Fibrosis progression rate for female age 50-59 (reference: Salomon JAMA 2003)
	double	pr_FibProg_F6069; //Fibrosis progression rate for female age 60-69 (reference: Salomon JAMA 2003)
	double	pr_FibProg_F7079; //Fibrosis progression rate for female age 70-79 (reference: Salomon JAMA 2003)
	double	pr_FibProg_F80plus; //Fibrosis progression rate for female age 80 plus (reference: Salomon JAMA 2003)
	double	pr_F0_F1;
	double	pr_F1_F2;
	double	pr_F2_F3;
	double	pr_F3_CoCirr;
	double	pr_F3_HCC;
	double	pr_F3SVR_HCC; // added 1/10/2017
	double	pr_regression; // added 1/10/2017
	double	pr_CoCirr_DeCirr;
	double	pr_CoCirr_HCC;
	double	pr_SVR_CoCirr_DeCirr;
	double	pr_SVR_CoCirr_HCC;
	double	pr_DeCirr_HCC;
	double	pr_DeCirr_LivTr;
	double	pr_DeCirr_DeathLiv;
	double	pr_DeCirr1yrPlus_DeathLiv;
	double	pr_HCC_DeathLiv ; //source: Sibert et al. 2009
	double	pr_HCC_LivTr;
	double	pr_LivTr_DeathLiv; // SA range: (6%-42%) probability of death in first year of liver transplant; Sibert et al. 2009
	double	pr_LivTr1yrPlus_DeathLiv;// SA range (2.4%-11%) probability of death after 1st year of liver transplant; Sibert et al. 2009

	double	pr_SVR_Delta_DAA; // decrement in SVR rates used for sensitivity analysis
	double	pr_SVR_Delta_oSOC; // decrement in SVR rates used for sensitivity analysis

	double	pr_SVR_BOC; // probability of achieving SVR (using data from Phase II trial) 0.748 for 48week P/R/B with lead-in and 0.375 for control arm
	double	pr_ETR_BOC;// probability of achieving ETR (using data from Phase II trial) 0.786 for 48week P/R/B with lead-in and 0.51 for control arm
	double	pr_SVR_SOC; // probability of achieving SVR (using data from Phase II trial) 0.748 for 48week P/R/B with lead-in and 0.375 for control arm
	double	pr_ETR_SOC;// probability of achieving ETR (using data from Phase II trial) 0.786 for 48week P/R/B with lead-in and 0.51 for control arm
	double	pr_ContinueTx_TxNaive_BOC; //probability of continuing treatment in a given week
	double	pr_ContinueTx_TxNaive_RGT; //probability of continuing treatment in a given week
	double	pr_ContinueTx_TxNaive_SOC; //probability of continuing treatment in a given week

	double	pr_SVR_TxNaive_RGT28;
	double	pr_ETR_TxNaive_RGT28;
	double	pr_SVR_TxNaive_RGT48;
	double	pr_ETR_TxNaive_RGT48;

	double	pr_TxFailTW24_SOC; // probability of detectable HCV-RNA at TW 24 (treatment failure) in SOC
	double	pr_TxFailTW24_PRB48; // probability of detectable HCV-RNA at TW 24 (treatment failure) in PRB48
	double	pr_TxFailTW24_RGT;
	double	pr_ResponseTW24_TxNaive; //probability of response between TW8-TW24

	double	pr_SVR_TxExp_PRB48; //  107.0/161.0  (Phase III)
	double	pr_ETR_TxExp_PRB48;// 124.0/161.0 (Phase III)
	double	pr_SVR_TxExp_SOC; // 17.0/80.0 (Phase III)
	double	pr_ETR_TxExp_SOC;// 25.0/161.0 (Phase III)
	double	pr_SVR_TxExp_RGT36;// = 64.0/74.0;// (Phase III)
	double	pr_ETR_TxExp_RGT36;// = 72.0/74.0;// (Phase III)
	double	pr_SVR_TxExp_RGT48;// = 29.0/72.0;// (Phase III)
	double	pr_ETR_TxExp_RGT48;// = 40.0/72.0;// (Phase III)

	double	pr_TxFailTW12_TxExp_SOC; // probability of detectable HCV-RNA at TW 12 (treatment failure) in SOC
	double	pr_TxFailTW12_TxExp_PRB48;// probability of detectable HCV-RNA at TW 12 (treatment failure) in BOC arm 3
	double	pr_ResponseTW8_TxExp;// = 74.0/162.0; // probability of detectable HCV-RNA at TW 12 (treatment failure) in RGT 36 week therapy
	double	pr_ResponseNoneTW8_TxExp;// = 72.0/162.0; // probability of detectable HCV-RNA at TW 12 (treatment failure) in RGT 36 week therapy
	double	pr_TxFailTW12_TxExp_RGT36;// = 0.0/74.0; // probability of detectable HCV-RNA at TW 12 (treatment failure) in RGT 36 week therapy
	double	pr_TxFailTW12_TxExp_RGT48;// = 32.0/72.0; // probability of detectable HCV-RNA at TW 12 (treatment failure) in RGT 48 week therapy

	double	pr_AEAnemia_SOC; // probability of getting Anemia on PEG+RIB = 0.34 (P03523:Table 27)
	double	pr_AEAnemia_BOC; // probability of getting Anemia on Boceprevir = 0.56 (P03523:Table 27)
	double	pr_AEAnemia_RGT; //
	double	pr_AERash_SOC; // probability of getting rash on PEG+RIB = 0.37 (P03523:Table 30)
	double	pr_AERash_BOC; // probability of getting rash on Bocepravir = 0.40 (P03523:Table 30)

	double	pr_AEAnemia_TxExp_SOC;// = 16.0/60.0; //
	double	pr_AEAnemia_TxExp_RGT;// = 70.0/162.; //
	double	pr_AEAnemia_TxExp_PRB48;// = 74.0/161.; //
	double	pr_ContinueTx_TxExp_BOC;
	double	pr_ContinueTx_TxExp_SOC;
	double	pr_ContinueTx_TxExp_RGT;

	double	pr_EpoUse_2W;
	double	pr_EpoUse_8W;
	double	pr_EpoUse_18W;
	double	pr_EpoUse_30W;
	double	pr_EpoUse_42W;


	vector<double> _fibDistr,_genoDistr, _mfDistr;
	double _p8Week_LDV, _pMale, _pTE, _pIFN_intol, _pBOC;
	



	// added Qiushi @ 11/28/2016 for Acute HCV treatment
	int _wks_acute_svr12;
	int _wks_acute_wait_and_see; 
	double _pr_acute_spont_clearance;
	double _pr_svr_for_SA, _pr_svr_for_SA_acute, _pr_svr_for_SA_F0;
	double _pr_acute_prob_complete_tx;
	double _pr_chronic_prob_follow_up;

	// added Qiushi @ 03/14/2017 for R1 India CEA manuscript
	bool _flag_override_pefct_tx;

	 // added Qiushi @ 03/31/2017, update the meta-regression progression risks
	double _ageHCVacquisition, _iduProp;
	bool _flagUseUpdatedMetaRegressionForFibProgr;

	// ------- added Qiushi @ 6/19/2017 for HCC screening model ---------------
	map<type_scr_fib_state, double> _map_inci_HCC;
	double _pr_monthly_hcc_small_med;
	double _pr_monthly_hcc_med_large;
	double _mort_hcc_large;
	double _mort_F4;
	double _mort_DC;

	double _mort_tx_transpl_1y, _mort_tx_transpl_1yPlus;
	double _mort_tx_res_small, _mort_tx_res_med, _mort_tx_res_large;
	double _mort_tx_abl_cc_small, _mort_tx_abl_cc_med, _mort_tx_abl_cc_large;
	double _mort_tx_abl_dc_small, _mort_tx_abl_dc_med, _mort_tx_abl_dc_large;
	double _mort_tx_pal;

	double _sens_surveillance;
	double _spec_surveillance;
	double _sens_diagnostic;

	double _pr_transpl_hcv;
	double _pr_transpl_hcc;
	
};


class costType{
public:
	costType();
	~costType(){};
	const double GetDrugCost(drugType theDrug);

	double c_acute;
	double c_F0;	//	costs assoc w\ F0 stage fibrosis"
	double c_F1;	//	cumulative costs assoc w\ subj in F1 stage"
	double c_F2;	//	cum costs assoc w\ F2 stage"
	double c_F3;	//	cum costs assoc w/ F3 stage"
	double c_CoCirr;	//	cost associated with compensated cirrhosis"
	double c_DeCirr;	//	first year cost associated with decompensated cirrhosis"
	double c_DeCirr1yrPlus; // subsequent year cost associated with decompensated cirrhosis
	double c_HCC;	//	cost associated with hepatocellular carcinoma"
	double c_LivTr;	//	cost associated with liver tranplantation"
	double c_PostLivTr;	//	cost associated with post-liver transplant"
	double c_ETR;	//	cost associated with ETR"
	double c_SE_Boc;	//	cost of treating side effects assoc with Bocepravir"
	double c_SVR;	//	cost associated with SVR state"
	double c_PEG; 	//  weekly cost of peginterferon;"
	double c_RBV; 	//  weekly cost of ribavirin;"
	double c_BOC;		//	weekly cost of Bocepravir"
	double c_SOF;		//	weekly cost of Sofosbuvir"
	double c_TEL;     //	weekly cost of Telaprevir"
	double c_SMV;     //  weekly cost of Simeprevir
	double c_Epo;		// weekly cost of treating Anemia with ESA"

	double c_LDV;		// weekly cost of "ladipasvir"

	double c_DCV;
	double c_PrOD;

	double c_ASV;


	// testing cost
	// added 11/23/2016
	double c_testing_RNA;
	double c_testing_genotype;
	double c_testing_antiHCV;
	double c_testing_preTx;
	double c_testing_postTx;

	// annual background cost
	// added 03/08/2017 (for R1 India CEA manuscript)
	double c_background;

	// for disease burden model
	map<int, double> map_drug_cost_discount;
	map<int, double> _table_tx_cost;

	double override_DAA_cost_value;
	int override_DAA_cost_since_when;
	double scalingCoef;
	double c_screening;



	// ------ for HCC screening analysis, added 6/20/2017 ---------
	double hccscr_c_tx_transpl_1y, hccscr_c_tx_transpl_1yPlus;
	double hccscr_c_tx_res, hccscr_c_tx_abl, hccscr_c_tx_pal;
	double hccscr_c_test_surveillance, hccscr_c_test_diagnostic;


};

class qolType
{
public:
	qolType();
	~qolType(){};

	double q_acute;

	double q_F0;
	double q_F1;
	double q_F2;
	double q_F3;
	double q_CoCirr; //	health related quality of life associated with compensated cirrhosis  Siebert, Gut 2003, 52:425-432
	double q_DeCirr; //	health related quality of life associated with decompensated cirrhosis  Siebert, Gut 2003, 52:425-432
	double q_SVR; //health related quality of life associated with post-SVR  Siebert, Gut 2003, 52:425-432
	double q_ETR; //health related quality of life associated with end of treatment response
	double q_HCC; //health related quality of life associated with hepatocellular carcinoma  Siebert, Gut 2003, 52:425-432
	double q_LivTr; //	health related quality of life associated with liver transplantation  Siebert, Gut 2003, 52:425-432
	double q_PostLivTr; //	health related quality of life associated with liver transplantation after 1 year Siebert, Gut 2003, 52:425-432
	double q_TX_oSOC; // treatment-related qol multiplier
	double q_TX_DAA; // treatment-related qol multiplier
	double q_Dec_Anemia; // QOL decrement factor associated with Anemia

	//		static double q_PEG_Rib_Boc; //	health related quality of life associated with triple drug therapy
	//	    static double q_PEG_Rib; //	health related quality of life associated with PEG + Rib drug therapy
	//	    static double q_Boc_noSE; //	QoL while on Bocepravir and not experiencing side effects
	//	    static double q_Boc_SE; //	QoL while on Bocepravir and experiencing side effects
	//  qolType(); //constructor

	// ------ for HCC screening analysis, added 6/20/2017 ---------
	double q_tx_transpl_1y, q_tx_transpl_1yPlus ;
	double qChange_tx_res_1y, qChange_tx_res_1yPlus;
	double qChange_tx_abl_1y, qChange_tx_abl_1yPlus;
	double qChange_tx_pal_1y, qChange_tx_pal_1yPlus;

};


class txProfileType{
public:
	txProfileType();
	~txProfileType(){};
	double GetTxCostWeek(costType & argCost, int txWeek);

	bool _oSOC; // if this treatment is old standard-of-care (or DAA)
	double _pr_TxComplete;
	double _pr_ContinueTx;
	double _pr_AEAnemia ;
	double _pr_ETRgivenEOT;
	double _pr_SVRgivenETR;
	int _txDuration;

	int _assignedEpoUse; // total assigned time for Epo use for anemia
	int _beginEpoUse; //week to begin treatment for anemia with Epo
	double _adjustEpoUse; //adjustment for Epo use to incorporate AE discontinuations

	double _pr_ETRgivenTxDiscontinue;
	double _pr_TxFailTW12;

	map<drugType,pair<int,int> > _tx_timeline; //[drugType --> (start week, end week)]

};
class counterType
{
public:
	counterType(); // constructor
	~counterType(){};

	vector<int>	aliveCount		;
	vector<int>	incidentDeathCount		;
	vector<int>	incidentDeathLiv		;
	vector<int>	txCountDouble		;
	vector<int>	txCountTriple		;
	vector<int>	txCountNone		;



	vector<int>	incidentDisTxDouble		;
	vector<int>	incidentDisTxTriple		;
	vector<int>	incidentStateCountF0		;
	vector<int>	incidentStateCountF1		;
	vector<int>	incidentStateCountF2		;
	vector<int>	incidentStateCountF3		;
	vector<int>	incidentStateCountCoCirr		;
	vector<int>	incidentStateCountDeCirr		;
	vector<int>	incidentStateCountDeCirr1yrPlus		;
	vector<int>	incidentStateCountHCC		;
	vector<int>	incidentStateCountLivTr		;
	vector<int>	incidentStateCountLivTr1yrPlus		;
	vector<int>	incidentStateCountCured		;
	vector<int>	incidentStateCountETR		;

	vector<int>	prevStateCountF0		;
	vector<int>	prevStateCountF1		;
	vector<int>	prevStateCountF2		;
	vector<int>	prevStateCountF3		;
	vector<int>	prevStateCountCoCirr		;
	vector<int>	prevStateCountDeCirr		;
	vector<int>	prevStateCountDeCirr1yrPlus		;
	vector<int>	prevStateCountHCC		;
	vector<int>	prevStateCountLivTr		;
	vector<int>	prevStateCountLivTr1yrPlus		;
	vector<int>	prevStateCountCured		;
	vector<int>	prevStateCountETR		;

	vector<int>	prevAEAnemia		;
	vector<int>	incidentAEAnemia		;


	//int aliveCount[MAXCYCLE]; // count alive cases in each cycle for WebModel ouput (max array size = 52*95, i.e. max years lived after entering model is 95 years. age at death can be max 110.)
	//int incidentDeathCount[MAXCYCLE]; // count death incidents in each cycle for WebModel ouput (max array size = 52*95)
	//int incidentDeathLiv[MAXCYCLE];
	//int txCountDouble[MAXCYCLE];
	//int txCountTriple[MAXCYCLE];
	//int txCountNone[MAXCYCLE];

	////   double q_AgeNormalAvg[MAXCYCLE];

	//int incidentDisTxDouble[MAXCYCLE];
	//int incidentDisTxTriple[MAXCYCLE];
	//int incidentStateCountF0[MAXCYCLE];
	//int incidentStateCountF1[MAXCYCLE];
	//int incidentStateCountF2[MAXCYCLE];
	//int incidentStateCountF3[MAXCYCLE];
	//int incidentStateCountCoCirr[MAXCYCLE];
	//int incidentStateCountDeCirr[MAXCYCLE];
	//int incidentStateCountDeCirr1yrPlus[MAXCYCLE];
	//int incidentStateCountHCC[MAXCYCLE];
	//int incidentStateCountLivTr[MAXCYCLE];
	//int incidentStateCountLivTr1yrPlus[MAXCYCLE];
	//int incidentStateCountCured[MAXCYCLE];
	//int incidentStateCountETR[MAXCYCLE];

	//int prevStateCountF0[MAXCYCLE];
	//int prevStateCountF1[MAXCYCLE];
	//int prevStateCountF2[MAXCYCLE];
	//int prevStateCountF3[MAXCYCLE];
	//int prevStateCountCoCirr[MAXCYCLE];
	//int prevStateCountDeCirr[MAXCYCLE];
	//int prevStateCountDeCirr1yrPlus[MAXCYCLE];
	//int prevStateCountHCC[MAXCYCLE];
	//int prevStateCountLivTr[MAXCYCLE];
	//int prevStateCountLivTr1yrPlus[MAXCYCLE];
	//int prevStateCountCured[MAXCYCLE]; 
	//int prevStateCountETR[MAXCYCLE]; 

	//int prevAEAnemia[MAXCYCLE];
	//int incidentAEAnemia[MAXCYCLE];

	int countDeCirr; //counter for the person with a history of decompensated cirrhosis
	int countHCC; // counter for history of hepatocellular carcinoma
	int countLivTr; // counter for number of liver transplants
	int countTxEx; // counter for previous treatment experienced patients
	int countNR; // counter for non-responders to drug (Bocepravir or MK-7009)
	int countReL; // counter for  person who attains ETR and then relapses
	int countDeathLiv; // counter for  deaths due to liver-related complications

};


class pplProfileType
{
public:
	pplProfileType() {};
	~pplProfileType() {};

	map<stateType, double> _distr_fib;
	map<string, double> _distr_genotype;
	map<char, double> _distr_gender;


};

class DALYType
{
public:
	DALYType();
	~DALYType() {};

	double _discountRate;	// discount rate usually 3%/year
	double _beta;			// age weighting constant, =0.04
	double _K;				// age-weighting modulation constant, K=0 for no age-weighting, 1 for age weighting
	double _C;				// adjustment constant for age-weights (C = 0.1658)

	map<int, double> _lifeExpct_male;
	map<int, double> _lifeExpct_female;

	//enum stateType {s_F0, s_F1, s_F2, s_F3, s_CoCirr, s_DeCirr, s_DeCirr1yrPlus, s_HCC, s_LivTr, s_LivTr1yrPlus, s_SVR, s_Death}; 
	double _dw_f0;
	double _dw_f1;
	double _dw_f2;
	double _dw_f3;
	double _dw_CoCirr;
	double _dw_DeCirr;
	double _dw_DeCirr1yrPlus;
	double _dw_HCC;
	double _dw_LivTr;
	double _dw_LivTr1yrPlus;
	double _dw_SVR;

	double GetYLD(double argAge, int argState, double argL);
	double GetYLL(double argAge, char argGender);
	double GetDisutilityWeight(int argState);

	int ReadLifeExpectancyData(string argFile);

};


class NonNHANESType {
public: 
	NonNHANESType(){
		_mapNonNHANESCohortName.clear();
		_mapNonNHANESCohortName.insert(pair<string,typeCohortDisBurdnModel>("Incarcerated",cohort_incarcerated));
		_mapNonNHANESCohortName.insert(pair<string,typeCohortDisBurdnModel>("Homeless",cohort_homeless));
		_mapNonNHANESCohortName.insert(pair<string,typeCohortDisBurdnModel>("Hospitalized",cohort_hospitalized));
		_mapNonNHANESCohortName.insert(pair<string,typeCohortDisBurdnModel>("NursingHome",cohort_nursingHome));
		_mapNonNHANESCohortName.insert(pair<string,typeCohortDisBurdnModel>("IndianRsv",cohort_indianReservation));

		_incarcerated_release_prob = 0.45;
		_incacerated_jail_pct = 0.33;
		_bcscr_coverage = 0.75; // assumption for now
	};
	~NonNHANESType(){};
	typeCohortDisBurdnModel GetCohortType(string argStr) const {
		map<string, typeCohortDisBurdnModel>::const_iterator it = _mapNonNHANESCohortName.find(argStr);
		if(it != _mapNonNHANESCohortName.end()){
			return it->second;
		}else{
			ExitWithMsg("Error @ nonNHANESType: unknown cohort type = " + argStr);
			
		}
	}
	double GetNonNHANESPplData(string argStr, typeCohortDisBurdnModel argCohort) const {
		map<string, map<typeCohortDisBurdnModel, double>>::const_iterator it = _mapNonNHANESPplData.find(argStr);
		if (it == _mapNonNHANESPplData.end()) {
			ExitWithMsg("ERROR @ GetNonNHANESPplData: unknown parameter name: " + argStr);
		}
		map<typeCohortDisBurdnModel, double>::const_iterator it2 = it->second.find(argCohort);
		if (it2 == it->second.end()) {
			ExitWithMsg("ERROR @ GetNonNHANESPplData: unknown cohort index: " + basicToStr((int)argCohort));
		}
		// if the return value = -1, then show the error message
		if ( abs(it2->second+1) < EPSILON) {
			ExitWithMsg("ERROR @ GetNonNHANESPplData: undefined value for parameter "+ argStr +" and population " + typeCohort_label[(int)argCohort]);
		}
		return it2->second;
	}


	map<string, typeCohortDisBurdnModel> _mapNonNHANESCohortName;
	map<string, map<typeCohortDisBurdnModel, double>> _mapNonNHANESPplData;
	map<typeCohortDisBurdnModel, vector<int>> _distr_age_initial_category_by_cohort;
	map<typeCohortDisBurdnModel, vector<int>> _distr_age_new_inci_category_by_cohort;
	map<typeCohortDisBurdnModel, vector<double>> _distr_age_initial_by_cohort;	
	map<typeCohortDisBurdnModel, vector<double>> _distr_age_new_inci_by_cohort;
	map<int, map<typeCohortDisBurdnModel,int>> _txCap_by_year_cohort;
	map<int, map<typeCohortDisBurdnModel,double>> _scrRate_by_year_cohort; 

	double _bcscr_coverage;

	map<typeCohortDisBurdnModel, map<int, double>> _distr_genotype_by_cohort;
	map<typeCohortDisBurdnModel, map<typeStateDisBurdnModel, double>> _distr_fib_by_cohort;
	map<typeCohortDisBurdnModel, map<typeInsurance, double>> _distr_insur_by_cohort;
	
	double _incarcerated_release_prob;
	double _incacerated_jail_pct;


	// ***** immigrant parameters **********
	map<int, map<typeStateDisBurdnModel, map<int, int>>> _immigrants_hcv_num_by_year_fib_genotype; // [year][fib][genotype]
	double _immigrants_male_pct;
	vector<int> _immigrants_age_category;
	vector<double> _immigrant_distr_age;
	int _immigrants_year_start;
};

class psaDistrType
{
public:
	psaDistrType(){};
	~psaDistrType(){};
	int ReadDistrForPSA(string argFileStr);

	vector<string> _distrPSA_varName;
	vector<string> _distrPSA_distrName;
	vector<double> _distrPSA_param1, _distrPSA_param2;

};


//==== patientProfilet: may not be used ===
struct patientProfile{
	float bAge;
	int bState;
	char bGender;
};
//struct patientProfilet{
//	float bAget;
//	int bStatet;
//	char bGendert;
//};
//extern patientProfilet patProfilet[400];


struct SParamTriplet{
public:
	double _base, _lb,_ub;
};


map<string,SParamTriplet> ReadOneWaySAParamValueRange(string argFile);
int ReadOneWaySAParamValueRange(string argFile, vector<string> & listVarName, vector<SParamTriplet> & listVal);


//////////////////////////////////////////////////////////////////////////
/////////   FOLLOWING IS DEFINED FOR THE HCC SCREENING MODEL ////////////
//////////////////////////////////////////////////////////////////////////
class HCCScreeningDataType
{
public:
	HCCScreeningDataType() {};
	~HCCScreeningDataType() {};

	double _initial_age;
	map<int, double> _distr_genotype;
	double _distr_male;
	map<type_scr_fib_state, double> _distr_initial_fib;

	map<int, double> _bgMort_male, _bgMort_female;
	map<int, double> _qol_male, _qol_female;

	vector<type_scr_absb_state> _list_hcc_tx;
	map<string, vector<double>> _distr_hcc_tx;
	
};


class HCCScrPatientType	// [QC] used in simulations (simulated patients)
{
public:
	HCCScrPatientType() {
		_currentAge = 40;
		_initialAge = 40;
		_gender = 'M';
		_genotype = 1;
		_stateFib = scr_fib_F3;
		_stateHCC = scr_hcc_none;
		_stateAbs = scr_abs_none;

		_stateFibWhenTreated = scr_fib_F3;
		_stateHCCWhenTreated = scr_hcc_small;
		_cycles_on_HCC_treatment = 0;
		_flag_started_transpl = false;



	}	// assign default value in the constructor
	~HCCScrPatientType() {};

	double _currentAge, _initialAge;
	char _gender;
	int _genotype; //values = {1,2,3}
	type_scr_fib_state _stateFib, _stateFibWhenTreated;
	type_scr_hcc_state _stateHCC, _stateHCCWhenTreated;
	type_scr_absb_state _stateAbs;

	int _cycles_on_HCC_treatment;
	bool _flag_started_transpl;
};



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

class modelParamType
{
public:
	modelParamType(string argName) { _armName = argName; };
	modelParamType(string argName, baseCohortType argBaseCohort) { _armName = argName; _cohortData = argBaseCohort; };
	modelParamType(baseCohortType argBaseCohort) { _cohortData = argBaseCohort; };
	modelParamType() {};
	~modelParamType() {};

	int ReduceDrugCost(double reductPerct);


	baseCohortType _cohortData;
	transitionType _transData;
	costType _costData;
	qolType _qolData;
	//txProfileType _txData;
	pplProfileType _pplData;
	DALYType _dalyData;


	string _armName;	// [QC] unique identify of the arm (treatment)



						/////// For global analysis //////////////////
	int ApplyIndiaParameters();


	/////// For acute HCV treatment //////////////
	string _armName_acute; // 2016/11/28 added for acute HCV phase 
	inline void SetAcuteTx(string argName) { _armName_acute = argName; }
	int _acute_tx_duration, _f0_tx_duration;
	inline void SetAcuteTxDuration(int argV) { _acute_tx_duration = argV; }
	inline void SetF0TxDuration(int argV) { _f0_tx_duration = argV; }

	// =======================================
	// definitions for disease burden model
	// ---------------------------------------
	disBurdnDataType _disBurdnData;
	string _disBurdnOutFile;

	NonNHANESType _nonNHANESPplData;

	HCCScreeningDataType _hccScreeningData;

};



#endif

