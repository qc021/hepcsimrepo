//#ifdef _DEBUG
//#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#define new DEBUG_NEW
//#endif

//#include<crtdbg.h>
//#define _CRTDBG_MAP_ALLOC

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>
#include<cassert> //header file for function "assert()"

//#include"project_doublecheck.h"
//#include"project_ageNoTx.h"
#include"project_psa_metamodeling.h"
//#include"project_priortzTx.h"
//#include"project_cea_prison.h"
//#include"project_hcvcalculator.h"
//#include"project_3wk.h" // added: 12/12/2015
//#include"project_tfs.h" // added: 4/4/2016
#include"project_india.h" // added: 7/9/2016
#include"project_acute.h" // added 11/22/2016
#include"project_disburdn_txfailure.h" // added 12/22/2016
#include"project_disburdn_eu5salvage.h" // added 03/07/2017
#include"project_disburdn_pakistan.h"  // added 04/23/2017
#include"project_disburdn_us_nonnhanes.h" // added 05/15/2017
#include"project_disburdn_statex.h"// added 06/7/2017
#include"project_cea_hcc_scr.h"	// added 6/16/2017
#include"project_global.h" // added 7/18/2017

using namespace std;
///*************** DON'T delete the following declaration ***********************/
///*************** MUST keep the declaration of following extern variables ******/
bool TREATMENT; // to swicth-off treatment, change value to 0
long SIM_SEED;
long NUM_PATIENTS;
long NUM_RUNS_PSA;
unsigned long PSA_SEED;
long NUM_RUNS_PSA_ONE_BATCH;
int PRINT_SCREEN;
int PROGRESSION_AGE_GENDER; //use 1 for fibrosis progression rates to vary with age and gender
int PSA_OPTION;
double TIME_HORIZON;
bool EXPECTED_LIFE_ONLY;
string FILE_background_mortality;
double DISCOUNT_FACT_C;
double DISCOUNT_FACT_Q;
double DISCOUNT_FACT_DALY;
double DISCOUNT_FACT_L;

bool ALLOW_RETREATMENT, OVERRIDE_NOTREATMENT, OVERRIDE_SYNC_RND_STRM;
int ONEWAYSA_OPTION_PROJ_ACUTE;

// global variables declared for disease burden variables
bool FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION;
bool FLAG_SMOOTH_RESULTS;

void Analysis_CEA_DoubleCheck() {
	//// =======================================================================================
	//// project 1
	//// =======================================================================================
	// double check the new code

	//Project_Doublecheck myProject;
	//myProject.SetDrugCostReduction(-0.11); //-0.25, -0.44
	//PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	//myProject.CEA_BaseResults();


	//myProject.SetDrugCostReduction(-0.23); //-0.25, -0.44
	//PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	//int idx = atoi(argv[1]);
	//myProject.SetOutFile("output_project_doublecheck/output_change_horizon" + basicToStr(idx) + ".txt");
	//myProject.CEA_Horizon(idx);

	//PROGRESSION_AGE_GENDER = 0; 
	//int idx=atoi(argv[1]);	// comparator index
	//myProject.CEA_OneWay(idx);


	//PROGRESSION_AGE_GENDER = 0; 
	//int idx=atoi(argv[1]);	// comparator index
	//int batch=atoi(argv[2]);
	//int bsize=atoi(argv[3]);
	//NUM_PATIENTS = 10000; //number of patients generated in simulation
	//myProject.CEA_PSA(idx,batch,bsize);
}

void Analysis_CEA_Prison() {
	// =======================================================================================
	// project 1.5 cost-effectiveness in prison
	// =======================================================================================

	//Project_CEA_Prison myProject;
	//myProject.SetDrugCostReduction(-0.23); //

	//PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	//myProject.CEA_BaseResults();


	//PROGRESSION_AGE_GENDER = 0; 
	//int idx=atoi(argv[1]);	// comparator index
	//myProject.CEA_OneWay(idx);


	//PROGRESSION_AGE_GENDER = 0; 
	//int idx=atoi(argv[1]);	// comparator index
	//int batch=atoi(argv[2]);
	//int bsize=atoi(argv[3]);
	//NUM_PATIENTS = 10000; //number of patients generated in simulation
	//myProject.CEA_PSA(idx,batch,bsize);


}

void Analysis_Prioritize_Age() {
	// =======================================================================================
	// project 2
	// =======================================================================================
	// Find the age threshold
	//NUM_PATIENTS = 10000; //number of patients generated in simulation
	////NUM_RUNS_PSA = 1000;
	//Project_AgeNoTx myProject_ageForTx;
	//myProject_ageForTx.SetInputScenarioFile("project_ageNoTx_input_scenarios.txt");

	// NOTE: SEQUENTIAL APPROACH IS TOO SLOW!!!!
	//double listDisc[]={};
	//double listDisc[]={0,-0.15,-0.25,-0.48,-0.25,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6};
	//for(int d=0; d<sizeof(listDisc)/sizeof(double); d++){
	//	myProject_ageForTx.ReduceDrugCost(listDisc[d]);
	//	myProject_ageForTx.SetOutputDetail(true);
	//	myProject_ageForTx.BaseCase_aggregated(); // default parameters: double argAgeStart=30.0, double argAgeEnd=86.0, double argAgeIncr=2.0
	//	myProject_ageForTx.SetOutputDetail(false);
	//}

	// ****************************
	// Base case results 
	// ****************************

	//// NOTE: USE THE PARALLEL APPROACH
	//double a=30+2*atof(argv[1]);	// start 0,1,...,28	
	//double listDisc[]={-0.48,-0.25,-0.1,-0.15,-0.2,-0.3,-0.4,-0.5,-0.6,0};
	//for(int d=0; d<sizeof(listDisc)/sizeof(double); d++){
	//	myProject_ageForTx.ReduceDrugCost(listDisc[d]);
	//	myProject_ageForTx.BaseCase_distributed(a,a+2,2);
	//	
	//}
	// ****************************
	// One-way sensitivity analysis
	// ****************************
	//myProject_ageForTx.ReduceDrugCost(-0.25);
	//PROGRESSION_AGE_GENDER = 0; 
	//myProject_ageForTx.DSA_OneWay_ICER_for_age(74);
	//myProject_ageForTx.DSA_OneWay_ICER_for_age(76);
	//myProject_ageForTx.DSA_OneWay_ICER_for_age(80);

	//myProject_ageForTx.ReduceDrugCost(-0.48);
	//PROGRESSION_AGE_GENDER = 0; 
	//myProject_ageForTx.DSA_OneWay_ICER_for_age(80);
	//myProject_ageForTx.DSA_OneWay_ICER_for_age(82);
	//myProject_ageForTx.DSA_OneWay_ICER_for_age(86);

	// ****************************
	// PSA: ICER vs. age
	// ****************************
	////-------------------
	//// 10,000MC 5,000PSA

	//NUM_PATIENTS = 10000; //number of patients generated in simulation
	//NUM_RUNS_PSA = 5000;
	//PROGRESSION_AGE_GENDER = 0; 

	//myProject_ageForTx.ReduceDrugCost(-0.25);
	//double a1=atof(argv[1]);				// 30, 35, 40, ..., 85
	//NUM_RUNS_PSA_ONE_BATCH=atol(argv[2]);	// 250 for each batch, 20 batches in total
	//int a2=atoi(argv[3]);		//	$(Process): 0-19
	//int a3=3*(2+atoi(argv[3]));
	//myProject_ageForTx.PSA_ICER_for_each_age(a1,a2,a3); //age, batch_id, seed


	// ****************************
	// PSA: cut-off age
	// ****************************
	////-------------------
	//NUM_PATIENTS = 10000; //number of patients generated in simulation
	//NUM_RUNS_PSA = 5000;
	//int i=atoi(argv[1]); // 3*(2+atoi(argv[1]));	// 0 - 99
	//NUM_RUNS_PSA_ONE_BATCH=50;
	//myProject_ageForTx.ReduceDrugCost(-0.25);
	//myProject_ageForTx.PSA_age_threshold(100,i);

}

void Analysis_Metamodeling(int argc, const char * argv[]) {
	//// =======================================================================================
	//// project 3
	//// =======================================================================================
	Project_PSA_MetaModeling myProject;
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1000;
	PSA_SEED = 1;
	myProject.SetDrugCostReduction(-0.11);

	//// ------------ base case run -------------------	
	//myProject.RunBaseCase();

	// ------------ PSA run -------------------
	PSA_OPTION = 1;	// turn on the PSA mode

	// [Note] Output files : "./output_project_metamodel/project_psa_output_[replcPSA]PSA_[replcMC]MC_SEED[PSA_SEED].txt"
	PSA_SEED = 1;

	int replcMC = 1;	// initialization	
	int replcPSA = 1;	// initialization

	//int replcPSA = 100;		NUM_PATIENTS = 10000;
	//int replcPSA = 1000;	NUM_PATIENTS = 1000;
	//int replcPSA = 10000;	NUM_PATIENTS = 100;
	//int replcPSA = 100000;	NUM_PATIENTS = 10;


	//myProject.SetPrintParamValues(false); // default value = false;

	//replcPSA = 100000;	NUM_PATIENTS = 1;		
	//myProject.RunPSA_Multivariate_OutputSingleMCRun(replcPSA, replcMC);

	//replcPSA = 10000;	NUM_PATIENTS = 10;
	//myProject.RunPSA_Multivariate_OutputSingleMCRun(replcPSA, replcMC);

	//replcPSA = 1000;	NUM_PATIENTS = 100;
	//myProject.RunPSA_Multivariate_OutputSingleMCRun(replcPSA, replcMC);

	//replcPSA = 100;	NUM_PATIENTS = 1000;
	//myProject.RunPSA_Multivariate_OutputSingleMCRun(replcPSA, replcMC);


	// ------------ Change a particular parameter -------------------
	//myProject.SetPrintParamValues(false);	// default value = false;
	//myProject.RunModelWithDiffParamValue();





	//PSA_OPTION=1;	// turn on the PSA mode
	//myProject.RunPSA_Univariate_OutputSingleMCRun(1000,1);
	//myProject.RunPSA_Univariate_OutputSingleMCRun(100,10);

	// ------------ Running gold standard -------------------
	// running in batch 1-1000
	myProject.SetPrintParamValues(true); // default value = false;
	PSA_SEED = 999 + atoi(argv[1]);
	replcPSA = 1;	NUM_PATIENTS = 100000;
	myProject.RunPSA_Multivariate_OutputSingleMCRun(replcPSA, replcMC);






	// ***************************** NOTE ************************************************
	// PSA_SEED: random seed of PSA samplers (independent of seed for microsimulation)
	// PSA_MAX_SAMPLE_IN_SINGLE_OUTPUTFILE: in "sim_config.h", default value is 100,000: start a new PSA output file for each 100,000 lines of results;
	// 
	// arguments of main function:
	// a1: PSA samples
	// a2: MC samples: total number of MC runs = a1 * a2
	// a3: can be used as index and random seed for PSA at the same time.
	// For example: 100,000PSA 10MC, we can use following settings:
	//   a1=2,000; 
	//	 a2=10;
	//	 a3=0,1,2,3,...,49 (50 different sessions are running simultaneously [parallel computing])
	//	 Output files: (each has a1 * a2 lines, = 20,000 lines in this example)
	//		"./output_PSA_metamodel/project_psa_output_2000PSA_10MC_SEED0.txt"
	//		"./output_PSA_metamodel/project_psa_output_2000PSA_10MC_SEED1.txt"
	//										...
	//		"./output_PSA_metamodel/project_psa_output_2000PSA_10MC_SEED49.txt"
	// ***********************************************************************************

	//a1=100;
	//a2=100;
	//a3=0;

	// [additional remark: to avoid (undesirable) correlation between MC runs (same patient ID, different drug), change:
	// USE_COMMON_RANDOM_NUMBER=false; 
	// and use all random number from stream 1


}

void Analysis_Treatment_Prioritization() {
	//// =======================================================================================
	//// project 4 Tx prioritization
	//// =======================================================================================

	//NUM_PATIENTS = 10000; //number of patients generated in simulation
	//PROGRESSION_AGE_GENDER = 0; //use 1 for fibrosis progression rates to vary with age and gender
	//Project_priortzTx myProjectPriotz;
	//myProjectPriotz.BaseCase();


	////myProjectPriotz.OneWaySA(50,52,TN_TOL,s_F3);
	////myProjectPriotz.OneWaySA(50,52,TN_TOL,s_CoCirr);
	////myProjectPriotz.OneWaySA(50,52,TE_NA,s_F3);
	////myProjectPriotz.OneWaySA(50,52,TE_NA,s_CoCirr);
}

void Analysis_HCV_calculator() {
	////// =======================================================================================
	////// project 5 HCV calculator
	////// =======================================================================================
	//int arg1=atoi(argv[1])*5+30;	// argv[1]=0-9 ==> arg1=30-75
	//int arg2=atoi(argv[2]);
	//Project_HCVCalculator myCalculator;
	//PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	//myCalculator.SetInitialAge(arg1);
	//
	//if(arg2==0) myCalculator.RunBaseCase();

	//myCalculator.RunPerturbations("project_hcvcalculator_input_pert_values.txt", arg2);


	////Project_HCVCalculator myCalculator;
	////PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	////myCalculator.RunBaseCase();
	////myCalculator.RunPerturbations("project_hcvcalculator_input_pert_values.txt");


}

void Analysis_HCV_3wk_Regimen() {
	// =======================================================================================
	// project 6: Cost-effectiveness: HCV 3wk treatment
	// =======================================================================================

	//PRINT_SCREEN=1;
	//Project_3Wk myProject;
	//myProject.SetDrugCostReduction(-0.23); //-0.25, -0.44

	//PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	//myProject.CEA_BaseResults();


	////PROGRESSION_AGE_GENDER = 0; 
	////int idx=atoi(argv[1]);	// comparator index
	////myProject.CEA_OneWay(idx);


	////PROGRESSION_AGE_GENDER = 0; 
	////int idx=atoi(argv[1]);	// comparator index
	////int batch=atoi(argv[2]);
	////int bsize=atoi(argv[3]);
	////NUM_PATIENTS = 10000; //number of patients generated in simulation
	////myProject.CEA_PSA(idx,batch,bsize);
}

void Analysis_Get_Transpl_free_survival() {
	// =======================================================================================
	// project 7: Get transplant free survival
	//	Female/Male: 90/10%
	//	Age 54.5
	//	Total time: 15 years
	//	Run for F1, F2, F3, F4.
	// =======================================================================================
	//// Remember to change:
	//// const double TIME_HORIZON = 15; // horizon = 15 years

	//Project_TFS myProject;

	//PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	//myProject.BaseCase();
}

void Analysis_CEA_India() {
	//=======================================================================================
	//project 8: Cost-effectiveness analysis of hepatitis C treatment in India
	//[1] use fixed fibrosis progression probability (not meta regression)

	//=======================================================================================
	PROGRESSION_AGE_GENDER = 0; //use 1 for fibrosis progression rates to vary with age and gender
	//FILE_background_mortality = "Input_mortality_male_female_India.in"; //WRONG, but used in the first submission
	FILE_background_mortality = "Input_mortality_male_female_India_correct.in";
	FILE_background_mortality = "Input_mortality_male_female_WHO.in"; // FOR DALY RESULTS ONLY
	Project_India myProject;
	myProject.SetDrugCostReduction(0.0); //-0.25, -0.44
	 //-------------------
	 //QALY results
	 //-------------------

	 //base case result [can run on desktop]	


	EXPECTED_LIFE_ONLY = 0;
	//myProject.SetOutFile("project_india_output_compare_QALY.txt");
	myProject.SetOutFile("project_india_output_compare_QALY_0602.txt");
	myProject.CEA_BaseResults();

	//EXPECTED_LIFE_ONLY = 1;
	//myProject.SetOutFile("project_india_output_compare_LY_0316.txt");
	//myProject.CEA_BaseResults();

	//// also evaluate the following settings:
	//// discount: 0%, 5%

	//DISCOUNT_FACT_C = 0; //annual discount factor for costs
	//DISCOUNT_FACT_Q = 0; //anual discount factor for QALYs
	//EXPECTED_LIFE_ONLY = 0;
	//myProject.SetOutFile("project_india_output_compare_QALY_0316_disc0.txt");
	//myProject.CEA_BaseResults();

	//EXPECTED_LIFE_ONLY = 1;
	//myProject.SetOutFile("project_india_output_compare_LY_0316_disc0.txt");
	//myProject.CEA_BaseResults();

	//DISCOUNT_FACT_C = 0.05; //annual discount factor for costs
	//DISCOUNT_FACT_Q = 0.05; //anual discount factor for QALYs
	//EXPECTED_LIFE_ONLY = 0;
	//myProject.SetOutFile("project_india_output_compare_QALY_0316_disc5.txt");
	//myProject.CEA_BaseResults();

	//EXPECTED_LIFE_ONLY = 1;
	//myProject.SetOutFile("project_india_output_compare_LY_0316_disc5.txt");
	//myProject.CEA_BaseResults();

	//DISCOUNT_FACT_C = 0.03; //annual discount factor for costs
	//DISCOUNT_FACT_Q = 0.03; //anual discount factor for QALYs
	//EXPECTED_LIFE_ONLY = 0;
	//myProject.SetOutFile("project_india_output_compare_QALY_0316_disc3.txt");
	//myProject.CEA_BaseResults();

	//EXPECTED_LIFE_ONLY = 1;
	//myProject.SetOutFile("project_india_output_compare_LY_0316_disc3.txt");
	//myProject.CEA_BaseResults();



	////
	//// analysis on horizon/age
	////
	//int idx = atoi(argv[1]);	// comparator index	
	//myProject.SetOutFile("output_project_india/project_india_output_change_horizon" + basicToStr(idx) + ".txt");
	//myProject.CEA_VaryTimeHorizon(idx);

	//int idx = atoi(argv[1]);	// comparator index	
	//myProject.SetOutFile("output_project_india/project_india_output_change_age_and_horizon" + basicToStr(idx) + ".txt");
	//myProject.CEA_VaryAgeAndTimeHorizon(idx);

	//EXPECTED_LIFE_ONLY = 1;
	//myProject.SetOutFile("output_project_india/project_india_output_change_age.txt");
	//myProject.CEA_VaryAge(0);



	//myProject.SetOutFile("output_project_india/project_india_output_change_drugcost.txt");
	//myProject.CEA_VaryDrugCost(0);


	////
	//// sensitivity analysis
	////

	//PROGRESSION_AGE_GENDER = 0; 
	//int idx=atoi(argv[1]);	// comparator index
	//myProject.CEA_OneWay(idx, "project_india_input_OneWay.txt");
	//myProject.CEA_OneWay(idx, "project_india_input_OneWay_DALYParam.txt");

	//// idx = 0,1,...,99; bsize=100
	//PROGRESSION_AGE_GENDER = 0; 
	//
	//PSA_OPTION=1;	// turn on the PSA mode
	//NUM_RUNS_PSA = 5000;
	//int idx=atoi(argv[1]);	// comparator index
	//int nBatches=atoi(argv[2]);
	//PSA_SEED = idx;
	//myProject.CEA_PSA(idx, nBatches);



	// ------------------ additional SA requested for R1 -------------
	//myProject.CEA_SA_BackgroundCost();

	//myProject.CEA_R1_LiverRelatedMortality();
	//myProject.CEA_R1_discount();
	//myProject.CEA_R1_LifeExpectancy_w_wo_HCV();
	//myProject.CEA_R1_progRisk();
	//myProject.CEA_R1_progRisk_and_Horizon();

	// -------------------
	// DALY results
	//FILE_background_mortality = "Input_mortality_male_female_WHO.in"; // FOR DALY RESULTS ONLY
	// -------------------

	//EXPECTED_LIFE_ONLY = 0;
	//myProject.SetOutFile("project_india_output_compare_DALY.txt");
	//myProject.CEA_BaseResults();




	//PROGRESSION_AGE_GENDER = 0; 
	//int idx=atoi(argv[1]);	// comparator index
	//myProject.CEA_OneWay(idx, "project_india_input_OneWay_DALYParam.txt");



	//// idx = 0,1,...,99; bsize=100
	//PROGRESSION_AGE_GENDER = 0; 
	//
	//PSA_OPTION=1;	// turn on the PSA mode
	//NUM_RUNS_PSA = 5000;
	//int idx=atoi(argv[1]);	// comparator index
	//int nBatches=atoi(argv[2]);
	//PSA_SEED = idx;
	//myProject.CEA_PSA(idx, nBatches);

}


void Analysis_CEA_treating_acute_HCV(int argc, const char * argv[]) {
	// =======================================================================================
	// project 9: Acute or delay? 
	OVERRIDE_SYNC_RND_STRM = true;
	DISCOUNT_FACT_C = 0.03; //annual discount factor for costs
	DISCOUNT_FACT_Q = 0.03; //anual discount factor for QALYs
	EXPECTED_LIFE_ONLY = 0;
	NUM_PATIENTS = 100000;
	ONEWAYSA_OPTION_PROJ_ACUTE = 0; // default val = 0, no oneway SA

	FILE_background_mortality = "Input_mortality_male_female.in";
	//FILE_background_mortality = "Input_mortality_male_female_2011.in";
	//FILE_background_mortality = "Input_mortality_male_female_noBgMort.in";
	//ALLOW_RETREATMENT = false; [default value = false, in the beginning of main file]
	SIM_SEED = 1;
	// =======================================================================================


	Project_Acute myProject;
	myProject.SetInitialAge(26);
	myProject.SetMaleRatio(0.64);
	myProject.SetDrugCostReduction(0); // no discount considered in base case
	PROGRESSION_AGE_GENDER = 0; //use 1 for fibrosis progression rates to vary with age and gender

	ALLOW_RETREATMENT = true;



	//// ------------------- experiment on MC runs --------------------------------------
	//// 50,000 * {1,2,3,..., 40}
	//int idx = atoi(argv[1]); //0, ... , 79
	////int idx = 80+atoi(argv[1]); //0, ... , 79

	//SIM_SEED = idx / 40 + 1;
	//NUM_PATIENTS = (idx % 40 + 1) * 5000;
	////NUM_PATIENTS =  (atoi(argv[1]) + 1) * 5000;
	//
	//myProject.SetOutFile("output_project_acute/project_acute_output_varyMC"+basicToStr(idx)+".txt");	
	//myProject.CEA_BaseResults();



	//------------------- base case --------------------------------------
	DISCOUNT_FACT_C = 0; //annual discount factor for costs
	DISCOUNT_FACT_Q = 0; //anual discount factor for QALYs
	//EXPECTED_LIFE_ONLY = 1;

	//DISCOUNT_FACT_C = 0.03; //annual discount factor for costs
	//DISCOUNT_FACT_Q = 0.03; //anual discount factor for QALYs
	EXPECTED_LIFE_ONLY = 0;

	NUM_PATIENTS = 100000;
	//SIM_SEED = 1;// atol(argv[1]) + 1;	
	myProject.SetOutFile("output_project_acute/project_acute_output_QALY_disc0_0504.txt");
	//myProject.SetOutFile("output_project_acute/project_acute_output_QALY_disc3_"+basicToStr(SIM_SEED)+".txt");	
	myProject.CEA_BaseResults();



	//cout<<endl<<endl<<"----- Change treatment durations ----------"<<endl;
	////DISCOUNT_FACT_C = 0; //annual discount factor for costs
	////DISCOUNT_FACT_Q = 0; //anual discount factor for QALYs
	////EXPECTED_LIFE_ONLY = 1;
	//myProject.SetOutFile("output_project_acute/project_acute_output_change_bothTxDuration.txt");
	//myProject.CEA_VaryTxDuration_Acute_ChronicF0(0);


	////
	//// analysis on horizon/age
	////
	//cout << endl << endl << "----- Change initial age ----------" << endl;
	//myProject.SetOutFile("output_project_acute/project_india_output_change_age.txt");
	//myProject.CEA_VaryAge(0);

	//cout<<endl<<endl<<"----- Change time horizons ----------"<<endl;
	//myProject.SetOutFile("output_project_acute/project_india_output_change_horizon.txt");
	//myProject.CEA_VaryTimeHorizon(0);


	cout << endl << endl << "----- Change access probability ----------" << endl;
	ALLOW_RETREATMENT = true;
	myProject.SetOutFile("output_project_acute/project_acute_output_change_accessProb_chronicArm.txt");
	myProject.CEA_Vary_AccessProbability(0);



	//cout << endl << endl << "----- Change SVR vs. Clearance rate; by cost discount----------" << endl;
	//myProject.SetOutFile("output_project_acute/project_acute_output_change_cost_clr_svr.txt");
	//myProject.CEA_Vary_TxCost_SVRacute_ClrRate(atoi(argv[1])-1); // with index: 1-11


	//cout << endl << endl << "----- One way sensitivity ----------" << endl;	
	//PROGRESSION_AGE_GENDER = 0;
	//myProject.CEA_OneWay(0, "project_acute_input_OneWay.txt", atoi(argv[1])-1); // with index 1-20


	//// ------------------- PSA --------------------------------------

	//cout<<endl<<endl<<"----- PSA ----------"<<endl;
	//PROGRESSION_AGE_GENDER = 0; 	
	//PSA_OPTION=1;	// turn on the PSA mode
	//NUM_PATIENTS = 10000;
	//NUM_RUNS_PSA = 10000;
	//int idx_psa=atoi(argv[1])-1;	// comparator index = 1-200 
	//int nBatches=atoi(argv[2]);		// = 200
	//PSA_SEED = idx_psa;
	//myProject.CEA_PSA(idx_psa, nBatches);


	// ************* unused analysis *******************************************************************************************
	//cout<<endl<<endl<<"----- Change Cost discount vs. SVR ----------"<<endl;
	//myProject.SetOutFile("output_project_acute/project_acute_output_change_cost_svr.txt");
	//myProject.CEA_Vary_TxCost_SVRacute(0);

	//cout << endl << endl << "----- Change Cost discount vs. Clearance rate ----------" << endl;
	//myProject.SetOutFile("output_project_acute/project_acute_output_change_cost_clr.txt");
	//myProject.CEA_Vary_TxCost_ClearanceRate(0);



	//cout<<endl<<endl<<"----- Change Duration of Acute Treatment & Prob of Self-clearance ----------"<<endl;
	//int k = 1;
	//myProject.SetOutFile("output_project_acute/project_acute_output_change_acuteTxDuration_prClr"+basicToStr(k)+".txt");
	//myProject.CEA_VaryAcuteTxDuration_probSelfClearance(k);


	//myProject.SetOutFile("output_project_acute/project_acute_output_change_acuteTxDuration_SVR" + basicToStr(k) + ".txt");
	//myProject.CEA_VaryAcuteTxDuration_SVR(k);


}

void Analysis_US_Salvage() {
	// ************ IMPORTANT NOTES ************************************************************
	// [1] to replicate the previous results of [output_project_burden/out_basecase_us_costburden_0605.txt]
	//     might want to comment out the following line in the sim_disburdn.cpp:
	//     newPt->_flagInterferonTol = _rnd_pt.GetU01() < argModelParam._disBurdnData._prob_IFNintol ? false : true; // added for cost burden model @ 6/5/2017
	// [2] to replicate the previous results of [out_basecase_us_costburden_0716.txt]
	//     just run the code as is



	//// =======================================================================================
	//// project 10: Treatment failure - US
	FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION = false;
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1;
	// double-check version: const bool SMOOTH_RESULTS = true;
	//// =======================================================================================

	Project_TxFailure prjTxFailure;
	prjTxFailure.SetDrugCostReduction(-0.11); //-0.25, -0.44
	prjTxFailure.SetScreeningOption(typeScreeningScenario::screen_birthCohort);
	//prjTxFailure.SetScreeningOption(typeScreeningScenario::screen_universal_rate);
	prjTxFailure.SetOutFile("output_project_burden/out_basecase_us_costburden_0720.txt");
	prjTxFailure.Initialize_US_Base();
	prjTxFailure.BaseCaseRun();
}

void Analysis_EU5_Salvage() {

	//// =======================================================================================
	//// project 11: Treatment Salvage - EU5
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1;
	//// =======================================================================================

	//string strCountry = "FR";
	//string strCountry = "DE";
	//string strCountry = "IT";
	//string strCountry = "ES";
	//string strCountry = "UK";

	string listCountry[] = { "FR","DE","IT","ES","UK" };
	int idx;
	cout << "select country: 0=France, 1=Germany, 2=Italy, 3=Spain, 4=UK: ";
	cin >> idx;
	//for (int idx = 0; idx < sizeof(listCountry) / sizeof(string); idx++) {

	string strCountry = listCountry[idx];

	Project_SalvageEU5 prjSalvageEU5;
	prjSalvageEU5.SetDrugCostReduction(0); //-0.25, -0.44
	//prjSalvageEU5.SetScreeningOption(typeScreeningScenario::screen_birthCohort);
	//prjSalvageEU5.SetScreeningOption(typeScreeningScenario::screen_universal_rate);
	prjSalvageEU5.SetScreeningOption(typeScreeningScenario::screen_no_extra);
	prjSalvageEU5.SetOutFile("output_project_burden/out_basecase_EU5_" + strCountry + ".txt");
	prjSalvageEU5.Initialize_Model_Parameter(strCountry);
	prjSalvageEU5.BaseCaseRun();

	//}

}


void Analysis_Disease_Burden_Pakistan(int argc, const char * argv[]) {
	//// =======================================================================================
	//// project 12: Disease burden in Pakistan

	// !!!!!!!!!!!!!!! REMINDER: MAKE SURE TO HAVE !!!!!!!!!!!!!!!!!!!!!!!!!!!
	// #define __PAKISTAN_SETTING__  in disburdn_config.h
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION = true;
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1;
	//// =======================================================================================

	Project_DisBurd_Pakistan prjPak;
	prjPak.SetDrugCostReduction(0); //-0.25, -0.44
	//prjSalvageEU5.SetScreeningOption(typeScreeningScenario::screen_birthCohort);
	//prjPak.SetScreeningOption(typeScreeningScenario::screen_universal_rate);
	prjPak.SetScreeningOption(typeScreeningScenario::screen_universal_capacity);
	//prjPak.SetScreeningOption(typeScreeningScenario::screen_no_extra);

	//// ---------- calibrate the incidence proportion and intiial popualtion size --------------
	//int idx = atoi(argv[1]);	// comparator index; run on ERIS ONE server
	//prjPak.CalibInitialSizeAndInciRate(idx);


	//// ---------- base case run --------------
	//prjPak.SetOutFile("output_project_pakistan/out_basecase.txt");	
	//prjPak.Initialize_Pakistan_Parameter();
	//prjPak.BaseCaseRun();

	//// ---------- identify sufficient intervention to achieve WHO target --------------
	//int idx = atoi(argv[1]);	// comparator index
	//prjPak.TestWHOTarget(idx);

	//// ---------- re-evaluate with WHO-target intervention, change incidence scenarios ----------------------------	 
	prjPak.SetOutFile("output_project_pakistan/out_basecase.txt");
	prjPak.Initialize_Pakistan_Parameter(type_pakistan_incidence::INCI_CONSTANT);
	//prjPak.Initialize_Pakistan_Parameter(type_pakistan_incidence::INCI_INCREASING);
	prjPak.BaseCaseRun();



}


void Analysis_Disease_Burden_US_nonNHANES(int argc, const char * argv[]) {
	//// =======================================================================================
	//// project 12: Disease burden in US, including NHANES and nonNHANES population

	// !!!!!!!!!!!!!!! REMINDER: MAKE SURE TO HAVE !!!!!!!!!!!!!!!!!!!!!!!!!!!
	// #define __DEFAULT_SETTING__  in disburdn_config.h
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION = true;
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1;
	//// =======================================================================================

	Project_DisBurd_UScomprh prjNonNHANES;
	prjNonNHANES.SetDrugCostReduction(0); //-0.25, -0.44

	prjNonNHANES.SetScreeningOption(typeScreeningScenario::screen_universal_rate);
	//prjNonNHANES.SetScreeningOption(typeScreeningScenario::screen_birthCohort);
	//prjNonNHANES.SetScreeningOption(typeScreeningScenario::screen_universal_rate);
	//prjNonNHANES.SetScreeningOption(typeScreeningScenario::screen_no_extra);


	//int idx = atoi(argv[1]);	// comparator index
	//prjPak.CalibInitialSizeAndInciRate(idx);


	//// ---------- base case run --------------
	prjNonNHANES.SetOutFile("output_project_nonnhanes/out_basecase.txt");
	prjNonNHANES.Initialize_UScomprh_Parameter();
	prjNonNHANES.BaseCaseRun();



}




void Analysis_Elimin_HCV_State_X(int argc, const char * argv[]) {
	//// =======================================================================================
	//// project 12: Disease burden in Pakistan

	// !!!!!!!!!!!!!!! REMINDER: MAKE SURE TO HAVE !!!!!!!!!!!!!!!!!!!!!!!!!!!
	// const bool SMOOTH_RESULTS = false;
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FLAG_SMOOTH_RESULTS = false;
	FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION = true;
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1;
	//// =======================================================================================

	Project_DisBurd_StateX prjStateX;
	prjStateX.SetDrugCostReduction(0); //-0.25, -0.44

	prjStateX.SetScreeningOption(typeScreeningScenario::screen_universal_capacity);
	//prjStateX.SetScreeningOption(typeScreeningScenario::screen_birthCohort);
	//prjPak.SetScreeningOption(typeScreeningScenario::screen_universal_rate);
	//prjPak.SetScreeningOption(typeScreeningScenario::screen_no_extra);




	//// ---------- base case run --------------
	//prjStateX.SetOutFile("output_project_statex/out_basecase.txt");
	////prjPak.Initialize_StateX_Parameter(type_pakistan_incidence::INCI_PROPORTIONAL);
	////prjPak.Initialize_StateX_Parameter(type_pakistan_incidence::INCI_CONSTANT);
	//prjStateX.Initialize_StateX_Parameter(type_statex_incidence::STATEX_INCI_CONSTANT);
	//
	//prjStateX.BaseCaseRun();




	//// ---------- identify sufficient intervention to achieve WHO target --------------
	prjStateX.ROI_Analysis(0);

	//// ---------- Compare costs -------------------------------------------------------
	//prjStateX.CompareCost();


}

void Analysis_CEA_HCC_screening(int argc, const char * argv[]) {
	//// =======================================================================================
	//// project 14: HCC screening: cost-effectiveness analysis

	FLAG_SMOOTH_RESULTS = false;
	FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION = true;
	PROGRESSION_AGE_GENDER = 0;
	SIM_SEED = 1;
	//// =======================================================================================

	NUM_PATIENTS = 100000;
	DISCOUNT_FACT_L = 0.03;
	Project_CEA_HCC_Screening prjHCC;
	prjHCC.SetDrugCostReduction(0); //-0.25, -0.44

	//// ---------- base case run --------------
	prjHCC.Initialize();
	prjHCC.SetOutFile("output_project_hcc_scr/out_basecase.txt");
	prjHCC.CEA_BaseResults();
}



void Analysis_CEA_Global(int argc, const char * argv[]) {
	////=======================================================================================
	////project 15: Cost-effectiveness analysis for global cost-effectiveness
	////=======================================================================================
	PROGRESSION_AGE_GENDER = 0; //use 1 for fibrosis progression rates to vary with age and gender
								//FILE_background_mortality = "Input_mortality_male_female_India.in"; //WRONG, but used in the first submission
	FILE_background_mortality = "Input_mortality_male_female_India_correct.in";
	//FILE_background_mortality = "Input_mortality_male_female_WHO.in"; // FOR DALY RESULTS ONLY
	Project_CEA_Global myProject;
	myProject.SetDrugCostReduction(0.0); //-0.25, -0.44
	
	myProject.SetVersionNote("v0");

	//int idx = 0;
	int idx = atoi(argv[1])-1;
	myProject.RunCEAForListedCountries(9);
	





}
///******************************************************************************/

int main(int argc, const char * argv[]) {


	////cout<<endl<<endl;
	//CRndVarGen rvGen1;
	//rvGen1.Initialize(0,UNIFREAL,0,1);
	//map<string,CRndVarGen> sampler;
	//sampler["SP1"]=rvGen1;
	//for(int k=0; k<20; k++){
	//	cout << sampler["SP1"].GetRV();
	//}


	// *************** MUST assign value to following extern variables ******
	TREATMENT = 1; // to swicth-off treatment, change value to 0
	SIM_SEED = 0;
	NUM_PATIENTS = 10000; //number of patients generated in simulation
	NUM_RUNS_PSA = 10000;
	PSA_SEED = 1;
	NUM_RUNS_PSA_ONE_BATCH = 100;
	PRINT_SCREEN = 0;
	PROGRESSION_AGE_GENDER = 1; //use 1 for fibrosis progression rates to vary with age and gender
	PSA_OPTION = 0;
	TIME_HORIZON = 150.0;
	EXPECTED_LIFE_ONLY = 0; // making it equal to 1 will calculate expected life years instead of QALYs
	FILE_background_mortality = "Input_mortality_male_female.in";

	DISCOUNT_FACT_C = 0.03; //annual discount factor for costs
	DISCOUNT_FACT_Q = 0.03; //anual discount factor for QALYs
	DISCOUNT_FACT_DALY = 0.00; //anual discount factor for QALYs
	DISCOUNT_FACT_L = 0.00;

	ALLOW_RETREATMENT = false;
	OVERRIDE_NOTREATMENT = false;
	OVERRIDE_SYNC_RND_STRM = false;

	// global variables for disease burden model
	FLAG_UNIF_INT_RND_GENERATOR_ROUND_VERSION = true;
	FLAG_SMOOTH_RESULTS = true;
	// **********************************************************************


	clock_t start = clock();
	cout << "==== Hello Hep C =====" << endl << endl;

	//// Project 1
	//Analysis_CEA_DoubleCheck();

	//// Project 1.5 cost-effectiveness in prison
	//Analysis_CEA_Prison();

	//// Project 2 age threshold
	//Analysis_Prioritize_Age();

	//// Project 3
	//Analysis_Metamodeling(argc, argv);


	//// Project 4: treatment prioritization
	//Analysis_Treatment_Prioritization();

	//// Project 5: HCV calculator
	//Analysis_HCV_calculator();

	//// Project 6: 3wek treatment
	//Analysis_HCV_3wk_Regimen();

	//// Project 7: Get Transplant free survival
	//Analysis_Get_Transpl_free_survival();


	//// Project 8: CEA India
	//Analysis_CEA_India();


	//// Project 9: Treatment failures in US
	//Analysis_US_Salvage();


	//// Project 10: acute HCV treatment
	//Analysis_CEA_treating_acute_HCV(argc, argv);


	//// Project 11: Treatment failures in EU5
	Analysis_EU5_Salvage();




	// Project 12: Disease burden in Pakistan
	// MAKE SURE TO CHANGE THIS: #define __PAKISTAN_SETTING__  @ disburdn_config.h
	//Analysis_Disease_Burden_Pakistan(argc, argv);


	// Project 12: Disease burden in Pakistan
	// MAKE SURE TO CHANGE THIS: #define __DEFAULT_SETTING__  @ disburdn_config.h
	//Analysis_Disease_Burden_US_nonNHANES(argc, argv);


	// Project 13: HCV elimination in state X
	//Analysis_Elimin_HCV_State_X(argc, argv);


	// Project 14: HCC screening cost-effectiveness analysis
	//Analysis_CEA_HCC_screening(argc, argv);

	// Project 15: HCV global CEA
	//Analysis_CEA_Global(argc, argv);


	cout << endl << endl << "Program running for " << Time(start) << " sec...";
	//getchar();

//	_CrtDumpMemoryLeaks();

	return 0;
}