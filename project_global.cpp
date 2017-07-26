#include "project_global.h"
using namespace std;



Project_CEA_Global::Project_CEA_Global()
{

	_age_initial = 35.0; // base case average age.
	_drugCostReduction = 0;
	//// read all arms for comparisons
	//ifstream inf;
	//inf.open("project_india_input_comparators.txt");
	//string word;
	//while (inf >> word) {
	//	_listCmp.push_back(word);
	//	inf >> word; _listArm1.push_back(word);
	//	inf >> word; _listArm2.push_back(word);
	//	int g;
	//	inf >> g;	_listGenotype.push_back(g);
	//}
	//inf.close();

	// read the parameter files for each country
	// read all arms for comparisons
	ifstream inf;
	inf.open("project_global_input_param_by_country.txt");
	if (inf.fail()) {
		ExitWithMsg("ERROR @ Project_CEA_Global(): can't open file project_global_input_param_by_country.txt");
	}
	string word;
	string line;
	while (std::getline(inf, line))
	{

		std::stringstream   linestream(line);
		getline(linestream, word, '\t');

		if ("//" == word) {		
			continue;
		}
		map<string, double> m;
		linestream >> m["distrF0"]
			>> m["distrF1"]
			>> m["distrF2"]
			>> m["distrF3"]
			>> m["distrF4"]
			>> m["distrG1"]
			>> m["distrG2"]
			>> m["distrG3"]
			>> m["distrG456"]
			>> m["distrMale"]
			>> m["c_DAA_wk"]
			>> m["c_F0"]
			>> m["c_F1"]
			>> m["c_F2"]
			>> m["c_F3"]
			>> m["c_F4"]
			>> m["c_DC"]
			>> m["c_HCC"]
			>> m["c_preTxTesting"]
			>> m["c_postTxTesting"];



		_map_countryParam[word] = m;
	}
	inf.close();


	// read background mortality	
	inf.open("project_global_input_life_table_by_country.txt");
	if (inf.fail()) {
		ExitWithMsg("ERROR @ Project_CEA_Global(): can't open file project_global_input_life_table_by_country.txt");
	}	
	// ======================================
	// NOTE: life table is for 5-year age groups
	// need to use .lowerbound() function, 
	// instead of .find() function to retrive the background mortality
	// ======================================
	_map_country_life_table_female.clear();
	_map_country_life_table_male.clear();
	_map_qol_female.clear();
	_map_qol_male.clear();

	map<int, double> table_mFemale, table_mMale, table_qFemale, table_qMale;

	
	// read table header, register countries
	vector<string> listCountry, listLabel;
	string t;
	getline(inf, t);
	istringstream iss(t);
	
	map<int, double> emptymap;
	while (iss >> word) {
		listCountry.push_back(word);
		_map_country_life_table_male[word] = emptymap;
		_map_country_life_table_female[word] = emptymap;
	}
	getline(inf, t);
	

	double val;
	while (inf >> val) {
		int age = (int)val;

		inf >> _map_qol_male[age];
		inf >> _map_qol_female[age];

		for (int k = 0; k < listCountry.size(); k=k+2) {
			inf >> val;
			_map_country_life_table_male[listCountry[k]][age] = val;
			inf >> val;
			_map_country_life_table_female[listCountry[k]][age] = val;
		}
	}

	inf.close();


	inf.open("project_global_input_country_list.txt");
	if (inf.fail()) {ExitWithMsg("ERROR @ Project_CEA_Global(): can't open file project_global_input_country_list.txt");	}
	_listCountries.clear();
	ReadList(_listCountries, inf);
	inf.close();


}

void Project_CEA_Global::RunCEAForOneCountry(string argCountry, int argIdx)
{
	_id_batch = argIdx;
	EXPECTED_LIFE_ONLY = 0;
	LoadParamByCountry(argCountry);

	// backup the model parameters for the country
	_baselineModelParam_country = _modelParam;



	ofstream outf;
	string strOutFile = "output_project_global/out_by_country/" + _strVerMark + "_" + argCountry + "_" + basicToStr(argIdx)+".txt";
	outf.open(strOutFile.c_str());
	if (outf.fail()) { ExitWithMsg("ERROR @ RunCEAForONeCOuntry: Can't open file: " + strOutFile); }
	

	string listCostParam[] = { "c_DAA_wk","c_F0","c_F1","c_F2","c_F3","c_F4","c_DC","c_HCC","c_preTxTesting","c_postTxTesting" };
	double listPertPct[] = { 0, -0.5, 0.5 };

	int counter = 0;

	
		for (int idxParam = 0; idxParam < sizeof(listCostParam) / sizeof(string); idxParam++) {
			for (int idxPert = 0; idxPert < sizeof(listPertPct) / sizeof(double); idxPert++) {
				// run base case results only at counter == 0
				if (idxPert == 0) {
					if (counter != 0) {
						continue;
					}
				}

				// reset the model parameters
				_modelParam = _baselineModelParam_country;
				// change one cost parameter
				ChangeValue(listCostParam[idxParam], _map_countryParam[argCountry][listCostParam[idxParam]] * (1.0 + listPertPct[idxPert]), _modelParam);


				FILE_background_mortality = _str_file_temp_life_table;
				SetOutFile("output_project_global/temp_out/" + _strVerMark + "_" + argCountry + "_QALY_details.txt");
				vector<double> val_QALY = GetAggregatedResults();

				// ********** need re-run the model to get DALY results ***********************************************
				// Skip this step for now
				// -----------------------
				//FILE_background_mortality = "Input_mortality_male_female_WHO.in"; // FOR DALY RESULTS ONLY	
				//SetOutFile("output_project_global/temp_out/" + _strVerMark + "_" + argCountry + "_DALY_details.txt");
				//vector<double> val_DALY = GetAggregatedResults();
				vector<double> val_DALY = val_QALY;


				//if (argIdx == 0 && counter == 0) {
				//	// print header of the table
				//	outf << "country\tidx\tage\thorizon\t";
				//	outf << "c_DAA_wk\tc_F0\tc_F1\tc_F2\tc_F3\tc_F4\tc_DC\tc_HCC\tc_preTxTesting\tc_postTxTesting\t";
				//	outf << "QALY1\tCost1\tDC1\tHCC1\tLivTr1\tDeathLiv1\tTxCost1\tYLL1\tYLD1\tDALY1\tTestCost1\tLY1\t";
				//	outf << "QALY2\tCost2\tDC2\tHCC2\tLivTr2\tDeathLiv2\tTxCost2\tYLL2\tYLD2\tDALY2\tTestCost2\tLY2\t";
				//	outf << endl;s
				//}
				outf << fixed << setprecision(3)
					<< argCountry <<"\t" << counter << "\t" << _age_initial << "\t" <<TIME_HORIZON<<"\t"
					// -------- print out the cost parameters --------------
					<< _modelParam._costData.c_SOF << "\t"
					<< _modelParam._costData.c_F0 << "\t"
					<< _modelParam._costData.c_F1 << "\t"
					<< _modelParam._costData.c_F2 << "\t"
					<< _modelParam._costData.c_F3 << "\t"
					<< _modelParam._costData.c_CoCirr << "\t"
					<< _modelParam._costData.c_DeCirr << "\t"
					<< _modelParam._costData.c_HCC << "\t"
					<< _modelParam._costData.c_testing_preTx << "\t"
					<< _modelParam._costData.c_testing_postTx << "\t";

				// -------- print out the model outputs ---------------
				////[ 0] r.push_back(mySim->GetAvgQALY());
				////[ 1] r.push_back(mySim->GetAvgCost());
				////[ 2] r.push_back(mySim->GetCounter().countDeCirr);
				////[ 3] r.push_back(mySim->GetCounter().countHCC);
				////[ 4] r.push_back(mySim->GetCounter().countLivTr);
				////[ 5] r.push_back(mySim->GetCounter().countDeathLiv);
				////[ 6] r.push_back(mySim->GetTxCost());
				////[ 7] r.push_back(mySim->GetDALY_YLL());
				////[ 8] r.push_back(mySim->GetDALY_YLD());
				////[ 9] r.push_back(mySim->GetDALY_Total());
				////[10] r.push_back(mySim->GetTestCost());
				////[11] r.push_back(mySim->GetAvgLY());

				for (int k = 0; k < val_QALY.size(); k++) {
					if ((k >= 7 && k <= 9)
						|| (k >= 7 + 11 && k <= 9 + 11)) {
						outf << val_DALY[k] << "\t";
					}
					else {
						outf << val_QALY[k] << "\t";
					}
				}
				outf << endl;

				counter++;
			}// end for - perturbation 
		}// end for cost parameters



	outf.close();

}

void Project_CEA_Global::RunCEAForListedCountries(int argIdx)
{

	double listAge[] = { 25,30,35,40,45,50,55,60,65,70 };	 // # = 10
	double listHorz[] = { 1,2,3,5,7,10,20,30,40,150 }; // # = 10
	int nAge = sizeof(listAge) / sizeof(double);
	int nHorz = sizeof(listHorz) / sizeof(double);

	int idxCountry = argIdx / (nAge * nHorz);
	int idxAge = (argIdx % (nAge * nHorz)) / nHorz;
	int idxHorz = (argIdx % (nAge * nHorz)) % nHorz;

	cout <<endl<< "============ Calculating CEA of HCV treatment for countr = " << _listCountries[idxCountry] << ", initial age = " << listAge[idxAge] << ", horizon = " << listHorz[idxHorz] << " ========== "<<endl<<endl;


	SetInitialAge(listAge[idxAge]);
	TIME_HORIZON = listHorz[idxHorz];
	RunCEAForOneCountry(_listCountries[idxCountry], argIdx % (nAge * nHorz));

}


int Project_CEA_Global::LoadParamByCountry(string argCountry)
{

	if (_map_countryParam.find(argCountry) == _map_countryParam.end()
		|| _map_country_life_table_female.find(argCountry) == _map_country_life_table_female.end()
		|| _map_country_life_table_male.find(argCountry) == _map_country_life_table_male.end()) {
		ExitWithMsg("ERROR @ LoadParamCountry(): Can't find data for country = " + argCountry);
	}
	


	// Generate the correct by-each-age life table.
	_str_file_temp_life_table = "output_project_global/temp_out/temp_" + argCountry + "_input_life_table"+basicToStr(_id_batch)+".txt";
	ofstream outf(_str_file_temp_life_table.c_str());
	if (outf.fail()) {
		ExitWithMsg("ERROR @ fail to write to file: "+_str_file_temp_life_table);
	}
	for (int age = 0; age <= 120; age++) {
		outf << age << "\t";
		double p_male = _map_country_life_table_male[argCountry].lower_bound(age)->second;
		double p_female = _map_country_life_table_female[argCountry].lower_bound(age)->second;

		if (age == 0) {
			// do  nothing
		}
		else if (age >= 1 && age <= 4) {
			p_male = 1.0 - pow(1.0 - p_male, 1.0 / 4.0);
			p_female = 1.0 - pow(1.0 - p_female, 1.0 / 4.0);
		}
		else {
			p_male = 1.0 - pow(1.0 - p_male, 1.0 / 5.0);
			p_female = 1.0 - pow(1.0 - p_female, 1.0 / 5.0);
		}
		outf<< p_male<<"\t"
			<< p_female<<"\t"
			<< _map_qol_male.lower_bound(age)->second << "\t"
			<< _map_qol_female.lower_bound(age)->second << endl;
	}
	outf.close();


	// reset the parameter values
	_modelParam = _baselineModelParam;



	// [standard across countries] 
	// used for calculating the fibrosis progression risks
	_modelParam._transData._flagUseUpdatedMetaRegressionForFibProgr = true;
	_modelParam._transData._ageHCVacquisition = 23.5; // Ahmad study, lowest age observed in F0 category
	_modelParam._transData._iduProp = 0.41; // Ahmad study, lowest age observed in F0 category


	// modifying cost parameters by country
	
	_modelParam._costData.c_testing_antiHCV = 0;
	_modelParam._costData.c_testing_RNA = 0;
	_modelParam._costData.c_testing_genotype = 0;
	_modelParam._costData.c_testing_preTx = _map_countryParam[argCountry]["c_preTxTesting"]; // = 98 + 435 + 35; // RNA, gentoype, antiHCV testing, respectively
	_modelParam._costData.c_testing_postTx = _map_countryParam[argCountry]["c_postTxTesting"]; // =98; // RNA testing

	// annual cost	
	_modelParam._costData.c_acute = 0;
	_modelParam._costData.c_F0 = _map_countryParam[argCountry]["c_F0"];
	_modelParam._costData.c_F1 = _map_countryParam[argCountry]["c_F1"];
	_modelParam._costData.c_F2 = _map_countryParam[argCountry]["c_F2"];
	_modelParam._costData.c_F3 = _map_countryParam[argCountry]["c_F3"];
	_modelParam._costData.c_CoCirr = _map_countryParam[argCountry]["c_F4"];
	_modelParam._costData.c_DeCirr = _map_countryParam[argCountry]["c_DC"];
	_modelParam._costData.c_DeCirr1yrPlus = _map_countryParam[argCountry]["c_DC"];
	_modelParam._costData.c_HCC = _map_countryParam[argCountry]["c_HCC"];
	_modelParam._costData.c_LivTr = 0; //								
	_modelParam._costData.c_PostLivTr = 0; //								
	_modelParam._costData.c_ETR = 0;
	_modelParam._costData.c_SE_Boc = 0; 
	_modelParam._costData.c_SVR = 0;

	// weekly cost
	// we only use SOF for simplifity
	// keep other drug cost = 0
	_modelParam._costData.c_SOF = _map_countryParam[argCountry]["c_DAA_wk"];  ////weekly cost
	_modelParam._costData.c_PEG = 0;//	588.00; //weekly cost								
	_modelParam._costData.c_RBV = 0;//	309.00; //weekly cost							
	_modelParam._costData.c_BOC = 0;// 1100.00; //weekly cost								
	_modelParam._costData.c_TEL = 0;// 4100.00;	// weekly cost									
	_modelParam._costData.c_SMV = 0;
	_modelParam._costData.c_LDV = 0;
	_modelParam._costData.c_DCV = 0;
	_modelParam._costData.c_PrOD = 0;
	_modelParam._costData.c_ASV = 0;
	


	// turn off transplant for all countries
	_modelParam._transData.pr_DeCirr_LivTr = 0;						
	_modelParam._transData.pr_HCC_LivTr = 0;
	

	// distribution in the population
	// fibrosis, genotype, and sex
	_modelParam._pplData._distr_fib[s_F0] = _map_countryParam[argCountry]["distrF0"];
	_modelParam._pplData._distr_fib[s_F1] = _map_countryParam[argCountry]["distrF1"];
	_modelParam._pplData._distr_fib[s_F2] = _map_countryParam[argCountry]["distrF2"];
	_modelParam._pplData._distr_fib[s_F3] = _map_countryParam[argCountry]["distrF3"];
	_modelParam._pplData._distr_fib[s_CoCirr] = _map_countryParam[argCountry]["distrF4"];
	_modelParam._pplData._distr_genotype["G1"] = _map_countryParam[argCountry]["distrG1"];
	_modelParam._pplData._distr_genotype["G2"] = _map_countryParam[argCountry]["distrG2"];
	_modelParam._pplData._distr_genotype["G3"] = _map_countryParam[argCountry]["distrG3"];
	_modelParam._pplData._distr_genotype["G4"] = _map_countryParam[argCountry]["distrG456"];
	_modelParam._pplData._distr_gender['M'] = _map_countryParam[argCountry]["distrMale"];
	_modelParam._pplData._distr_gender['F'] = 1 - _map_countryParam[argCountry]["distrMale"];

	// ************************************************
	// ------- do NOT modify the following ------------
	// parameters for calculating DALYS
	// standard across all countries
	_modelParam._dalyData._discountRate = 0.0;
	_modelParam._dalyData._beta = 0.04;
	_modelParam._dalyData._C = 0.1658;
	_modelParam._dalyData._K = 0;
		
	_modelParam._dalyData.ReadLifeExpectancyData("./Input_life_expectancy_WHO.txt");
	
	_modelParam._dalyData._dw_f0 = 0;
	_modelParam._dalyData._dw_f1 = 0;
	_modelParam._dalyData._dw_f2 = 0;
	_modelParam._dalyData._dw_f3 = 0;
	_modelParam._dalyData._dw_CoCirr = 0;
	_modelParam._dalyData._dw_DeCirr = 0.194;
	_modelParam._dalyData._dw_DeCirr1yrPlus = 0.194;
	_modelParam._dalyData._dw_HCC = 0.508;
	_modelParam._dalyData._dw_LivTr = 0;
	_modelParam._dalyData._dw_LivTr1yrPlus = 0;
	_modelParam._dalyData._dw_SVR = 0;
	// ------- end of DALY parameter settings  -------
	// ************************************************

	return 0;
}

int Project_CEA_Global::CEA_BaseResults()
{


	Compare("NoTx", "G1_SOF_LDV12", 1);
	//Compare("G1_SOF_LDV12", "NoTx", 1);
	//Compare("NoTx","G1_SOF_DCV12_24",1);

	//Compare("NoTx","TN_G1_LDP_SOF12",11);	
	//Compare("NoTx","TN_TOL_G1a_DCV_SOF12_24",11);
	//Compare("NoTx","TN_TOL_G1a_PrOD_RBV12_24",11);

	//Compare("NoTx","TN_G1_LDP_SOF12",12);
	//Compare("NoTx","TN_TOL_G1b_DCV_SOF12_24",12);
	//Compare("NoTx","TN_TOL_G1b_PrOD12",12);
	//Compare("NoTx","TN_TOL_G1b_SOF_SMV12",12);

	Compare("NoTx", "TN_TOL_G2_SOF_RBV12_16", 2);


	Compare("NoTx", "G3_SOF_DCV12_24", 3);
	//Compare("NoTx","TN_TOL_G3_DCV_SOF12_24",3);
	//Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_RBV24",3);
	//Compare("TN_TOL_G3_PEG_RBV24","TN_TOL_G3_SOF_PEG_RBV12",3);

	Compare("NoTx", "G4_SOF_LDV12", 4);
	//Compare("NoTx","TN_TOL_G4_LDV_SOF12",4);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_PrOD_RBV12",4);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_RBV24",4);
	//Compare("TN_TOL_G4_PEG_RBV48","TN_TOL_G4_SOF_PEG_RBV12",4);

	return 0;
}

vector<double> Project_CEA_Global::GetAggregatedResults() {
	vector<double> outVec;
	vector<double> v1, v2, v3, v4;
	v1 = Compare("NoTx", "G1_SOF_LDV12", 1);
	v2 = Compare("NoTx", "TN_TOL_G2_SOF_RBV12_16", 2);
	v3 = Compare("NoTx", "G3_SOF_DCV12_24", 3);
	v4 = Compare("NoTx", "G4_SOF_LDV12", 4);
	for (int k = 0; k < v1.size(); k++) {
		outVec.push_back(_modelParam._pplData._distr_genotype["G1"] * v1[k]
			+ _modelParam._pplData._distr_genotype["G2"] * v2[k]
			+ _modelParam._pplData._distr_genotype["G3"] * v3[k]
			+ _modelParam._pplData._distr_genotype["G4"] * v4[k]);
	}	
	return outVec;
}

vector<double> Project_CEA_Global::Compare(string strArm1, string strArm2, int argGenotype)
{

	vector<double> outVec;
	_outFileSummary.open(_outFile, ios::app);
	_outFileSummary << fixed << showpoint;
	
	
	// Naming rule: arm_[XX]_[genotype]_[drug name(s)]_[treatment duration, weeks]
	/********************************** Define trial arms **************************************************************/
	// === scenarios ===
	int nScenarios = 10;
	stateType listState[] = { s_F0,s_F0,s_F1,s_F1,s_F2,s_F2,s_F3,s_F3,s_CoCirr,s_CoCirr };
	char listGender[] = { 'M','F','M','F','M','F','M','F','M','F' };

	double listAge[10];
	for (int k = 0; k < sizeof(listAge) / sizeof(double); k++) {
		listAge[k] = _age_initial;
	}

	clock_t t_start = clock();
	
	for (int k = 0; k < nScenarios; k++) {
		
		double wtScenario = _modelParam._pplData._distr_fib[listState[k]] * _modelParam._pplData._distr_gender[listGender[k]];

		cout << setprecision(0)<<fixed<<"[" << Time(t_start) << "\tsec]\t";
		cout << "Comparing\t" << strArm1 << "-" << strArm2 << ":\tFibrosis=" << listState[k] << "\tAge=" << listAge[k]
			<< "\tGender=" << listGender[k] << "\tGenotype=" << argGenotype << "\t" << endl;

		int genotype_noSubtype = (argGenotype == 11 || argGenotype == 12) ? 1 : argGenotype;

		// define cohort
		// create patient profile:
		// [state=s_F0, age=50, gender=male, race=white, genotype=1, il28b=none,priorTxRes=all]
		//baseCohortType testCohort(listState[k],listAge[k],listGender[k], 'W', argGenotype, 'N', 'A');	
		baseCohortType testCohort(listState[k], listAge[k], listGender[k], 'W', genotype_noSubtype, 'N', 'A');

		//vector<double> r1=EvaluateOneArm(strArm1, testCohort);
		//vector<double> r2=EvaluateOneArm(strArm2, testCohort);

		// changed, 7/9/2016
		//vector<double> r1=EvaluateOneArm_CombineLDV8Or12WkOnline(strArm1, testCohort);
		//vector<double> r2=EvaluateOneArm_CombineLDV8Or12WkOnline(strArm2, testCohort);
		vector<double> r1 = EvaluateOneArm_Averaged(strArm1, testCohort);
		vector<double> r2 = EvaluateOneArm_Averaged(strArm2, testCohort);

		assert(r1.size() == r2.size());
		_outFileSummary << _age_initial << "\t" << TIME_HORIZON << "\t"
			<< r1[0] << "\t" << r2[0] << "\t" << r2[0] - r1[0] << "\t"
			<< r1[1] << "\t" << r2[1] << "\t" << r2[1] - r1[1] << "\t"
			<< (r2[1] - r1[1]) / (r2[0] - r1[0]) << "\t";

		for (int l = 2; l < r1.size(); l++) {
			_outFileSummary << r1[l] << "\t" << r2[l] << "\t";
		}

		_outFileSummary << _modelParam._costData.c_SOF << "\t";
		_outFileSummary << endl;


		// initialize outVec
		if (k == 0) {
			for (int l = 0; l < 2 * r1.size(); l++) {
				outVec.push_back(0.0);
			}
		}

		for (int l = 0; l < r1.size(); l++) {
			outVec[l] = outVec[l] + wtScenario * r1[l];
			outVec[r1.size() + l] = outVec[r1.size() + l] + wtScenario * r2[l];
		}



	}

	_outFileSummary << endl;
	_outFileSummary.close();



	return outVec;
}


vector<double> Project_CEA_Global::EvaluateOneArm_Averaged(string txArm, const baseCohortType & testCohort)
{
	vector<double> r;
	if ("TN_TOL_G1_BOC_TEL" == txArm) {
		vector<double> r_boc = EvaluateOneArm("TN_TOL_G1_BOC", testCohort);
		vector<double> r_tel = EvaluateOneArm("TN_TOL_G1_TEL", testCohort);
		r = LinearCombTwoVec(r_boc, r_tel, _modelParam._transData._pBOC);

	}
	else if ("TE_G1_BOC_TEL" == txArm) {
		vector<double> r_boc = EvaluateOneArm("TE_G1_BOC", testCohort);
		vector<double> r_tel = EvaluateOneArm("TE_G1_TEL", testCohort);
		r = LinearCombTwoVec(r_boc, r_tel, _modelParam._transData._pBOC);

	}
	else if ("TN_G1_LDP_SOF" == txArm) {
		vector<double> r_8week = EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
		vector<double> r_12week = EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
		r = LinearCombTwoVec(r_8week, r_12week, _modelParam._transData._p8Week_LDV);
		//assert(r_8week.size()==r_12week.size());
		//for(int l=0; l<r_8week.size(); l++){
		//	r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1-_modelParam._transData._p8Week_LDV)*r_12week[l]);
		//}
	}
	else {
		r = EvaluateOneArm(txArm, testCohort);
	}
	return r;
}



vector<double> Project_CEA_Global::EvaluateOneArm(string strArm, const baseCohortType & testCohort)
{
	modelParamType arm = _modelParam;
	arm._armName = strArm;
	arm._cohortData = testCohort;
	arm.ReduceDrugCost(_drugCostReduction);
	//arm.ApplyIndiaParameters();

	HepCSim * mySim;
	mySim = new HepCSim;
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
	r.push_back(mySim->GetDALY_YLL());
	r.push_back(mySim->GetDALY_YLD());
	r.push_back(mySim->GetDALY_Total());
	r.push_back(mySim->GetTestCost());
	r.push_back(mySim->GetAvgLY());

	FreeMem(mySim);
	return r;

}

//vector<double> Project_CEA_Global::EvaluateOneArm_CombineLDV8Or12WkOnline(string txArm, const baseCohortType & testCohort )
//{
//
//	if("TN_G1_LDP_SOF" == txArm){
//		vector<double> r;
//		vector<double> r_8week=EvaluateOneArm("TN_G1_LDP_SOF8", testCohort);
//		vector<double> r_12week=EvaluateOneArm("TN_G1_LDP_SOF12", testCohort);
//		assert(r_8week.size()==r_12week.size());
//		for(int l=0; l<r_8week.size(); l++){
//			r.push_back(_modelParam._transData._p8Week_LDV*r_8week[l] + (1-_modelParam._transData._p8Week_LDV)*r_12week[l]);
//		}
//		return r;
//	}else{
//		return EvaluateOneArm(txArm, testCohort);
//	}
//
//}

//int Project_CEA_Global::ReadPSASampledValues( ifstream & inf, modelParamType & argModelParam )
//{
//	//update the model parameters:
//	//Fibrosis score (F0-F4)	sex	age	
//	// q_F0	q_F1	q_F2	q_F3	q_CoCirr	q_DeCirr	q_HCC	q_LivTr	q_PostLivT	q_SVR	q_Anemia	q_TX_oSOC	q_Tx_DAA	
//	inf>>argModelParam._qolData.q_F0;
//	inf>>argModelParam._qolData.q_F1;
//	inf>>argModelParam._qolData.q_F2;
//	inf>>argModelParam._qolData.q_F3;
//	inf>>argModelParam._qolData.q_CoCirr;
//	inf>>argModelParam._qolData.q_DeCirr;
//	inf>>argModelParam._qolData.q_HCC;
//	inf>>argModelParam._qolData.q_LivTr;
//	inf>>argModelParam._qolData.q_PostLivTr;
//	inf>>argModelParam._qolData.q_SVR;
//	inf>>argModelParam._qolData.q_Dec_Anemia;
//	inf>>argModelParam._qolData.q_TX_oSOC;
//	inf>>argModelParam._qolData.q_TX_DAA;
//
//	//c_F0	c_F1	c_F2	c_F3	c_CoCirr	c_DeCirr	c_DeCirr1yrPlus	c_HCC	c_LivTr	c_PostLivTr	
//	inf>>argModelParam._costData.c_F0;
//	inf>>argModelParam._costData.c_F1;
//	inf>>argModelParam._costData.c_F2;
//	inf>>argModelParam._costData.c_F3;
//	inf>>argModelParam._costData.c_CoCirr;
//	inf>>argModelParam._costData.c_DeCirr;
//	inf>>argModelParam._costData.c_DeCirr1yrPlus;
//	inf>>argModelParam._costData.c_HCC;
//	inf>>argModelParam._costData.c_LivTr;
//	inf>>argModelParam._costData.c_PostLivTr;
//
//	//pF0_F1_SA	pF1_F2_SA	pF2_F3_SA	pF3_F4_SA	pF4_DC_SA	pF4_HCC_SA	pDC_HCC_SA	pDC_Liv_Transpl_SA	
//	inf>>argModelParam._transData.pr_F0_F1;
//	inf>>argModelParam._transData.pr_F1_F2;
//	inf>>argModelParam._transData.pr_F2_F3;
//	inf>>argModelParam._transData.pr_F3_CoCirr;
//	inf>>argModelParam._transData.pr_CoCirr_DeCirr;
//	inf>>argModelParam._transData.pr_CoCirr_HCC;
//	inf>>argModelParam._transData.pr_DeCirr_HCC;
//	inf>>argModelParam._transData.pr_DeCirr_LivTr;
//
//	//pMort_dc_cyc_1_SA	pMort_dc_cyc_2_SA	pHCC_Liv_Transpl_SA	pMort_hcc_cyc_SA	pMort_Liv_Transpl_SA	pMort_Post_Liv_Transpl_SA	
//	inf>>argModelParam._transData.pr_DeCirr_DeathLiv;
//	inf>>argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv;
//	inf>>argModelParam._transData.pr_HCC_LivTr;
//	inf>>argModelParam._transData.pr_HCC_DeathLiv;
//	inf>>argModelParam._transData.pr_LivTr_DeathLiv;
//	inf>>argModelParam._transData.pr_LivTr1yrPlus_DeathLiv;
//
//	//pr_SVR_CoCirr_DeCirr	pr_SVR_CoCirr_HCC	pSVR_Delta_oSOC	pSVR_Delta_DAA
//	inf>>argModelParam._transData.pr_SVR_CoCirr_DeCirr;
//	inf>>argModelParam._transData.pr_SVR_CoCirr_HCC;
//	inf>>argModelParam._transData.pr_SVR_Delta_oSOC;
//	inf>>argModelParam._transData.pr_SVR_Delta_DAA;
//
//
//
//	return 0;
//}

int Project_CEA_Global::CEA_PSA(int argIdx, int argNBatches) {
	ofstream outf_psa;
	outf_psa.open("output_project_india/PSA_" + basicToStr(argIdx) + ".txt");
	outf_psa << fixed << showpoint;

	psaDistrType psaDistr;
	psaDistr.ReadDistrForPSA("project_india_intput_PSA_distr.txt");

	HepCSim mySim_psaSampler;		//used for sampling
	mySim_psaSampler.PSA_initializeSampler(psaDistr);

	int nIterOfThisBatch = NUM_RUNS_PSA / argNBatches;
	for (int n = 0; n < nIterOfThisBatch; n++) {
		if ((n + 1) % (nIterOfThisBatch / 10) == 0) {
			cout << (n + 1) << " PSA iterations finished...(" << (n + 1)*1.0 / (1.0 * nIterOfThisBatch) * 100 << "%)" << endl;
		}
		// sample
		mySim_psaSampler.PSA_sampleModelParamValue(_modelParam);	// sample cost/probability/quality value
																	// solve


		//if (argIdx == 0 && n == 0) {
		//	// print table head
		//	outf_psa << "iter\t"
		//		<< "QALY1\tCost1\t"
		//		<< "QALY2\tCost2\t"
		//		<< endl;
		//}

		vector<double> vecResult;
		outf_psa << argIdx * nIterOfThisBatch + n<<"\t";
		vecResult = GetAggregatedResults();
		for (int k = 0; k < vecResult.size(); k++) {
			outf_psa << vecResult[k] << "\t";
		}
		outf_psa << endl;



	}
	outf_psa.close();

	return 0;


}

//	ofstream outf_psa;
//	outf_psa.open("output_project_doublecheck/PSA_"+_listCmp[argCmpIdx]+"_"+basicToStr(argBatch)+".txt");
//	outf_psa<< fixed << showpoint;
//
//	ifstream inf;
//	inf.open("project_doublecheck_input_PSA.txt");
//	string line;
//	getline(inf,line);//skip the firstline
//	
//	//skip the first argBatch*1000 lines
//	for(int i=0; i<argBatchSize*argBatch; i++) getline(inf,line);
//
//	int sa_counter=1+argBatchSize*argBatch;
//
//	int fib;
//	char gender;
//	double initialAge;
//	for(int i=0; i<argBatchSize; i++){
//		inf>>fib;
//		inf>>gender;
//		inf>>initialAge;
//		// create cohort
//		baseCohortType testCohort(fib,initialAge,gender, 'W', _listGenotype[argCmpIdx], 'N', 'A');	// first three charateristics (fib, age, gender) will be updated in the next line
//		// read sampled parameter values
//		ReadPSASampledValues(inf,_modelParam);
//	
//
//		vector<double> r1=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm1[argCmpIdx], testCohort);
//		vector<double> r2=EvaluateOneArm_CombineLDV8Or12WkOnline(_listArm2[argCmpIdx], testCohort);
//
//		assert(r1.size()==r2.size());
//		outf_psa<<"F"<<fib<<" "<<gender<<" "<<initialAge<<"\t"<<sa_counter<<"\t"
//			<<r1[0]<<"\t"<<r2[0]<<"\t"<<r2[0]-r1[0]<<"\t"
//			<<r1[1]<<"\t"<<r2[1]<<"\t"<<r2[1]-r1[1]<<"\t"
//			<<(r2[1]-r1[1])/(r2[0]-r1[0])<<"\t";		
//		for(int i=2;i<r1.size(); i++){
//			outf_psa<<r1[i]<<"\t"<<r2[i]<<"\t";
//		}
//		outf_psa<<endl;
//
//		sa_counter++;
//	}
//
//	inf.close();
//	outf_psa.close();
//
//	return 0;
//}

int Project_CEA_Global::ChangeValue(string varName, double argVal, modelParamType & argModelParam) {
	if ("q_F0" == varName) {
		argModelParam._qolData.q_F0 = argVal;
	}
	else if ("q_F1" == varName) {
		argModelParam._qolData.q_F1 = argVal;
	}
	else if ("q_F2" == varName) {
		argModelParam._qolData.q_F2 = argVal;
	}
	else if ("q_F3" == varName) {
		argModelParam._qolData.q_F3 = argVal;
	}
	else if ("q_CoCirr" == varName) {
		argModelParam._qolData.q_CoCirr = argVal;
	}
	else if ("q_DeCirr" == varName) {
		argModelParam._qolData.q_DeCirr = argVal;
	}
	else if ("q_HCC" == varName) {
		argModelParam._qolData.q_HCC = argVal;
	}
	else if ("q_LivTr" == varName) {
		argModelParam._qolData.q_LivTr = argVal;
	}
	else if ("q_PostLivT" == varName) {
		argModelParam._qolData.q_PostLivTr = argVal;
	}
	else if ("q_SVR" == varName) {
		argModelParam._qolData.q_SVR = argVal;
	}
	else if ("q_Anemia" == varName) {
		argModelParam._qolData.q_Dec_Anemia = argVal;
	}
	else if ("q_TX_oSOC" == varName) {
		argModelParam._qolData.q_TX_oSOC = argVal;
	}
	else if ("q_Tx_DAA" == varName) {
		argModelParam._qolData.q_TX_DAA = argVal;
	}
	else if ("c_F0" == varName) {
		argModelParam._costData.c_F0 = argVal;
	}
	else if ("c_F1" == varName) {
		argModelParam._costData.c_F1 = argVal;
	}
	else if ("c_F2" == varName) {
		argModelParam._costData.c_F2 = argVal;
	}
	else if ("c_F3" == varName) {
		argModelParam._costData.c_F3 = argVal;
	}
	else if ("c_CoCirr" == varName || "c_F4" == varName) {
		argModelParam._costData.c_CoCirr = argVal;
	}
	else if ("c_DC" == varName) {
		argModelParam._costData.c_DeCirr = argVal;
		argModelParam._costData.c_DeCirr1yrPlus = argVal;
	}

	else if ("c_DeCirr" == varName) {
		argModelParam._costData.c_DeCirr = argVal;
	}
	else if ("c_DeCirr1yrPlus" == varName) {
		argModelParam._costData.c_DeCirr1yrPlus = argVal;
	}
	else if ("c_HCC" == varName) {
		argModelParam._costData.c_HCC = argVal;
	}
	else if ("c_LivTr" == varName) {
		argModelParam._costData.c_LivTr = argVal;
	}
	else if ("c_PostLivTr" == varName) {
		argModelParam._costData.c_PostLivTr = argVal;
	}
	else if ("pF0_F1_SA" == varName) {
		argModelParam._transData.pr_F0_F1 = argVal;
	}
	else if ("pF1_F2_SA" == varName) {
		argModelParam._transData.pr_F1_F2 = argVal;
	}
	else if ("pF2_F3_SA" == varName) {
		argModelParam._transData.pr_F2_F3 = argVal;
	}
	else if ("pF3_F4_SA" == varName) {
		argModelParam._transData.pr_F3_CoCirr = argVal;
	}
	else if ("pF4_DC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_DeCirr = argVal;
	}
	else if ("pF4_HCC_SA" == varName) {
		argModelParam._transData.pr_CoCirr_HCC = argVal;
	}
	else if ("pDC_HCC_SA" == varName) {
		argModelParam._transData.pr_DeCirr_HCC = argVal;
	}
	else if ("pDC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_DeCirr_LivTr = argVal;
	}
	else if ("pMort_dc_cyc_1_SA" == varName) {
		argModelParam._transData.pr_DeCirr_DeathLiv = argVal;
	}
	else if ("pMort_dc_cyc_2_SA" == varName) {
		argModelParam._transData.pr_DeCirr1yrPlus_DeathLiv = argVal;
	}
	else if ("pHCC_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_HCC_LivTr = argVal;
	}
	else if ("pMort_hcc_cyc_SA" == varName) {
		argModelParam._transData.pr_HCC_DeathLiv = argVal;
	}
	else if ("pMort_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr_DeathLiv = argVal;
	}
	else if ("pMort_Post_Liv_Transpl_SA" == varName) {
		argModelParam._transData.pr_LivTr1yrPlus_DeathLiv = argVal;
	}
	else if ("pr_SVR_CoCirr_DeCirr" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_DeCirr = argVal;
	}
	else if ("pr_SVR_CoCirr_HCC" == varName) {
		argModelParam._transData.pr_SVR_CoCirr_HCC = argVal;
	}
	else if ("pSVR_Delta_oSOC" == varName) {
		argModelParam._transData.pr_SVR_Delta_oSOC = argVal;
	}
	else if ("pSVR_Delta_DAA" == varName) {
		argModelParam._transData.pr_SVR_Delta_DAA = argVal;


	}

	else if ("dw_F0" == varName) {
		argModelParam._dalyData._dw_f0 = argVal;
	}
	else if ("dw_F1" == varName) {
		argModelParam._dalyData._dw_f1 = argVal;
	}
	else if ("dw_F2" == varName) {
		argModelParam._dalyData._dw_f2 = argVal;
	}
	else if ("dw_F3" == varName) {
		argModelParam._dalyData._dw_f3 = argVal;
	}
	else if ("dw_CoCirr" == varName) {
		argModelParam._dalyData._dw_CoCirr = argVal;
	}
	else if ("dw_DeCirr" == varName) {
		argModelParam._dalyData._dw_DeCirr = argVal;
		argModelParam._dalyData._dw_DeCirr1yrPlus = argVal;
	}
	else if ("dw_HCC" == varName) {
		argModelParam._dalyData._dw_HCC = argVal;
	}
	else if ("dw_LivTr" == varName) {
		argModelParam._dalyData._dw_LivTr = argVal;
	}
	else if ("dw_PostLivT" == varName) {
		argModelParam._dalyData._dw_LivTr1yrPlus = argVal;
	}

	else if ("c_testing_preTx" == varName || "c_preTxTesting" == varName) {
		argModelParam._costData.c_testing_preTx = argVal;
	}
	else if ("c_testing_postTx" == varName || "c_postTxTesting" == varName) {
		argModelParam._costData.c_testing_postTx = argVal;
	}
	else if ("c_DAA_wk" == varName) {
		argModelParam._costData.c_SOF = argVal;
	}
	else {
		ExitWithMsg("[Error] Project_Doublecheck::ChangeValue(string varName, double argVal): Unknown parameter " + varName);
	}
	return 0;
}



int Project_CEA_Global::CEA_VaryAge(int argIdx)
{

	//double listInitialAge[] = { 20,25,30,35,40,45,50,55,60,65,70 };
	//for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {

	//	SetInitialAge(listInitialAge[k]);

	//	CEA_BaseResults();
	//}

	//return 0;


	ofstream outf;
	outf.open("output_project_india/out_change_age.txt");
	outf << fixed << showpoint;

	double listInitialAge[] = { 20,25,30,35,40,45,50,55,60,65,70 };


	for (int k = 0; k < sizeof(listInitialAge) / sizeof(double); k++) {

		SetInitialAge(listInitialAge[k]);
		vector<double> rs = GetAggregatedResults();



		outf << listInitialAge[k] << "\t" << TIME_HORIZON << "\t";
		for (int l = 0; l < rs.size(); l++) { outf << rs[l] << "\t"; }
		outf << endl;

	}

	outf.close();
	return 0;
}

