// StatisticUtil.h
// Description:
// 1. Some functions to obtain some useful statistics from a set of data
// 2. Integrate the boost library for "statistics distribution"
// Qiushi Chen, April 1, 2010.

#ifndef STATUTIL_H
#define STATUTIL_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;


// calculate the mean and variance simultaneously
// valid only for PRIMATIVE TYPE (int, double, long, primative)
template <typename T>
inline void MeanAndVariance(const vector<T> & dataArray, double & mean, double & variance){
	if(dataArray.size()==1){
		//cout<<"[Warning][MeanAndVariance]: The data array has only one element, set variance = bigM ("<<big_M<<")"<<endl;
		mean=dataArray[0];
		variance=big_M;
		return;
	}

	double oldmean;
	for(int i=0;i<dataArray.size();i++){
		if(0==i){
			mean=static_cast<double>(dataArray[i]);
			oldmean=mean;
			variance=0.0;
		}else{
			mean=mean*i/(i+1);
			mean=mean+static_cast<double>(dataArray[i])/(i+1);
			variance=variance*(i-1)/i;
			variance=variance+(mean-oldmean)*(mean-oldmean)*(i+1);
			oldmean=mean;
		}
	}
};

template<typename T>
inline int Min(const vector<T> & theVec, T& theMin){
	if(0==theVec.size()){
		cout<<endl<<"[Warning][Min()]: null vector... returning big_M as the min value"<<endl;
		theMin=big_M;
		return 0;
	}
	else
		theMin=theVec[0];
	
	int idx=0;
	
	if(1==theVec.size()){
		return idx;
	}
	
	//for (int i=1;i<theVec.size();i++){
	//	if( abs(theVec[i]-theMin)<1e-8){
	//		cout<<setprecision(9)<<"-- ties in vector: ["<<idx<<"]="<<theMin<<" and ["<<i<<"]="<<theVec[i]<<endl;
	//	}else if(theMin - theVec[i] > EPSILON){
	//		theMin=theVec[i];
	//		idx=i;
	//	}
	//}

	for (int i=0;i<theVec.size();i++){
		if (theMin>theVec[i]){
			theMin=theVec[i];
			idx=i;
		}
	}
	return idx;
};

template<typename T>
inline int Max(const vector<T> & theVec, T& theMax){
	if(0==theVec.size()){
		cout<<endl<<"[Warning][Max()]: null vector... returning big_M as the max value"<<endl;
		theMax=big_M;
		return 0;
	}
	else
		theMax=theVec[0];
	int idx=0;
	for (int i=0;i<theVec.size();i++){
		if (theMax<theVec[i]){
			theMax=theVec[i];
			idx=i;
		}
	}
	return idx;
};

// return multiple maximum.
template<typename T>
inline vector<int> GetMaxIndicies(const vector<T> & theVec){
	vector<int> indicies;
	if(theVec.size()>0){

		T maxval=theVec[0];
		indicies.push_back(0);
		for(int i=1;i<theVec.size();i++){
			if(theVec[i]-maxval>EPSILON){
				maxval=theVec[i];
				indicies.clear();
				indicies.push_back(i);
			}else if(abs(theVec[i]-maxval)<=EPSILON){
				indicies.push_back(i);
			}else{

			}
			
		}
	}
	return indicies;
}


#endif