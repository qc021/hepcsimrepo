//: UtilPackage.h
// 	Declarations and definitions of commonly used utility macros,
//	functions, and templates, etc.
// 	Author:	Qiushi Chen @ IE Dept., Tsinghua University.
//  Version: 1.0
//	Updated: Feb 4, 2010, Qiushi Chen

#ifndef __UTIL_H__
#define __UTIL_H__

//	Standard C++ Headers
#include <cassert>
#include <ctime> 	// for clock_t, CLOCKS_PER_SEC.
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>

#include <iomanip>
#include <map>
#include <iterator>
#include <cmath>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <sstream>
//#include <boost/iostreams/stream.hpp>
//#include <boost/iostreams/tee.hpp>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/variate_generator.hpp>

using namespace std;


// global variables:
// by Qiushi Chen, March, 2010

//typedef boost::iostreams::tee_device<std::ostream, std::ofstream> type_dual_device;
//typedef boost::iostreams::stream<type_dual_device> type_dual_out_stream;

#define EPSILON 1e-6
#define BIGM 1e20

//	Constants.
const double util_zero_tol = 1.0e-6;
const double big_M = 1.0e+20;

// =====================================================================
// Macros.
// From: Lei, Zhao @ IE Dept., Tsinghua University.

#ifdef DEBUG
#define _ERROR(msg) { cout << "ERROR: " << msg << endl; assert(0); }
#else
#define _ERROR(msg) {}
#endif

#ifdef DEBUG
#define _CASSERT(x, y) { if (!(x)) _ERROR((y)); }
#else
#define _CASSERT(x, y) {}
#endif

#ifdef DEBUG
#define _PRTVAR(var) { cout << #var " = " << (var) << endl; }
#else
#define _PRTVAR(var) { cout << #var " = " << (var) << endl; }
#endif




//	Inline funcitons.
//
// Traces the function visited in the program.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
inline void trace_func(const char * func_name)
{
    cout << "in function: " << func_name << "() ... " << endl;
}
//
// A simple function to mark the successful completion of a program.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
inline void prt_done(string str = "Congratulations!!! DONE!!!")
{
    int extra = 8;
    int	len = (int) str.size();

    for (int i = 0; i < len + extra; i++)
        cout << "~";
    cout << "\n^_^ " << str << " @_@\n";
    for (int i = 0; i < len + extra; i++)
        cout << "~";
    cout << "\n" << endl;
}

//	A simple funtion to mark two part of the output.
//	Author: Peng Xiaoshan
//
inline void prt_break(char mark='*')
{
    int length=30;
    for (int i = 0; i < length; i++)
        cout << mark ;
    cout <<endl;
}

//
// Calculates CPU time in seconds.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
inline double calc_time(clock_t start_time, clock_t end_time)
{
    return static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;
}


//
// Prints out the elements stored in a container.
// From: Lei, Zhao
template<typename CPP_Type>
inline void prt_cont(CPP_Type *vec, int len, const string name)
{
    // Five elements per line.
    int mod = 5;
    int	cnt = 0;

    cout << "Printing container ... " << name << endl;
    typedef typename CPP_Type::iterator itertype;
    itertype iter = (*vec).begin();
    for (; iter != (*vec).end(); iter++)
    {
        cout << *iter << "\t";
        if ((cnt % mod) == (mod - 1))
            cout << endl;
        cnt++;
    }

    // Avoids duplicated 'return' at the end.
    if ((cnt % mod) != 0)
        cout << endl;
}



//
// This function appends an arry, _append_vec_, to another array,
// _orig_vec_.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline void append_vect(CPP_Type **orig_vec, CPP_Type **append_vec,
                        int orig_sz, int append_sz)
{
    inflate_vect(orig_vec, orig_sz, orig_sz + append_sz);
    for (int i = 0; i < append_sz; i++)
        (*orig_vec)[orig_sz + i] = (*append_vec)[i];
}

//
// This function compares if the two (int or char) vectors are equal. It
// returns true if all the element pairs of the two vectors are equal.
// Note:
// The starting position of the elements under comparison of
// _v1_ is specified by _beg1_, while that of _v2_ is specified
// by _beg1 + shift2_. That is, we are comparing [beg1, beg1 + len]
// of _v1_ and [beg1 + shift2, beg1 + shift2 + len] of _v2_.
// The default of is to compare [0, len] of _v1_ and _v2_.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline bool cmp_icvects(const CPP_Type *v1, const CPP_Type *v2, int len,
                        int beg1 = 0, int shift2 = 0)
{
    _CASSERT(len >= 0, "Length must be nonnegative");
    _CASSERT(beg1 >= 0, "Nonnegative index needed for v1");

    int beg2 = beg1 + shift2;
    _CASSERT(beg2 >= 0, "Nonnegative index needed for v2");

    if (0 == len)
        return true;

    for (int i = 0; i < len; i++)
    {
        if (v1[beg1 + i] != v2[beg2 + i])
            return false;
    }

    return true;
}

//
// This function copys an array, _copy_vec_, from another array,
// _orig_vec_.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline void copy_vect(CPP_Type **copy_vec, const CPP_Type *orig_vec,
                      int orig_sz)
{
    for (int i = 0; i < orig_sz; i++)
        (*copy_vec)[i] = orig_vec[i];
}

//
// This function increase the size of an array from the old
// size, _old_sz_, to a new size, _new_sz_, with the original
// elements copied to the new array.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline void inflate_vect(CPP_Type **vec, int old_sz, int new_sz)
{
    CPP_Type *tmp_vec = new CPP_Type[new_sz];
    memset(tmp_vec, 0, new_sz * sizeof(CPP_Type));
    if (*vec)
    {
        memcpy(tmp_vec, *vec, old_sz * sizeof(CPP_Type));
        delete [] *vec;
    }
    (*vec) = tmp_vec;
}

//
// A function to print a vector of a given length and name.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline void prt_vect(const CPP_Type *vec, int len, const char *vec_name,
                     ostream & out = cout)
{
    // Five elements per line.
    int mod = 10;

    out << "value of " << vec_name << "are as follows:"<<endl<<"\t";
    for (int i = 0; i < len; i++)
    {
        out << vec[i] << "\t";
        if ((i % mod) == (mod - 1))
            out << endl;
    }

    // Avoids duplicated 'return' at the end.
    if ((len % mod) != 0)
        out << endl;
}

//
// These function loops through the given vector (NOT an array) of pointers,
// frees them, and removes them from the vector.
// NOTE:
// 		The given vector was pointed to by the pointer _vect_p_.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline void free_ptrvect(vector<CPP_Type *> *vect_p)
{
#ifdef TRACE
    trace_func("free_ptrvect");
#endif

    typedef typename vector<CPP_Type *>::iterator vectorIterType;
    vectorIterType iter;

    //	for (iter = vect.begin(); iter != vect.end(); iter) {
    for (iter = vect_p->begin(); iter != vect_p->end(); iter)
    {
        if (*iter)
        {
            //			cout << (*(*iter)).get_id() << endl;
            // Frees memory that the pointer pointing to.
            delete *iter;
            (*iter) = NULL;
        }
        // Removes the element from the vector.
        //		vect.erase(iter);
        vect_p->erase(iter);
    }

    //	_PRTVAR(vect.size());
}

//
// These two functions deal with dynamic memory allocations.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline void mem_malloc(CPP_Type **ppt, int length)
{
    (*ppt) = new CPP_Type[length];
}

template<class CPP_Type>
inline void mem_free(CPP_Type **ppt)
{
    if (*ppt)
        delete [] (*ppt);
    (*ppt) = NULL;
}

// Function prototypes.

bool cmp_dvects(const double *v1, const double *v2, int len, int beg1 = 0,
                int shift2 = 0);

// =====================================================================

//
//Given the time point of start running the program,
//return the current "simulation clock" in seconds
//type clock_t has the unit of "ms".
//
// By Qiushi Chen, Feb 5, 2010
inline double Time(clock_t start)
{
    return (double)(clock()-start)/CLOCKS_PER_SEC;
}

//
// This function converts a basic type to a string.
// From: Lei, Zhao @ IE Dept., Tsinghua University.
// Origin: util.h
template<class CPP_Type>
inline std::string basicToStr(CPP_Type basic)
{
    using namespace std;
    ostringstream stream;
    stream << basic << flush;
    string str(stream.str());

    return str;
}



// free the memory of a pointer
template<typename ptrType>
inline int FreeMem(ptrType & p)
{
    if(NULL!=p)
    {
        delete p;
        p=NULL;
        return 0;
    }
    return 1;
}

template<typename pointerVectorType>
inline int FreeMemPtrVec(pointerVectorType & ptrVec)
{
    for(int i=0; i<ptrVec.size(); i++)
    {
        FreeMem(ptrVec[i]);
    }
    return 0;
};


// SHow the error message, and exit the program
inline int ExitWithMsg(string errMsg,int flag=1)
{
    cout<<endl<<errMsg<<endl;
    getchar();
    exit(flag);
}

template<typename T>
inline int OutputVector(const vector<T> & theVec,string fileName)
{
    ofstream outf(fileName.c_str());
    if(outf.fail())
    {
        ExitWithMsg("Error! [OutputVector()]: Can't open the file..."+fileName);
    }
    outf<<setiosflags(ios::fixed)<<setprecision(3);
    for(int i=0; i<(int)theVec.size(); i++)
    {
        outf<<i<<"\t"<<theVec[i]<<endl;
    }
    outf.close();
    return 0;
}

//
template<typename T1, typename T2>
inline void FindSameKeyAndLocateIter(typename map<const T1* ,T2>::iterator & theIter, const map<const T1* ,T2> & theMap, T1*  pTarget)
{
    for(; theIter!=theMap.end(); theIter++)
    {
        if((*(theIter->first))==(*pTarget))
        {
            // find the "same" region
            break;
        }
    }
}

template<typename T>
inline void FindSameKeyAndLocateIter(typename vector<T*>::const_iterator & theIter,  const vector<T*> & theVec, T*  pTarget)
{
    for(; theIter!=theVec.end(); theIter++)
    {
        if((**theIter)==(*pTarget))
        {
            // find the "same" region
            break;
        }
    }
}

template<typename T>
inline void PrintTable(vector<vector<T> > temptable, int pres=6)
{
    //print the table
    cout<<fixed<<setprecision(pres);
    for(vector<vector<double> >::iterator iter1=temptable.begin(); iter1!=temptable.end(); iter1++)
    {
        for(vector<double>::iterator iter2=iter1->begin(); iter2!=iter1->end(); iter2++)
        {
            cout<<*iter2<<"\t";
        }
        cout<<endl;
    }
}

template<typename T>
inline string VectorToString(const vector<T> & v,string sep="")
{
    string s="[\t";

    for(typename vector<T>::const_iterator iter=v.begin(); iter!=v.end(); iter++)
    {
        s+=(basicToStr(*iter))+sep+"\t";

    }
    s+="]";
    return s;
}

template<typename T>
inline string VectorToStringNoSep(const vector<T> & v)
{
    string s="[\t";
    for(typename vector<T>::const_iterator iter=v.begin(); iter!=v.end(); iter++)
    {
        s+=(basicToStr(*iter))+"\t";

    }
    s+="]";
    return s;
}

template<typename T>
inline bool AreSameVector(const vector<T> & v1, const vector<T> & v2)
{
    if(v1.size()!=v2.size())
        ExitWithMsg("Inconsistent size of two vectors...");

    if(v1.size()==0)
        ExitWithMsg("Empty vector...");

    for(int n=0; n<v1.size(); n++)
    {
        //if(v1[n]!=v2[n]){
        if(abs(v1[n]-v2[n])>EPSILON)
        {
            return false;
        }
    }
    return true;
}

template<typename T>
inline bool Has(const vector<T> & theVec, const T & theElem)
{
    for(typename vector<T>::const_iterator i=theVec.begin(); i!=theVec.end(); i++)
    {
        //if(*i == theElem) return true;
        if(abs(*i-theElem)<EPSILON) return true;
    }
    return false;
}

template<typename T>
inline bool Has(const vector<vector<T> > & theVec,const vector<T> & theElem)
{
    for(typename vector<vector<T> >::const_iterator i=theVec.begin(); i!=theVec.end(); i++)
    {
        if(AreSameVector(*i,theElem)) return true;
    }
    return false;
}

inline bool IsZeroVec(const vector<double> & theVec)
{
    for(vector<double>::const_iterator it=theVec.begin();
            it!=theVec.end(); it++)
    {
        double x=abs(*it - 0);
        if(x > EPSILON)
            return false;
    }
    return true;
}

inline int CountZeros(const vector<double> & b)
{

    int count=0;
    for(int i=0; i<b.size(); i++)
    {
        if(abs(b[i]-0)<=EPSILON)
        {
            count++;
        }
    }
    return count;
}


inline int ReadPairsReverse(map<double,double> & dist, ifstream & inf)
{

    string stra;
    double num1,num2;
    inf>>stra;
    while("END"!=stra)
    {
        num1=atof(stra.c_str());
        inf>>num2;
        dist.insert(pair<double,double>(num2,num1));
        inf>>stra;
    }
    return 0;
}

inline int ReadPairs(map<double,double> & dist, ifstream & inf)
{

    string stra;
    double num1,num2;
    inf>>stra;
    while("END"!=stra)
    {
        num1=atof(stra.c_str());
        inf>>num2;
        dist.insert(pair<double,double>(num1,num2));
        inf>>stra;
    }
    return 0;
}

template<typename T>
inline int ReadList(vector<T> & theList, ifstream & inf)
{
    theList.clear();
    T stra;
    
    while(inf>>stra)
    {
        theList.push_back(stra);
    }
    return 0;
}




// begIdx, endIdx are counted from 1. (not 0).
template<typename T>
inline void TruncVec(const vector<T> & orgVec, vector<T> & tgtVec, int begIdx, int endIdx )
{
    // include begIdx, not include endIdx
    tgtVec.clear();
    for(int i=begIdx-1; i<endIdx; i++)
    {
        tgtVec.push_back(orgVec[i]);
    }
}

template<typename T>
inline void CartesianProduct(const vector<vector<T> > & vecA, const vector<T> & vecB, vector<vector<T> > & vecResult)
{
    // include begIdx, not include endIdx
    for(typename vector<vector<T> >::const_iterator it=vecA.begin(); it!=vecA.end(); it++)
    {
        for(typename vector<T>::const_iterator i=vecB.begin(); i!=vecB.end(); i++)
        {
            vector<T> tempVec=*it;
            tempVec.push_back(*i);
            vecResult.push_back(tempVec);
        }
    }
}

template<typename T>
inline void CartesianProduct(const vector<T> & vecA, const vector<T> & vecB, vector<vector<T> > & vecResult, bool appending=true)
{
    // include begIdx, not include endIdx
    if(!appending)
    {
        vecResult.clear();
    }

    for(typename vector<T>::const_iterator i = vecB.begin(); i != vecB.end(); i++)
    {
        vector<T> tempVec = vecA;
        tempVec.push_back(*i);
        vecResult.push_back(tempVec);
    }

}

inline void GenSequence(double lb, double ub, double intv, vector<double> & seq, double ref=0)
{
    double val=max(lb+ref,0.0);
    seq.clear();
    while(val<=min(ub+ref,1.0))
    {
        seq.push_back(val);
        val=val+intv;
    }

}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
inline const std::string currentDateTime()
{
    time_t     now = time(0);
	
	//// [[[ alternative version ]]]
	//tm *ltm = localtime(&now);	
	//return  basicToStr(1900 + ltm->tm_year)+"-"
	//	+ basicToStr(1 + ltm->tm_mon) + "-"
	//	+ basicToStr(ltm->tm_mday) + ". "
	//	+ basicToStr(1 + ltm->tm_hour)+":"
	//	+ basicToStr(1 + ltm->tm_min)+":"
	//	+ basicToStr(1 + ltm->tm_sec);

	struct tm  tstruct;
	tstruct = *localtime(&now);
	//localtime_s(&tstruct, &now);
    char       buf[80];
    
    // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	//time_t _tm = time(NULL);
	//struct tm * curtime = localtime(&_tm);
	//cout << "The current date/time is:" << asctime(curtime);

    return buf;
}

inline vector<vector<vector<double> > > ZeroVec3d(int d1,int d2,int d3)
{
    vector<double> v3;
    for(int i=0; i<d3; i++) v3.push_back(0);
    vector<vector<double> > v2;
    for(int i=0; i<d2; i++) v2.push_back(v3);
    vector<vector<vector<double> > > v1;
    for(int i=0; i<d1; i++) v1.push_back(v2);
    return v1;
}

inline vector<vector<double> > ZeroVec2d(int d1,int d2)
{
    vector<double> v2;
    for(int i=0; i<d2; i++) v2.push_back(0);
    vector<vector<double> > v1;
    for(int i=0; i<d1; i++) v1.push_back(v2);
    return v1;
}

inline vector<vector<int> > ZeroVec2dInt(int d1,int d2)
{
    vector<int> v2;
    for(int i=0; i<d2; i++) v2.push_back(0);
    vector<vector<int> > v1;
    for(int i=0; i<d1; i++) v1.push_back(v2);
    return v1;
}
inline vector<double> ZeroVec1d(int d1)
{
    vector<double> v1;
    for(int i=0; i<d1; i++) v1.push_back(0);
    return v1;
}

inline vector<double> OneVec1d(int d1)
{
    vector<double> v1;
    for(int i=0; i<d1; i++) v1.push_back(1);
    return v1;
}

template<typename T>
inline int OutputMtx(const vector<T> & theMtx,string fileName)
{
    ofstream outf(fileName.c_str());

    if(outf.fail())
    {
        ExitWithMsg("Error! [OutputMatrx()]: Can't open the file..."+fileName);
    }

    for(int i=0; i<(int)theMtx.size(); i++)
    {
        for(int j=0; j<(int)theMtx[i].size(); j++)
        {
            outf<<theMtx[i][j]<<"\t";
        }
        outf<<endl;
    }

    outf.close();
    return 0;
}

// Added: 7/3/2014




template<typename T1, typename T2>
inline int ReadTable( string argFile, map<T1,T2> & argMap, bool argHeader=false )
{
	argMap.clear();
    cout<<"- Reading global parameters from: "<<argFile<<" ..."<<endl;
    ifstream inf;
    inf.open(argFile.c_str());
    if(inf.fail())      ExitWithMsg("CMarkovModel::ReadTable(): Fail to read "+argFile);

    string header;
    if (argHeader) getline(inf,header);

    T1 word1;
    T2 word2;
    while (inf >> word1) {
        inf>>word2;
        argMap.insert(pair<T1,T2>(word1,word2));
    }
    inf.close();
    return 0;
}


template<typename T1, typename T2, typename T3>
inline int ReadTable( string argFile, map<T1,map<T2,T3>> & argMap)
{
	argMap.clear();
    cout<<"- Reading global parameters from: "<<argFile<<" ..."<<endl;
    ifstream inf;
    inf.open(argFile.c_str());
    if(inf.fail())      ExitWithMsg("ReadTable(): Fail to read "+argFile);

    string word;
    inf>>word; // skip the first word
    // read the first line (table header)
    vector<T1> colName;
    T1 next;
    while(inf>>next){
    if (next == '\n'){  // If the file has been opened in
        break;
    }
    colName.push_back(next);
    map<T2,T3> newMap;
    argMap[next]=newMap;
    }


     T2 rowlabel;
     T3 elem;

    while(inf>>rowlabel){
        for(int colIdx=0; colIdx<colName.size(); colIdx++){
                inf>>elem;
                argMap.find(colName[colIdx])->second.insert(pair<T2,T3>(rowlabel,elem));
        }


    }

    return 0;
}


template<typename T1>
inline int ReadTable( string argFile, map<string,vector<T1>> & argMap)
{
	argMap.clear();
    cout<<"- Reading global parameters from: "<<argFile<<" ..."<<endl;
    ifstream inf;
    inf.open(argFile.c_str());
    if(inf.fail())      ExitWithMsg("ReadTable(): Fail to read "+argFile);
	// read the first line (table header)
	string line;
	getline(inf,line);
	vector<string> colName;
	istringstream iss(line);
	string sub;
	while(iss>>sub){
		if(""==sub) 
			break;
		colName.push_back(sub);
		vector<T1> newVec;
		argMap.insert(pair<string,vector<T1>>(sub,newVec));
	}
    
    T1 val;
    while(inf>>val){

        for(int colIdx=0; colIdx<colName.size(); colIdx++){
			if (colIdx>0) 
				inf>>val;
			argMap.find(colName[colIdx])->second.push_back(val);         
        }
    }
	inf.close();
    return 0;
}


template<typename T1>
inline int ReadTable(string argFile, map<int, vector<T1>> & argMap)
{
	argMap.clear();
	cout << "- Reading global parameters from: " << argFile << " ..." << endl;
	ifstream inf;
	inf.open(argFile.c_str());
	if (inf.fail())      ExitWithMsg("ReadTable(): Fail to read " + argFile);
	// read the first line (table header)
	string line;
	while (getline(inf, line)) {
		istringstream iss(line);
		int ky;
		iss >> ky;
		vector<T1> newVec;
		T1 sub;
		while (iss >> sub) {
/*			if ("" == sub)
				break;*/			
			newVec.push_back(sub);
		}
		argMap.insert(pair<int, vector<T1>>(ky, newVec));
	}
	inf.close();
	return 0;
}

template<typename T1, typename T2>
inline int ReadTable_2Col(string argFile, map<T1, T2> & argMap, bool argWithHeader=false)
{
	argMap.clear();
	ifstream inf;
	inf.open(argFile.c_str());
	if (inf.fail())      ExitWithMsg("ReadTable_2Col(): Fail to read " + argFile);
	// read the first line (table header)
	if (argWithHeader) {
		string line;
		getline(inf, line);
	}
	
	T1 idx;
	T2 val;
	while (inf >> idx) {
		inf >> val;
		pair<T1,T2> y(idx,val);
		argMap.insert(y);		
	}
	inf.close();
	return 0;
}


inline vector<int> ReadRowAsIntVector(ifstream & inFile){
	string t;
	getline(inFile, t);
	istringstream iss(t);
	int word;
	vector<int> theVec;
	while (iss >> word) {		
		theVec.push_back(word);
	}
	return theVec;
}


inline vector<double> ReadRowAsDoubleVector(ifstream & inFile){
	string t;
	getline(inFile, t);
	istringstream iss(t);
	double word;
	vector<double> theVec;
	while (iss >> word) {		
		theVec.push_back(word);
	}
	return theVec;
}


template<typename T>
inline int ReadMatrix( string argFile, vector<vector<T>> & argMtx, bool argHeader=true )
{
	argMtx.clear();
    ifstream inf;
    inf.open(argFile.c_str());
    if(inf.fail())      ExitWithMsg("ReadMatrix(): Fail to read "+argFile);

    string header;
    if (argHeader) getline(inf,header);

    T next;
    vector<T> row;
    while(inf>>next){
    if (next == '\n'){  // If the file has been opened in
        argMtx.push_back(row);
        row.clear();
    }
    row.push_back(next);
    }
    inf.close();
    return 0;
}

inline double Round(double argV, int argDigits){
    double coef=pow(10.0,argDigits);
    return floor(argV*coef+0.5)/coef;
}

inline double Interpolate (double x, const std::map<int, double> &table)
{
    assert(table.size() > 0);

    if(x<table.begin()->first){
        return table.begin()->second;        // lower bound
    }else if(x>table.rbegin()->first){
        return table.rbegin()->second;       // upper bound
    }else{
        // interpolate
        std::map<int, double>::const_iterator it = table.lower_bound((int)x); 
        double x1 = it->first;
		double y1 = it->second; 
		if (it == table.begin()) return y1;

		
        it--;
        double x2 = it->first;
        double y2 = it->second;
        double p = (x - x1) / (x2 - x1);		
        
		//return y2; // ONLY for DALY testing (compre with WHO spresheet).
		return (1 - p) * y1 + p * y2;
    }

}

inline int Sum(const vector<int> & argV){
	int val = 0;
	for (vector<int>::const_iterator it = argV.begin(); it != argV.end(); it++){
		val = val + (*it);
	}
	return val;

}


inline double Sum(const vector<double> & argV){
	double val = 0;
	for (vector<double>::const_iterator it = argV.begin(); it != argV.end(); it++){
		val = val + (*it);
	}
	return val;

}

//inline string GetTimeCompactStr()
//{
//	time_t timeObj;
//	time(&timeObj);
//	tm *pTime = gmtime(&timeObj);
//	char buffer[100];
//	sprintf(buffer, "%d%d%d%d", pTime->tm_yday , pTime->tm_hour, pTime->tm_min, pTime->tm_sec);
//	return buffer;
//}

template<typename T>
inline T InnerProd(const vector<T> & v1, const vector<T> & v2){
	assert(v1.size()==v2.size());
	T rs=0;
	for(int k=0; k<v1.size(); k++){
		rs=rs+v1[k]*v2[k];
	}
	return rs;
}


inline vector<double> Renormalize(const vector<double> & v, int idxStart, int idxEnd){
	
	vector<double> rs;
	assert(idxEnd<v.size());
	double sum=0;
	for(int k=idxStart; k<=idxEnd; k++){
		sum=sum+v[k];
	}
	for(int k=idxStart; k<=idxEnd; k++){
		rs.push_back(v[k]/sum);
	}
	return rs;
}

template<typename T>
inline vector<T> SubVector(const vector<T> & v,  int idxStart, int idxEnd){
	vector<T> newVec;
	for(int k=idxStart; k<=idxEnd; k++){
		newVec.push_back(v[k]);
	}
	return newVec;
}


inline vector<double> LinearCombTwoVec(const vector<double> & v1, const vector<double> & v2, double weightV1){
	assert(v1.size()==v2.size());
	vector<double> r;
	for(int k=0;k<v1.size(); k++){
		r.push_back(weightV1* v1[k]+(1-weightV1)*v2[k]);
	}
	return r;
}


inline vector<string> SplitVec(string str, char delimiter) {
	vector<string> internal;
	stringstream ss(str); // Turn the string into a stream.
	string tok;

	while (getline(ss, tok, delimiter)) {
		internal.push_back(tok);
	}

	return internal;
}


#endif	//_UTIL_H_///:~
