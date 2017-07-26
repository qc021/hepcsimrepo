#ifndef __RANDGEN_H__
#define __RANDGEN_H__

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>
#include<ctime>
#include<stdlib.h>

using namespace std;

double lcgrand(int stream);
void lcgrandst (long zset, int stream); // Set the current zrng for stream "stream" to zset. 
long lcgrandgt (int stream);


#endif