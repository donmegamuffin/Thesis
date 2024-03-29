// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

//CONSTS
const double q = 1.6021766208e-19;		//Charge on electron
const double e_0 = 8.854187817e-10;		//Vacuum permittivity F cm^-1
const double T = 300;					//Temperature
const double mu = 1000;					//mu, GaN electron mobility (cm^2 V^-1 s^-1)
const double kB = 1.38064852e-19;		//Boltzmann constant cm^2 kg s-2 K-1
const double D = kB * T*mu / q;			//Diffussion coefficient
const int LOOP_CAP = 10000;
// TODO: reference additional headers your program requires here
