// FebCapSim.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <Shlwapi.h>

//Experiment Constants
const double I_n = 6.24e18; //Number of electrons in 1 amp of current
const double B = 2.4e-17;	//Recombination constant (m^3/s)
const double Vol = 1e-13;	//Approx Volume (m^3)
const double Vramp = 4e9;	//Ramping voltage (V/s)
const double R = 5;			//Ohmic Resistance of device (Ohms)
const double C = 1.602e-19;	//1 Coulomb
const double QE = 0.01;		//Quantum Effeciency guesstimate

//Time related 
const double t_start = 0;
const double t_step = 1e-10;		//Step time	(s)
const double t_trans_max = 10e-9;	//Maximum transition time
const double t_trans_step = 0.5e-9;	//Transition time step


const std::string filename = "FebCap.csv";
const std::string fileLoc = "C:\\Users\\UoLzh\\Documents\\Year 4 Project\\Code\\Feb\\pn_CapSim\\FebCapSim\\FebCapSim\\";

double calcV(double time, double transition_time);



int main()
{	//Begin program here
	//Set output stream + Open file
	std::ofstream file;
	file.open(filename,'w');
	//Print Heading
	file << "t,V,J,Rrad,RradCum,dndt,n\n";
	
	double t_trans = 0.5e-9;	//Transition time to being ramping down (s)

	//Main body loop
	//Cycle through transition times
	for (std::size_t i = 0; (i*t_trans_step) < t_trans_max; i++)
	{
		//Update transition time
		t_trans = i*t_trans_step;

		//Set initial conditions
		double V = 0;		//Voltage applied (V)
		double J = 0;		//Current from external voltage (A)
		double Jn = 0;		//Curent in terms of 
		double n = 0;		//number of additional electrons
		double Rrad = 0;	//Radiative recombination
		double Rrad_cum = 0;//Cumulative photons
		double dndt = 0;	//Initialise change in electrons
		double t_now = 0;	//Time of experiment at that point
						
		//Print initial conditions
		file << t_now << "," << V << "," << Jn << "," << Rrad << "," << Rrad_cum << "," << dndt << "," << n << std::endl;
		//Whilst the number of electrons isnt lower than 0
		do
		{
			//Increment Time
			t_now += t_step;
			//Calculate V(t) + Print
			V = calcV(t_now, t_trans);
			//Calculate Current + Print
			J = V / R;
			Jn = (J*t_step) / C;
			//Set new n
			n += Jn;
			//Calculate Rrad + Print
			Rrad = pow(n, 2)*(B / Vol)*t_step;
			//Add onto cumulative radiation
			Rrad_cum += Rrad;
			n -= sqrt(Rrad);
			//Calculate dndt + Print
			dndt = Jn - sqrt(Rrad);
			file << t_now << "," << V << "," << Jn << "," << Rrad << "," << Rrad_cum << "," << dndt << "," << n << std::endl;
		} while (!(n < 0));
		file << "Total photons likely detected (Quantum efficiency factored in):," << Rrad_cum*QE;
		file << std::endl;
	}
	
	//stop writing, and open with excel
	file.close();
	ShellExecute(NULL, _T("Open"), _T("C:\\Users\\UoLzh\\Documents\\Year 4 Project\\Code\\Feb\\pn_CapSim\\FebCapSim\\FebCapSim\\FebCap.csv"), NULL, NULL, SW_SHOWNORMAL);
	return 0;
}
//Effectively Calculates the voltage tent map
double calcV(double time, double transition_time)
{
	if (time < transition_time)
		return Vramp*time;
	else
		return Vramp*((2 * transition_time) - time);
}
