// FebCapSim.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <Shlwapi.h>
#include <vector>
#include <string>

//Experiment Constants
const double I_n = 6.24e18; //Number of electrons in 1 amp of current
const double B = 2.4e-17;	//Recombination constant (m^3/s)
const double Vol = 1e-11;	//Approx Volume (m^3)
const double Vramp = 4e9;	//Ramping voltage (V/s)
const double R = 5;			//Ohmic Resistance of device (Ohms)
const double C = 1.602e-19;	//1 Coulomb
const double QE = 0.01;		//Quantum Effeciency guesstimate

const std::string filename = "FebCap.csv";
const std::string fileLoc = "C:\\Users\\UoLzh\\Documents\\Year 4 Project\\Code\\Feb\\pn_CapSim\\FebCapSim\\FebCapSim\\";

//User input handling
bool InputToBool(std::string UserInput);

//Tent map Voltage function
double calcV(double time, double transition_time);
int findPeakBin(std::vector<double> RadVector);
int FWHMfindFirstBinBelow(std::vector<double> RadVector, int PeakIndex);
int FWHMfindFirstBinAbove(std::vector<double> RadVector, int PeakIndex);
double GetRadFWHM(std::vector<double> RadVector, std::vector<double> TimeVector);

int main(int argc, char **argv)
{	
	//Begin program here

	//Time related runtime constants
	const double t_start = 0;
	double t_step;		//Step time	(s)
	double t_trans_max;	//Maximum transition time
	double t_trans_step;	//Transition time step


	//Bool for full data print-out
	bool bFullPrintout = false;
	std::string UserInput_temp = "";

	//Set output stream + Open file
	std::ofstream file;
	file.open(filename,'w');
	//Print Heading
	file << "t,V,J,Rrad,RradCum,dndt,n\n";
	
	double t_trans = 0.5e-9;			//Transition time to being ramping down (s)
	std::vector<double> t_transVector;	//Will contain a vector of all use transition times
	std::vector<double> RradCumVector;		//A vector containing all the different Cumulative rads
	std::vector<double> FWHMVector;		//Vector to store all FWHMs

	//Console start user interface control
	std::cout << "What is the maximum transition time wished for (in seconds):\t (Format: 1e-9 for 1 nanosecond): \t";
	std::cin >> t_trans_max;
	std::cout << std::endl;
	std::cout << "What is the step in transition time step increase wished for (in seconds):\t (Format: 1e-10 for 0.1 nanoseconds): \t";
	std::cin >> t_trans_step;
	std::cout << std::endl;
	std::cout << "What time step do you wish to run the simulation at:\t (Format: 1e-12 for 1 picosecond):\t";
	std::cin >> t_step;
	std::cout << std::endl;
	std::cout << "Do you wish for full data printout?:\t(Y/N):\t";
	std::cin >> UserInput_temp;
	bFullPrintout = InputToBool(UserInput_temp);	//Converts the user's input to a bool
	std::cout << std::endl;

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

		std::vector<double> t_vector = { 0 };
		std::vector<double> RradVector = { 0 };	//Stores all the radiations at each time interval

		//Print initial conditions
		if (bFullPrintout) { file << t_now << "," << V << "," << Jn << "," << Rrad << "," << Rrad_cum << "," << dndt << "," << n << std::endl; }
		//Whilst the number of electrons isnt lower than 0
		do
		{
			//Increment Time
			t_now += t_step;
			t_vector.emplace_back(t_now);	//Add time point to end of time vector
			//Calculate V(t) + Print
			V = calcV(t_now, t_trans);
			//Calculate Current + Print
			J = 2*V / R;
			Jn = (J*t_step) / C;
			//Set new n
			n += Jn;
			//Calculate Rrad + Print
			Rrad = pow(n, 2)*(B / Vol)*t_step;
			RradVector.emplace_back(Rrad);	//Add Rrad point at end of Rad vector 
			//Add onto cumulative radiation
			Rrad_cum += Rrad;
			n -= Rrad;
			//Calculate dndt + Print
			dndt = Jn - Rrad;
			if (bFullPrintout) { file << t_now << "," << V << "," << Jn << "," << Rrad << "," << Rrad_cum << "," << dndt << "," << n << std::endl; }
		} while (!(n < 0));
		//Write cumulative photons to end of data spread
		if (bFullPrintout) { file << "Total photons likely detected (Quantum efficiency factored in):," << Rrad_cum*QE << std::endl; }
		//Print FWHM at end of data
		double FWHM = GetRadFWHM(RradVector, t_vector);
		if (bFullPrintout) { file << "Rad FWHM:," << FWHM << std::endl; }
		FWHMVector.emplace_back(FWHM);

		//Write time to array
		t_transVector.emplace_back(t_trans);
		//Add the new cumulative radiative (*QE) onto the end of the data array
		RradCumVector.emplace_back(Rrad_cum*QE);
	}
	//Print the arrays at end of file
	//Time
	file << "Transition time array: " << std::endl;
	for (auto a : t_transVector)
	{	//For each element, write it to file
		file << a << ",";
	}
	file << std::endl;
	//Cumulative radiation detected
	file << "RradCumQE array: " << std::endl;
	for (auto a : RradCumVector)
	{	//For each element, write it to file
		file << a << ",";
	}
	file<<std::endl;
	//FWHM Vector Print
	file << "FWHM array: " << std::endl;
	for (auto a : FWHMVector)
	{	//For each element, write it to file
		file << a << ",";
	}
	file << std::endl;

	//stop writing, and open with excel
	file.close();
	ShellExecute(NULL, _T("Open"), _T("C:\\Users\\UoLzh\\Documents\\Year 4 Project\\Code\\Feb\\pn_CapSim\\FebCapSim\\FebCapSim\\FebCap.csv"), NULL, NULL, SW_SHOWNORMAL);
	return 0;

	

}
//Takes input string and converts it into
//Associated logical bool.
//If no logical bool found, will throw an error
bool InputToBool(std::string UserInput)
{
	if (UserInput == "Y"   ||
		UserInput == "y"   ||
		UserInput == "Yes" ||
		UserInput == "yes" ||
		UserInput == "true"||
		UserInput == "True"||
		UserInput == "1"){
		return true;}

	else if (UserInput == "N" ||
		UserInput == "n"	  ||
		UserInput == "No"	  ||
		UserInput == "no"	  ||
		UserInput == "false"  ||
		UserInput == "False"  ||
		UserInput == "0") {
		return false;
	}
	else
	{
		throw std::logic_error("Please respond with a boolean response.\n");
		
	}
}

//Effectively Calculates the voltage tent map
double calcV(double time, double transition_time)
{
	if (time < transition_time)
		return Vramp*time;
	else
		return Vramp*((2 * transition_time) - time);
}

//Goes through all the curve datapoints to retrieve highest value one
int findPeakBin(std::vector<double> RadVector)
{
	double PeakValue = 0;
	int PeakValueAssociatedIndex = 0;
	for (std::size_t i = 0; i < RadVector.size(); i++)
	{
		if (RadVector[i] > PeakValue)
		{
			PeakValue = RadVector[i];
			PeakValueAssociatedIndex = i;
		}
	}
	return PeakValueAssociatedIndex;
}
//Goes through all datapoints from peak backwards to retrieve first bin below half Peak
int FWHMfindFirstBinBelow(std::vector<double> RadVector, int PeakIndex)
{
	int Index = 0;
	for (std::size_t i = PeakIndex; i >= 0; i--)
	{
		if (RadVector[i] < (RadVector[PeakIndex]/2))
		{
			Index = i;
			break;
		}
	}
	return Index;
}


//Goes through all datapoints from peak forwards to retrieve first bin below half Peak
int FWHMfindFirstBinAbove(std::vector<double> RadVector, int PeakIndex)
{
	int Index = 0;
	for (std::size_t i = PeakIndex; i < RadVector.size(); i++)
	{
		if (RadVector[i] < (RadVector[PeakIndex]/2))
		{
			Index = i;
			break;
		}
	}
	return Index;
}

/*Takes in Vector of radiation values, and associated time vector.
* finds FWHM of radiation curve, and then applies that to associated
* time vector to find how much time elapsed between these two points.
* Returns duration of the Full width half maximum in whatever unit time vector is in.
*/
double GetRadFWHM(std::vector<double> RadVector, std::vector<double> TimeVector)
{
	//Get the peak of the curve
	int PeakBin = findPeakBin(RadVector);
	//Get the first bin below half max
	int LowBin = FWHMfindFirstBinBelow(RadVector, PeakBin);
	//Get the first bin after half max
	int HighBin = FWHMfindFirstBinAbove(RadVector, PeakBin);
	//Return the difference in time between the associated time vector bins
	return (TimeVector[HighBin] - TimeVector[LowBin]);
}