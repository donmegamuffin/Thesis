/* This class organises the dimensions 
* and interactions between different nodes 
* (i.e. currents of carriers between each other)
*/

#pragma once
#include <vector>

class Node;
class Device_1D
{
public:
	/*Constructors & Destructors*/
	//DEPRECIATED: Old constructor when building simple pn junction from small numbers. Depreciated but helpful
	Device_1D( double Length, std::size_t noOfNodes, std::size_t transitionNode, double NdDensity, double NaDensity);
	//Construct device from file
	Device_1D(std::ifstream &instream);
	//Blank Device Constructor.
	Device_1D();
	//Destructor
	~Device_1D();
	
private:
	double length;				//Length of device in m
	double nodeWidth;			//Calculated length of each node
	std::vector<Node> nAry;		//Data storage for all nodes
	const double B = 2.4e-17;	//Rad recombine constant
	double QE;
	double A;

public:
//Printing functions
	//Print data to console in n, p, V order (type-grouped)
	//Can print a heading if needed
	void PrintDeviceData(bool bPrintHeader);
	//Print data to csv in n, p, V order (type-grouped)
	//Can print a heading if needed
	void fPrintDeviceData(bool bPrintHeader, std::ofstream& outStream);
	
//Voltage calculation functions
	void calculateVoltages();
//Charge transfer functions
	//Full Jn/Jp calculations
	/*
	Does full current and charge density exchange for ELECTRONS
	for all nodes in the device from Left to Right
	@param exchangeScale - aproximation scaling modifier. Higher = faster but less precise.
	*/
	double calculateJnLEC(double exchangeScale);
	/*
	Does full current and charge density exchange for HOLES
	for all nodes in the device from Left to Right
	@param exchangeScale - aproximation scaling modifier. Higher = faster but less precise.
	*/
	double calculateJpLEV(double exchangeScale);

//Equilibrium Methods
	/* Runs current calculations on device until combined current is less than the tolerance
		@param Tolerance:	total current considered to be EQM
		@param exchangeScale: How approximate want the calculation. Large = More approx but faster.
				A good value for a simple device about ~1000000
		@param bPrintCurrent: if true, will print the total current to console.
	*/
	void bringToEqm(double Tolerance, double exchangeScale, bool bPrintCurent);
	
	/* Runs current calculations on device until combined current is less than the tolerance and prints to file
		@param Tolerance:	total current considered to be EQM
		@param exchangeScale: How approximate want the calculation. Large = More approx but faster.
			A good value for a simple device about ~1000000
		@param fileOutput:	The output stream used to print results to.
		@param bPrintCurrent: if true, will print the total current to console.
		
	*/
	void bringToEqm(double Tolerance, double exchangeScale, std::ostream &fileOutput, bool bPrintConsole);

	/* This function goes along each node and cancels out n from p and vice versa 
		leaving only the majority.
	*/
	void cancelCharges();

	void injectCharges(double CurrentDensity, double injectionDuration);


	/* Loads all device stats from a given filename
	___________________________________________________________________
		@param fileName: Input the filename (without .csv) to be loaded
		@return bool : Returns true if successful
	*/
	bool loadState(std::string fileName);

	/* SEE calculateJnLEC FOR DETAILS
		Same as normal function, but runs from R to L
	*/
	void calculateJnREC(double exchangeScale);

	/* Calculates equation R = Bnp on each node
		@param timescale: Over what duration this is happening
		@return double: Returns all rad recombination across all nodes on device
	*/
	double calcRadRecombine(double timeScale);

	
	/* Saves all details of device to same file standard as can be loaded (.csv)
		@param &saveStream: reference to file output stream to be sending the data.
	*/
	void fSaveDevice(std::ofstream &saveStream);

private:
	//Converts .CSV file to an internal .dic file
	void csv2dic(std::ifstream &csvFileStream, std::string fileName);

	//Effectively Calculates the voltage tent map [Copied from my cap sim]
	double calcV(double time, double transition_time);
};

