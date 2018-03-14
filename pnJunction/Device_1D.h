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

	Device_1D( double Length, std::size_t noOfNodes, std::size_t transitionNode, double NdDensity, double NaDensity);
	//Construct device from file
	Device_1D(std::ifstream &instream);


	~Device_1D();
	double length;
	double nodeWidth;
	std::vector<Node> nAry;
	const double B = 2.4e-17;
	
	//EXPERIMENTAL
	//std::vector<std::thread> threadArray;


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
	/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	DEPRECIATED AND BUGGED!
	USE "calculateJnLEC(...)" INSTEAD!
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	Calculate the left direction current for electrons
	DON'T CALL WITH FALSE FALSE, SILLY.
	@param bDiff : Enable Diffusion in calculation
	@param bDrift : Enable Drift in calculation
	@param exchangeScale: Applies calculation this many times without full state recalc.
	*/
	void calculateJnL(bool bDiff, bool bDrift,int exchangeScale);
	
	/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	DEPRECIATED AND BUGGED!
	USE "calculateJpLEV(...)" INSTEAD!
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	Calculate the left direction current for holes
	DON'T CALL WITH FALSE FALSE, SILLY.
	@param bDiff : Enable Diffusion in calculation
	@param bDrift : Enable Drift in calculation
	@param exchangeScale: Applies calculation this many times without full state recalc.
	*/
	void calculateJpL(bool bDiff, bool bDrift,int exchangeScale);
	
	//Full Jn/Jp calculations
	/*
	Does full current and charge density exchange for ELECTRONS
	for all nodes in the device from Left to Right
	@param exchangeScale - aproximation scaling modifier. Higher = faster but less precise.
	*/
	void calculateJnLEC(double exchangeScale);
	/*
	Does full current and charge density exchange for HOLES
	for all nodes in the device from Left to Right
	@param exchangeScale - aproximation scaling modifier. Higher = faster but less precise.
	*/
	void calculateJpLEV(double exchangeScale);

	//EXPERIMENTAL
	void calculateJnLEC_MT(double exchangeScale);

	void injectCharges(double CurrentDensity, double injectionDuration);

	void cancelCharges();

	//Same as normal function, but runs from R to L
	void calculateJnREC(double exchangeScale);

	double calcRadRecombine(double timeScale);

private:
	//Don't call me
	Device_1D();	
};

