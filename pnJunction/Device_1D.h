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
	~Device_1D();
	double length;
	double nodeWidth;
	std::vector<Node> nAry;



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
	/*Calculate the left direction current for electrons
	DON'T CALL WITH FALSE FALSE, SILLY.
	@param bDiff : Enable Diffusion in calculation
	@param bDrift : Enable Drift in calculation
	*/
	void calculateJnL(bool bDiff, bool bDrift);
	/*Calculate the left direction current for holes
	DON'T CALL WITH FALSE FALSE, SILLY.
	@param bDiff : Enable Diffusion in calculation
	@param bDrift : Enable Drift in calculation
	*/
	void calculateJpL(bool bDiff, bool bDrift);

private:
	//Don't call me
	Device_1D();	
};
