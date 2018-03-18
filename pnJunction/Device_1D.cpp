#include "stdafx.h"
#include "Device_1D.h"
#include "Node.h"
#include <thread>

//------------------------------CONSTRUCTORS--------------------------------------------
Device_1D::Device_1D()
{
	length = 1;	//prevents any weird div by 0 errors
	nodeWidth = 0;
	QE = 0;
	A = 0;
}

Device_1D::Device_1D(double Length, std::size_t noOfNodes, std::size_t transitionNode, double NdDensity, double NaDensity)
{
	length = Length;
	nodeWidth = length / noOfNodes;
	nAry.resize(noOfNodes);
	/*For loops for before and after the pn junction transition from p to n
	*Loop 1:
	*		Loop upto the transition node filling with p type
	*/
	for (std::size_t i = 0; i < transitionNode; i++)
	{
		nAry[i] = { 0, NaDensity, 0, NaDensity };
	}
	//Loop 2:
	//		Loop from transition block and fill with n type
	for (std::size_t i = transitionNode; i < noOfNodes; i++)
	{
		nAry[i] = { NdDensity, 0, NdDensity, 0 };
	}
}

Device_1D::Device_1D(std::ifstream &instream)
{
	//Takes file input for number of nodes to store
	std::size_t numberOfNodes = 0;
	//Read the length and number of nodes from top of file
	instream >> length >> numberOfNodes >> A >> QE >> R;
	//Resize the Node array to correct number
	nAry.resize(numberOfNodes);
	//Calculate the nodewidth
	nodeWidth = length / nAry.size();
	//Begin going line by line to input the different node variables
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].n;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].p;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].Nd;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].Na;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].V;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].Ec;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		instream >> nAry[i].Ev;
	}
}

Device_1D::~Device_1D()
{
}

//-----------------------------FILE SAVE & LOAD ----------------------

bool Device_1D::loadState(std::string fileName)
{
	std::ifstream csvInStream;
	csvInStream.open(fileName + ".csv");
	csv2dic(csvInStream, fileName);
	csvInStream.close();

	std::ifstream dicInStream;
	dicInStream.open(fileName + ".dic");
	//Takes file input for number of nodes to store
	std::size_t numberOfNodes = 0;
	//Read the length and number of nodes from top of file
	dicInStream >> length >> numberOfNodes >> QE >> A >> R;
	//Resize the Node array to correct number
	nAry.resize(numberOfNodes);
	//Calculate the nodewidth
	nodeWidth = length / nAry.size();
	//Begin going line by line to input the different node variables
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].n;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].p;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].Nd;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].Na;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].V;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].Ec;
	}
	for (std::size_t i = 0; i < numberOfNodes; i++)
	{
		dicInStream >> nAry[i].Ev;
	}
	dicInStream.close();
	//TODO
	//For now, find a way to return true or false depending on whether it loaded the file correctly
	return true;
}

void Device_1D::fSaveDevice(std::ofstream & saveStream)
{
	saveStream << length << "," << nAry.size() << "," << QE << "," << A << "," << R << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.n << ",";
	}
	saveStream << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.p << ",";
	}
	saveStream << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.Nd << ",";
	}
	saveStream << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.Na << ",";
	}
	saveStream << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.V << ",";
	}
	saveStream << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.Ec << ",";
	}
	saveStream << std::endl;
	for (auto a : nAry)
	{
		saveStream << a.Ev << ",";
	}
	saveStream << std::endl;
}

//-----------J current calculations ------------------------------------------------------

double Device_1D::calculateJnLEC(double exchangeScale)
{
	double JnL_cum = 0;	//Cumulative current from all nodes calculations
	//Loop through each of the nodes(except first)
	for (std::size_t i = 1; i < nAry.size(); i++)
	{
		double dn = (nAry[i].n - nAry[i - 1].n) / nodeWidth;
		double dV = (nAry[i].V - nAry[i - 1].V) / nodeWidth;
		double dEc = (nAry[i].Ec - nAry[i - 1].Ev) / nodeWidth;
		double n = (nAry[i].n + nAry[i - 1].n) / 2;
		double JnL = (-q * mu*n*dV) - (mu*n*dEc) + (q*D*dn);
		//double JnL = mu*(kB*T*((nAry[i].n - nAry[i - 1].n) / nodeWidth) - (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n) ) + ((nAry[i].n + nAry[i-1].n)/2)*((nAry[i].Ec-nAry[i-1].Ec)/nodeWidth);	//o
		nAry[i].n = nAry[i].n - (JnL*exchangeScale);
		nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
		JnL_cum += abs(JnL);
	}
	return JnL_cum;
}

double Device_1D::calculateJnREC(double exchangeScale)
{
	double JnL_cum = 0;	//Cumulative current from all nodes calculations
	//Loop through each of the nodes(except first)
	for (std::size_t i = nAry.size() - 1; i > 0; i--)
	{
		double dn = (nAry[i].n - nAry[i - 1].n) / nodeWidth;
		double dV = (nAry[i].V - nAry[i - 1].V) / nodeWidth;
		double dEc = (nAry[i].Ec - nAry[i - 1].Ev) / nodeWidth;
		double n = (nAry[i].n + nAry[i - 1].n) / 2;
		double JnL = (-q * mu*n*dV) - (mu*n*dEc) + (q*D*dn);
		//double JnL = mu*(kB*T*((nAry[i].n - nAry[i - 1].n) / nodeWidth) - (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n)) + ((nAry[i].n + nAry[i - 1].n) / 2)*((nAry[i].Ec - nAry[i - 1].Ec) / nodeWidth); //o
		nAry[i].n = nAry[i].n - (JnL*exchangeScale);
		nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
		JnL_cum += abs(JnL);
	}
	return JnL_cum;
}

double Device_1D::calculateJpLEV(double exchangeScale)
{
	double JpL_cum = 0;	//Cumulative current from all nodes calculations
	//Loop through each of the nodes(except first)
	for (std::size_t i = 1; i < nAry.size(); i++)
	{
		double dp = (nAry[i].p - nAry[i - 1].p) / nodeWidth;
		double dV = (nAry[i].V - nAry[i - 1].V) / nodeWidth;
		double dEv = (nAry[i].Ev - nAry[i - 1].Ev) / nodeWidth;
		double p = (nAry[i].p + nAry[i - 1].p) / 2;
		double JpL = (q * mu*p*dV) - (mu*p*dEv) + (q*D*dp);
		//double JpL = mu*(kB*T*((nAry[i].p - nAry[i - 1].p) / nodeWidth) + (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].p + nAry[i - 1].p)) + ((nAry[i].p + nAry[i - 1].p) / 2)*((nAry[i].Ev - nAry[i-1].Ev) / nodeWidth); //o
		nAry[i].p = nAry[i].p - (JpL*exchangeScale);
		nAry[i - 1].p = nAry[i - 1].p + (JpL*exchangeScale);
		JpL_cum += abs(JpL);
	}
	return JpL_cum;
}

double Device_1D::calculateJpREV(double exchangeScale)
{
	double JpL_cum = 0;	//Cumulative current from all nodes calculations
						//Loop through each of the nodes(except first)
	for (std::size_t i = nAry.size() - 1; i > 0; i--)
	{
		double dp = (nAry[i].p - nAry[i - 1].p) / nodeWidth;
		double dV = (nAry[i].V - nAry[i - 1].V) / nodeWidth;
		double dEv = (nAry[i].Ev - nAry[i - 1].Ev) / nodeWidth;
		double p = (nAry[i].p + nAry[i - 1].p) / 2;
		double JpL = (q * mu*p*dV) - (mu*p*dEv) + (q*D*dp);
		//double JpL = mu * (kB*T*((nAry[i].p - nAry[i - 1].p) / nodeWidth) + (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].p + nAry[i - 1].p)) + ((nAry[i].p + nAry[i - 1].p) / 2)*((nAry[i].Ev - nAry[i - 1].Ev) / nodeWidth); //o
		nAry[i].p = nAry[i].p - (JpL*exchangeScale);
		nAry[i - 1].p = nAry[i - 1].p + (JpL*exchangeScale);
		JpL_cum += abs(JpL);
	}
	return JpL_cum;
}

void Device_1D::calculateVoltages()
{
	//Loop through each of the nodes
	for (std::size_t i = 0; i < nAry.size(); i++)
	{
		//If node is at a cap position change calculation
		if (i == 0)
		{
			//TODO Do a different calculation
			nAry[i].V = pow(nodeWidth, 2)*((q / (2 * e_0))*(nAry[i].p + nAry[i].Nd - nAry[i].n - nAry[i].Na) + nAry[i + 1].V);
		}
		else if (i == (nAry.size() - 1))
		{
			nAry[i].V = pow(nodeWidth, 2)*((q / (2 * e_0))*(nAry[i].p + nAry[i].Nd - nAry[i].n - nAry[i].Na) + nAry[i - 1].V);
		}
		else	//Do normal calculation
		{
			double a2 = pow(nodeWidth, 2);
			double C = (q*a2) / e_0;
			double netCharge = nAry[i].p + nAry[i].Nd - nAry[i].n - nAry[i].Na;
			nAry[i].V = 0.5*(C*(netCharge)+nAry[i - 1].V + nAry[i + 1].V);
			//nAry[i].V = pow(nodeWidth, 2)*((q / (2 * e_0))*(nAry[i].p + nAry[i].Nd - nAry[i].n - nAry[i].Na) + nAry[i + 1].V + nAry[i - 1].V); //O
		}
	}
}

void Device_1D::cancelCharges()
{
	//For each node in the array, take ALL of the minor charge carrier away from both major and minor
	//Loop through each of the nodes
	for (std::size_t i = 0; i < nAry.size(); i++)
	{
		if (nAry[i].n < nAry[i].p)
		{
			nAry[i].p = nAry[i].p - nAry[i].n;
			nAry[i].n = 0;
		}
		else
		{
			nAry[i].n = nAry[i].n - nAry[i].p;
			nAry[i].p = 0;
		}
	}
}

//----------------------SIMULATION STUFF-------------------------------------------------

void Device_1D::injectCharges(double CurrentDensity, double injectionDuration)
{
	//Assumes p-n style device. Inject n in right, p on left.
	nAry[0].n += (CurrentDensity*injectionDuration);
	nAry[nAry.size() - 1].p += (CurrentDensity*injectionDuration);
	return;
}

double Device_1D::calcRadRecombine(double timeScale)
{
	double RadCum = 0;
	for (std::size_t i = 0; i < nAry.size(); i++)
	{
		RadCum += B * nAry[i].n*nAry[i].p;
		nAry[i].n -= sqrt(B)*nAry[i].n;
		nAry[i].p -= sqrt(B)*nAry[i].p;
	}
	return RadCum;
}

double Device_1D::inputV(double time, double transition_time)
{
	double Vramp = 4e9;
	if (time < transition_time)
		return Vramp * time;
	else
		return Vramp * ((2 * transition_time) - time);
}

void Device_1D::simulateDevice_pn(double & outFWHM, double & outRrad, double t_step, double t_trans)
{
	double t_now = 0;	//current device simulation time
	std::vector<double> t_Vec;	//Stores all the times
	double Rrad_now;
	std::vector<double> Rrad_Vec;	//Current device radiation
	double Rrad_cum = 0;	//Cumulative radiation
	double inV = 0;			//Input voltage
	double J = 0;			//J

	//Main function loop
	do
	{
		t_now += t_step;
		t_Vec.emplace_back(t_now);

		inV = inputV(t_now, t_trans);
		J = inV / (A * R);

		injectCharges(J, t_step);
		calculateVoltages();
		calculateJpREV(t_step);
		calculateJnLEC(t_step);

		Rrad_now = calcRadRecombine(t_step);
		Rrad_Vec.emplace_back(Rrad_now);
		Rrad_cum += Rrad_now;
	} while (Rrad_now>0);

	outFWHM = GetRadFWHM(Rrad_Vec, t_Vec);
	outRrad = Rrad_cum;
}

void Device_1D::simulateDevice(double t_trans, double t_step, std::string radOutFileName)
{
	//Open file with designated name
	std::ofstream radOutFile;

	//DEBUG REMOVE ME LATER
	std::ofstream debugFile;
	debugFile.open("DEBUG.csv");

	radOutFile.open(radOutFileName + ".csv");
	double t_now = 0;
	double J = 0;
	double inV = 0;
	radOutFile << "t,V,Rrad\n";
	for (int i = 0; i < 500; i++)
	{
		t_now = i * t_step;
		inV = inputV(t_now, t_trans);
		J = inV / (A * R);
		injectCharges(J, t_step);
		calculateVoltages();
		//Assuming P - N device
		//In p in from left, n from right
		calculateJpREV(t_step);
		calculateJnLEC(t_step);
		radOutFile << t_now << "," << inV << "," << calcRadRecombine(t_step) << std::endl;

		//Print device state to debug
		for (auto a : nAry)
		{
			debugFile << a.n << ",";
		}
		for (auto a : nAry)
		{
			debugFile << a.p << ",";
		}
		debugFile << std::endl;
	}
	//Make sure to close file after doing calculations
	radOutFile.close();

	debugFile.close();
}

void Device_1D::bringToEqm(double Tolerance, double exchangeScale, bool bPrintCurent)
{
	//Initialise total current of device
	double Jtot = 1e18;
	//Whilst the current is higher than the current equilibrium limit
	//Loop until Jtot is below user defined limit
	while (Jtot > Tolerance)
	{
		//Find total current of charges being shifted across entire device
		//by running the two current calculations across entire device
		Jtot = calculateJnLEC(exchangeScale) + calculateJpLEV(exchangeScale);
		//Rebalance voltages
		calculateVoltages();
		//Make sure to cancel out any charges that would annihilate each other in device
		cancelCharges();
		//OPTIONAL : Print the current found to console
		if (bPrintCurent) 
		{ 
			//Bit of a bodge, clears line and updates it with new Jtot
			std::cout << "\r" << Jtot << "                 ";
		}
	}
	if (bPrintCurent) { std::cout << "\n Equilibrium calculation finished." << std::endl; }
	return;
}

void Device_1D::bringToEqm(double Tolerance, double exchangeScale, std::ostream &fileOutput, bool bPrintConsole)
{
	//Initialise total current of device
	double Jtot = 1e13;
	double J_prev = 1e8;
	//Whilst the current is higher than the current equilibrium limit
	//Loop until Jtot is below user defined limit
	while (Jtot > Tolerance && abs(Jtot - J_prev)>(0.1*Tolerance))
	{
		J_prev = Jtot;	//So the program won't enter endless loop unable to get below tolerance if grinds to convergance
		//Find total current of charges being shifted across entire device
		//by running the two current calculations across entire device
		Jtot = calculateJnREC(exchangeScale) + calculateJpLEV(exchangeScale);
		//Rebalance voltages
		calculateVoltages();
		//Make sure to cancel out any charges that would annihilate each other in device
		cancelCharges();
		//OPTIONAL : Print the current found to console
		if (bPrintConsole) 
		{	//Bit of a bodge, clears line and updates it with new Jtot
			std::cout << "\r" << Jtot << "                 ";
		}
		//Print current to file
		fileOutput << Jtot << std::endl;
	}
	if (bPrintConsole) { std::cout << "\n Equilibrium calculation finished." << std::endl; }
	return;
}

void Device_1D::fullSim(std::string eqmFileName, double timeStep, double transStep, double transMax)
{
	//Simulation Data arrays for output
	std::vector<double> t_trans_Vec;
	std::vector<double> FWHM_Vec;
	std::vector<double> Rcum_Vec;

	double FWHM_now = 0;	//Holds the loop FWHM
	double Rcum_now = 0;	//Holds the loop Rcum
	
	//Open filestream for printing results
	std::ofstream resultsFile;
	resultsFile.open(eqmFileName + "_Rad_results.csv",'w');

	for (int i = 0; i*transStep < transMax; i++)
	{
		//Set up device
		loadState(eqmFileName);

		//Add the transition time to trans vector
		t_trans_Vec.emplace_back(i*transStep);

		simulateDevice_pn(FWHM_now, Rcum_now, timeStep, i*transStep);

		FWHM_Vec.emplace_back(FWHM_now);
		Rcum_Vec.emplace_back(Rcum_now);
	}
	//Print Results to file
	for (auto a : t_trans_Vec)
	{
		resultsFile << a << ",";
	}
	resultsFile << std::endl;
	for (auto a : FWHM_Vec)
	{
		resultsFile << a << ",";
	}
	resultsFile << std::endl;
	for (auto a : Rcum_Vec)
	{
		resultsFile << a << ",";
	}
	resultsFile << std::endl;

	resultsFile.close();
}
//----------------------MISC TOOLS STUFF-------------------------------------------------

void Device_1D::PrintDeviceData(bool bPrintHeader)
{
	/*Prints to console headings in format: n1, n2..., until V[end]
	* in order n, p, V. All n, p, and V are together
	* for length of entire of entire device.
	*/
	if (bPrintHeader)
	{	//TODO must be a more elegant way to do this?
		//Print n labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			std::cout << "n" << a << ",";
		}
		//Print p labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			std::cout << "p" << a << ",";
		}
		//Print V labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			std::cout << "V" << a << ",";
		}
		std::cout << std::endl;
	}
	//Print the data to console
	for (auto a : nAry)
	{
		std::cout << a.n << ",";
	}
	for (auto a : nAry)
	{
		std::cout << a.p << ",";
	}
	for (auto a : nAry)
	{
		std::cout << a.V << ",";
	}
	std::cout << std::endl;
}

void Device_1D::fPrintDeviceData(bool bPrintHeader, std::ofstream& outStream)
{
	if (bPrintHeader)
	{	//TODO must be a more elegant way to do this?
		//Print n labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			outStream << "n" << a << ",";
		}
		//Print p labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			outStream << "p" << a << ",";
		}
		//Print V labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			outStream << "V" << a << ",";
		}
		outStream << std::endl;
	}

	//Print data to csv in n, p, V order (type-grouped still)
	for (auto a : nAry)
	{
		outStream << a.n << ",";
	}
	for (auto a : nAry)
	{
		outStream << a.p << ",";
	}
	for (auto a : nAry)
	{
		outStream << a.V << ",";
	}
	outStream << std::endl;
}

void Device_1D::csv2dic(std::ifstream &csvFileStream, std::string fileName)
{
	std::ofstream dicFileStream;

	dicFileStream.open(fileName + ".dic");

	char currentChar;

	while (csvFileStream.get(currentChar))
	{
		if (currentChar == ',')
		{
			dicFileStream << ' ';
		}
		else if (csvFileStream.eof())
		{
			dicFileStream.close();
		}
		else
		{
			dicFileStream << currentChar;
		}
	}
	dicFileStream.close();
}

//Goes through all the curve datapoints to retrieve highest value one
int Device_1D::findPeakBin(std::vector<double> RadVector)
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
int Device_1D::FWHMfindFirstBinBelow(std::vector<double> RadVector, int PeakIndex)
{
	int Index = 0;
	for (std::size_t i = PeakIndex; i >= 0; i--)
	{
		if (RadVector[i] < (RadVector[PeakIndex] / 2))
		{
			Index = i;
			break;
		}
	}
	return Index;
}


//Goes through all datapoints from peak forwards to retrieve first bin below half Peak
int Device_1D::FWHMfindFirstBinAbove(std::vector<double> RadVector, int PeakIndex)
{
	int Index = 0;
	for (std::size_t i = PeakIndex; i < RadVector.size(); i++)
	{
		if (RadVector[i] < (RadVector[PeakIndex] / 2))
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
double Device_1D::GetRadFWHM(std::vector<double> RadVector, std::vector<double> TimeVector)
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
