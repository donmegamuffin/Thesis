#include "stdafx.h"
#include "Device_1D.h"
#include "Node.h"
#include <thread>




Device_1D::Device_1D()
{
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
	int numberOfNodes = 0;
	//Read the length and number of nodes from top of file
	instream >> length >> numberOfNodes;
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

//Print data to console in n, p, V order (type-grouped)
//Can print a heading if needed
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
//Print data to csv in n, p, V order (type-grouped)
//Can print a heading if needed
void Device_1D::fPrintDeviceData(bool bPrintHeader, std::ofstream& outStream)
{
	if (bPrintHeader)
	{	//TODO must be a more elegant way to do this?
		//Print n labels
		for (std::size_t a = 0; a < nAry.size(); a++)
		{
			outStream<< "n" << a << ",";
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


//Implementation of V(i) = (a^2q/2er*e0)(pi-ni+Ndi-Nai)+a^2Vi+1 + a^2Vi-1
void Device_1D::calculateVoltages()
{
	//Loop through each of the nodes
	for (std::size_t i = 0; i < nAry.size(); i++)
	{
		//If node is at a cap position change calculation
		if (i == 0 || i == (nAry.size() - 1))
		{
			//TODO Do a different calculation
			nAry[i].V = 0;
		}
		else	//Do normal calculation
		{
			nAry[i].V = pow(nodeWidth, 2)*((q / (2 * e_0))*(nAry[i].p + nAry[i].Nd - nAry[i].n - nAry[i].Na) + nAry[i + 1].V + nAry[i - 1].V);
		}
	}
}

void Device_1D::calculateJnL(bool bDiff, bool bDrift, int exchangeScale)
{
	//If both drift and diffusion are requested
	if (bDiff&&bDrift)
	{//Loop through each of the nodes
		for (std::size_t i = 0; i < nAry.size(); i++)
		{
			//If node is at a cap position change calculation
			if (i == 0 || i == (nAry.size() - 1))
			{
				//TODO Do a different calculation
				nAry[i].n = 0;
			}
			else	//Do normal calculation
			{
				double JnL = kB*T*mu*(((nAry[i].n - nAry[i - 1].n) / nodeWidth) - (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n));
				nAry[i].n = nAry[i].n - (JnL*exchangeScale);
				nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
			}
		}
	}
	//Only if you want Diffusion
	else if (bDiff)
	{
		for (std::size_t i = 0; i < nAry.size(); i++)
		{
			//If node is at a cap position change calculation
			if (i == 0 || i == (nAry.size() - 1))
			{
				//TODO Do a different calculation
				nAry[i].n = 0;
			}
			else	//Do normal calculation
			{
				double JnL = kB*T*mu*((nAry[i].n - nAry[i - 1].n) / nodeWidth);
				nAry[i].n = nAry[i].n - (JnL*exchangeScale);
				nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
			}
		}
	}
	//Only use Drift
	else if (bDrift)
	{
		for (std::size_t i = 0; i < nAry.size(); i++)
		{
			//If node is at a cap position change calculation
			if (i == 0 || i == (nAry.size() - 1))
			{
				//TODO Do a different calculation
				nAry[i].n = 0;
			}
			else	//Do normal calculation
			{
				double JnL = kB*T*mu*-(q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n);
				nAry[i].n = nAry[i].n - (JnL*exchangeScale);
				nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
			}
		}
	}
	//Incase some idiot decided to do nothing with the function by accident
	else
	{
		std::cout << "calculateJnL was called with False False. Go fix.\n";
	}
}

void Device_1D::calculateJpL(bool bDiff, bool bDrift, int exchangeScale)
{
	//If both drift and diffusion are requested
	if (bDiff&&bDrift)
	{
		//Loop through each of the nodes
		for (std::size_t i = 0; i < nAry.size(); i++)
		{
			//If node is at a cap position change calculation
			if (i == 0 || i == (nAry.size() - 1))
			{
				//TODO Do a different calculation
				nAry[i].p = 0;
			}
			else	//Do normal calculation
			{
				double JpL = kB*T*mu*(((nAry[i].p - nAry[i - 1].p) / nodeWidth) + (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].p + nAry[i - 1].p));
				nAry[i].p = nAry[i].p - (JpL*exchangeScale);
				nAry[i - 1].p = nAry[i - 1].p + (JpL*exchangeScale);
			}
		}
	}
	//Only if you want Diffusion
	else if (bDiff)
	{
		//Loop through each of the nodes
		for (std::size_t i = 0; i < nAry.size(); i++)
		{
			//If node is at a cap position change calculation
			if (i == 0 || i == (nAry.size() - 1))
			{
				//TODO Do a different calculation
				nAry[i].p = 0;
			}
			else	//Do normal calculation
			{
				double JpL = kB*T*mu*(((nAry[i].p - nAry[i - 1].p) / nodeWidth));
				nAry[i].p = nAry[i].p - (JpL*exchangeScale);
				nAry[i - 1].p = nAry[i - 1].p + (JpL*exchangeScale);
			}
		}
	}
	//Only use Drift
	else if (bDrift)
	{
		//Loop through each of the nodes
		for (std::size_t i = 0; i < nAry.size(); i++)
		{
			//If node is at a cap position change calculation
			if (i == 0 || i == (nAry.size() - 1))
			{
				//TODO Do a different calculation
				nAry[i].p = 0;
			}
			else	//Do normal calculation
			{
				double JpL = kB*T*mu*(q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].p + nAry[i - 1].p);
				nAry[i].p = nAry[i].p - (JpL*exchangeScale);
				nAry[i - 1].p = nAry[i - 1].p + (JpL*exchangeScale);
			}
		}
	}
	//Incase some idiot decided to do nothing with the function by accident
	//Dont use this form
	else
	{
		std::cout << "calculateJpL was called with False False. Go fix.\n";
	}
}

/*------------------------TEST---------------------------------------*/
double Device_1D::calculateJnLEC(double exchangeScale)
{
	double JnL_cum = 0;	//Cumulative current from all nodes calculations
	//Loop through each of the nodes(except first)
	for (std::size_t i = 1; i < nAry.size(); i++)
	{
		double JnL = mu*(kB*T*((nAry[i].n - nAry[i - 1].n) / nodeWidth) - (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n) ) + ((nAry[i].n + nAry[i-1].n)/2)*((nAry[i].Ec-nAry[i-1].Ec)/nodeWidth);
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
		double JpL = mu*(kB*T*((nAry[i].p - nAry[i - 1].p) / nodeWidth) + (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].p + nAry[i - 1].p)) + ((nAry[i].p + nAry[i - 1].p) / 2)*((nAry[i].Ev - nAry[i-1].Ev) / nodeWidth);
		nAry[i].p = nAry[i].p - (JpL*exchangeScale);
		nAry[i - 1].p = nAry[i - 1].p + (JpL*exchangeScale);
		JpL_cum += abs(JpL);
	}
	return JpL_cum;
}

//EXPERIMENTAL
void Device_1D::calculateJnLEC_MT(double exchangeScale)
{	//MT = Multithreaded
	
	//Loop through each of the nodes(except last)
	for (std::size_t i = 1; i < nAry.size(); i++)
	{
		double JnL = mu*(kB*T*((nAry[i].n - nAry[i - 1].n) / nodeWidth) - (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n)) + ((nAry[i].n + nAry[i - 1].n) / 2)*((nAry[i].Ec - nAry[i - 1].Ec) / nodeWidth);
		nAry[i].n = nAry[i].n - (JnL*exchangeScale);
		nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
	}
}

void Device_1D::injectCharges(double CurrentDensity, double injectionDuration)
{
	//Assumes p-n style device. Inject p in left, n on right.
	nAry[0].n += (CurrentDensity*injectionDuration);
	nAry[nAry.size() - 1].p += (CurrentDensity*injectionDuration);
	return;
}

void Device_1D::cancelCharges()
{
	//For each node in the array, take ALL of the minor charge carrier away from both major and minor
	//Loop through each of the nodes
	for (std::size_t i = 0; i < nAry.size(); i++)
	{
		if(nAry[i].n < nAry[i].p)	
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

void Device_1D::calculateJnREC(double exchangeScale)
{
	//Loop through each of the nodes(except first)
	for (std::size_t i = nAry.size()-1; i > 0; i--)
	{
		double JnL = mu*(kB*T*((nAry[i].n - nAry[i - 1].n) / nodeWidth) - (q / (2 * nodeWidth))*(nAry[i].V - nAry[i - 1].V)*(nAry[i].n + nAry[i - 1].n)) + ((nAry[i].n + nAry[i - 1].n) / 2)*((nAry[i].Ec - nAry[i - 1].Ec) / nodeWidth);
		nAry[i].n = nAry[i].n - (JnL*exchangeScale);
		nAry[i - 1].n = nAry[i - 1].n + (JnL*exchangeScale);
	}
}

double Device_1D::calcRadRecombine(double timeScale)
{
	double RadCum = 0;
	for (std::size_t i = 0; i < nAry.size(); i++)
	{
		RadCum += B*nAry[i].n*nAry[i].p; 
		nAry[i].n -= sqrt(B)*nAry[i].n;
		nAry[i].p -= sqrt(B)*nAry[i].p;
	}
	return RadCum;
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
		if (bPrintCurent) { std::cout << Jtot << std::endl; }
	}
	return;
}