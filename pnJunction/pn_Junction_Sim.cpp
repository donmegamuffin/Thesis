// pn_Junction_Sim.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Node.h"
#include "Device_1D.h"
#include <random>

#include <shlwapi.h>



int main()
{
	std::size_t noOfIterations = 10;
	std::size_t noPerIteration = 100000;
	std::size_t perIterationDeci = noPerIteration / 10;

	std::size_t simSections = 0;

	//Ask user how many simulation sections the want
	std::cout << "Please enter the number of section for simulations:\n";
	std::cin >> simSections;

	std::ofstream testStream;
	testStream.open("test.csv");
	Device_1D testDevice(1e-6,simSections,simSections/2,1e12,1e12);
	testDevice.PrintDeviceData(true);
	testDevice.fPrintDeviceData(true,testStream);
	/*
	for (std::size_t i = 0; i<testDevice.nAry.size(); i++)
	{
		testDevice.nAry[i].n = rand() % 1000;
		testDevice.nAry[i].p = rand() % 1000;
	}*/

	for (std::size_t j = 0; j < noOfIterations; j++)
	{
		//Print iteration number
		std::cout << "Iteration " << j << "/" << noOfIterations << ":";
		for (std::size_t i = 0; i < noPerIteration; i++)
		{
			testDevice.calculateVoltages();
			testDevice.calculateJnL(true,true);
			testDevice.calculateJpL(true,true);
			//testDevice.fPrintDeviceData(false, testStream);
			//Print a . at each 10%
			if (i%perIterationDeci == 0)
				std::cout << ".";
		}
		testDevice.fPrintDeviceData(false, testStream);
		std::cout << std::endl;
	}
		
	


	//Closing statements
	testStream.close();
	ShellExecute(NULL, _T("Open"), _T("C:\\Users\\UoLzh\\Documents\\Year 4 Project\\Code\\Feb\\pn_Junction_Sim\\pn_Junction_Sim\\test.csv"), NULL, NULL, SW_SHOWNORMAL);
    return 0;
}

