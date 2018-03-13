// pn_Junction_Sim.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Node.h"
#include "Device_1D.h"
#include <random>
#include <thread>

#include <shlwapi.h>

void csv2dic(std::ifstream &csvFileStream);

int main()
{
	//Instantiate basic simulation params
	std::size_t noOfIterations = 0;
	std::size_t noPerIteration = 0;
	std::size_t simSections = 0;
	double exchangeScale = 1;

	//Ask user how many simulation sections the want
	std::cout << "\nPlease enter the number of main loops:\n";
	std::cin >> noOfIterations;
	std::cout << "\nPlease enter the number of iterations per loop (in thousands):\n";
	std::cin >> noPerIteration;
	std::cout << "\nPlease enter the Current Scaling (higher = faster, less Precise):\n";
	std::cin >> exchangeScale;
	noPerIteration = noPerIteration * 1000;

	std::size_t perIterationDeci = noPerIteration / 10;

	std::ofstream testStream;
	testStream.open("test.csv");

	/*CONVERT AND BUILD DEVICE FROM DATA*/
	//Open CSV
	std::ifstream csvInStream;
	csvInStream.open("Device.csv");
	//Convert to dic
	csv2dic(csvInStream);
	csvInStream.close();
	//Build device from dic data
	std::ifstream dicInStream;
	dicInStream.open("Device.dic");
	//Build Device
	Device_1D testDevice(dicInStream);

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
		std::cout << "Iteration " << j+1 << "/" << noOfIterations << ":";
		for (std::size_t i = 0; i < noPerIteration; i++)
		{
			testDevice.calculateVoltages();
			testDevice.calculateJnLEC(exchangeScale);
			testDevice.calculateJpLEV(exchangeScale);
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

void csv2dic(std::ifstream &csvFileStream)
{
	std::ofstream dicFileStream;

	dicFileStream.open("Device.dic");

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
}
