// pn_Junction_Sim.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Node.h"
#include "Device_1D.h"
#include <random>
#include <memory>
//#include <thread>
//#include <shlwapi.h>
#include <string>

//Bool for program to keep looping
bool bShouldContinue = true;
//Enums of all the different commands that can be called
enum class eInCmd { LOAD, SAVE, EQM, SIM, EXIT, HELP, INVALID };

//Instantiate pointer to Device
Device_1D Device;

eInCmd readCommand(std::string command);

double calcV(double time, double transition_time);

void callSave(std::string &csvSaveFileName, std::ofstream &csvSaveStream);

void callLoad(std::string &csvLoadFileName);

void callSim();

void loadFile(std::string FileName);

void callEqm();

//Print methods
void printWelcome();
void printHelp();
void printExit();

inline void flushCin();

int main()
{
	//Create in/output file streams for use later on
	std::ifstream csvLoadStream;		//Used when loading csv files
	std::string csvLoadFileName = "";	//Name of the load file
	std::ofstream csvSaveStream;		//Used when saving csv files
	std::string csvSaveFileName = "";	//Name of the save file



	//Print Welcome screen
	printWelcome();
	//Loop program whilst bShouldContinue is true
	while (bShouldContinue)
	{
		//Initialise input strings to blanks
		std::string userInput = "";
		csvLoadFileName = "";
		csvSaveFileName = "";

		std::cout << "CMD :>\t";
		std::cin >> userInput;
		eInCmd command = readCommand(userInput);
		
		//Executes code depending on command given
		switch (command)
		{
		case eInCmd::LOAD:
			callLoad(csvLoadFileName);
			flushCin();
			break;
		case eInCmd::SAVE:
			callSave(csvSaveFileName, csvSaveStream);
			flushCin();
			break;
		case eInCmd::EQM:
			callEqm();
			flushCin();
			break;
		case eInCmd::SIM:
			callSim();
			flushCin();
			break;
		case eInCmd::EXIT:
			printExit();
			bShouldContinue = false;
			flushCin();
			break;
		case eInCmd::HELP:
			printHelp();
			flushCin();
			break;
		case eInCmd::INVALID:
			std::cout << "Command not recognised! Please type \"help\" for command list!\n";
			flushCin();
			break;
		default:
			std::cout << "Something went terribly wrong for you to get this message...\n";
			flushCin();
			break;
		}
	}
    return 0;
}

void callLoad(std::string &csvLoadFileName)
{
	std::cin >> csvLoadFileName;
	if (csvLoadFileName != "")
	{
		loadFile(csvLoadFileName);
		std::cout << "Lad successful! Filename: \"" << csvLoadFileName << ".csv\". Please check to ensure it's correct!\n";
	}
	else
	{
		std::cout << "Filename not recognised!";
	}
}

void callSim()
{
	double t_trans = 0;
	double t_step = 0;
	std::string radOutFileName;

	Device.simulate(t_trans,t_step,radOutFileName);
}

void callSave(std::string &csvSaveFileName, std::ofstream &csvSaveStream)
{
	std::cin >> csvSaveFileName;
	if (csvSaveFileName != "")
	{
		csvSaveStream.open(csvSaveFileName + ".csv");
		Device.fSaveDevice(csvSaveStream);
		csvSaveStream.close();
		std::cout << "Save successful! Filename: \"" << csvSaveFileName << ".csv\". Please check to ensure it's correct!\n";
	}
	else
	{
		std::cout << "Filename not recognised!";
	}
}



eInCmd readCommand(std::string command)
{
	if (command == "load") { return eInCmd::LOAD; }
	else if (command == "save") { return eInCmd::SAVE; }
	else if (command == "exit") { return eInCmd::EXIT; }
	else if (command == "eqm") { return eInCmd::EQM; }
	else if (command == "sim") { return eInCmd::SIM; }
	else if (command == "help") { return eInCmd::HELP; }
	else { return eInCmd::INVALID; }
}

void printWelcome()
{
	std::cout << "//////////////////////////////////////////////////////\n";
	std::cout << "//     Zac's Semiconductor Device Sim               //\n";
	std::cout << "//     Author  :  Zachary Humphreys                 //\n";
	std::cout << "//                                                  //\n";
	std::cout << "//     Feel Free to distribute as you please.       //\n";
	std::cout << "//     Completely open source. Do whatever!         //\n";
	std::cout << "//////////////////////////////////////////////////////\n";
}

void printHelp()
{
	std::cout << "____________________________________________________________________________________________" << std::endl;
	std::cout << "Command List:\n";
	std::cout << "\t\"load\" <FILENAME>:\t\t This will load the device state from the filename written after load.\n\t\t\t\t\t   DO NOT INCLUDE \".CSV\" AT END!\n";
	std::cout << "\t \"save\" <FILENAME>:\t\t This will save the device state from the filename written after load.\n\t\t\t\t\t   DO NOT INCLUDE \".CSV\" AT END!\n";
	std::cout << "\t \"eqm\" <TOLERANCE> <MULTIPLIER>: This will bring the loaded device to equilibrium. Will be considered \n\t\t\t\t\t  equilibrium when total currents less than tolerance. Multiplier speeds \n\t\t\t\t\t  calculation up at cost of accuracy.\n";
	std::cout << "\t \"sim\" <FILENAME>:\t\t This will run the simulation on the loaded device.\n";
	std::cout << "\t \"exit\":\t\t\t This will exit the program.\n";
}

void printExit()
{
	std::cout << "Thank you so much for using my program! - Zac\n";
}

void loadFile(std::string FileName)
{
	Device.loadState(FileName);
}

void callEqm()
{
	double Tolerance = 0;
	double ExchangeScale = 0;
	std::cin >> Tolerance >> ExchangeScale;
	Device.bringToEqm(Tolerance, ExchangeScale, true);
}

inline void flushCin()
{
	//Flush buffer of further commands
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return;
}
