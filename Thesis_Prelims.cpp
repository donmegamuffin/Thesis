// Thesis_Prelims.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>

//Transfer constant

const double TRANSFER_CONST = 0.3;
const int NUMBER_OF_NODES = 10;

class Node
{
public:
	Node();
	Node(int x, double v, double n=0, double p=0);

/*private:*/
	int positionLabel;
	double voltage;
	double electrons;
	double holes;
	
};

Node::Node()
{	//Only use for initialising arrays
	positionLabel = 0;
	voltage = 0;
	electrons = 0;
	holes = 0;
}

Node::Node(int x, double v, double n, double p)
{
	positionLabel = x;
	voltage = v;
	electrons = n;
	holes = p;
}

class NodeInterfaceHandler
{
public:
	NodeInterfaceHandler();
	void VoltageTransfer2(Node& LHS, Node& RHS);
	void ElectronTransfer2(Node& LHS, Node& RHS);
	void HoleTransfer2(Node& LHS, Node& RHS);
};

NodeInterfaceHandler::NodeInterfaceHandler()
{

}

void NodeInterfaceHandler::VoltageTransfer2(Node& LHS, Node& RHS)
{
	//TODO
	//Filler code below
	double Vdiff = LHS.voltage - RHS.voltage;
	//If working on the first node (wire) connection
	if (LHS.positionLabel == 0)
	{
		RHS.voltage += TRANSFER_CONST*Vdiff;

	}
	//If working with last node (wire) connection
	else if(RHS.positionLabel== NUMBER_OF_NODES-1)
	{
		LHS.voltage -= TRANSFER_CONST*Vdiff;
	}
	//between semiconductor parts
	else
	{
		RHS.voltage += TRANSFER_CONST*Vdiff;
		LHS.voltage -= TRANSFER_CONST*Vdiff;
	}
	return;
}

void NodeInterfaceHandler::ElectronTransfer2(Node& LHS, Node& RHS)
{
	//TODO
	//Filler code below
	double eDiff = LHS.electrons - RHS.electrons;
	//If working on the first node (wire) connection
	if (LHS.positionLabel == 0)
	{
		RHS.electrons += TRANSFER_CONST*eDiff;
	}
	//If working with last node (wire) connection
	else if (RHS.positionLabel == NUMBER_OF_NODES-1)
	{
		LHS.electrons -= TRANSFER_CONST*eDiff;
	}
	//between semiconductor parts
	else
	{
		RHS.electrons += TRANSFER_CONST*eDiff;
		LHS.electrons -= TRANSFER_CONST*eDiff;
	}
	return;
}

void NodeInterfaceHandler::HoleTransfer2(Node& LHS, Node& RHS)
{
	//TODO
	//Filler code below
	double pDiff = LHS.holes - RHS.holes;
	//If working on the first node (wire) connection
	if (LHS.positionLabel == 0)
	{
		RHS.holes += TRANSFER_CONST*pDiff;
	}
	//If working with last node (wire) connection
	else if (RHS.positionLabel == NUMBER_OF_NODES-1)
	{
		LHS.holes -= TRANSFER_CONST*pDiff;
	}
	//between semiconductor parts
	else
	{
		RHS.holes += TRANSFER_CONST*pDiff;
		LHS.holes -= TRANSFER_CONST*pDiff;
	}
	return;
}

void printNodes(Node* nArray, int arrayLength)
{
	
	for (int n = 0; n < arrayLength; n++)
	{		
		std::cout << nArray[n].positionLabel
			<< "," << std::setprecision(5) << nArray[n].voltage
			<< "," << nArray[n].electrons
			<< "," << nArray[n].holes
			<< std::endl;			
	}
}

int main()
{
	Node nodes[NUMBER_OF_NODES];
	NodeInterfaceHandler nih;

	//SetData for +V wire
	nodes[0].positionLabel = 0;
	nodes[0].voltage = 4;
	nodes[0].electrons = 0;
	nodes[0].holes = 1e10;

	//Set data for ground wire
	nodes[NUMBER_OF_NODES-1].positionLabel = NUMBER_OF_NODES-1;
	nodes[NUMBER_OF_NODES-1].voltage = 0;
	nodes[NUMBER_OF_NODES-1].electrons = 1e10;
	nodes[NUMBER_OF_NODES-1].holes = 1e10;

	//Loop around entire system
	//Set the position vals
	for (int i = 1; i < (NUMBER_OF_NODES-1); i++)
	{
		nodes[i].positionLabel = i;
	}

	//Output to system
	std::cout << "x,V,n,p" << std::endl;
	for (int i = 0; i < 100; i++)
	{
		printNodes(nodes,NUMBER_OF_NODES);
		for (int n = 1; n < NUMBER_OF_NODES; n++)
		{
			nih.VoltageTransfer2(nodes[n - 1], nodes[n]);
			nih.ElectronTransfer2(nodes[n - 1], nodes[n]);
			nih.HoleTransfer2(nodes[n - 1], nodes[n]);			
		}
	}
    return 0;
}