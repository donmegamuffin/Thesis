#pragma once
//Basically just a container for the data 
//needed at each point in the device
class Node
{
public:

	//Constructor with Ec,Ev terms = 0
	Node(double n, double p, double Nd, double Na);
	//Full constructor
	Node(double n, double p, double Nd, double Na, double Ec, double Ev);

	Node();				//Default constructor creates 0'd out node
	~Node();
	double n;			//density of -ve charge carriers per cm^-3
	double p;			//density of +ve charge carriers per cm^-3
	double Nd;			//density of Donor Ions	(+ve) per cm^-3
	double Na;			//density of Acceptor Ions	(-ve) per cm^-3
	double V;			//Potential at that node
	double Ec;			//Conduction Band Energy
	double Ev;			//Valence Band Energy
private:			
};

