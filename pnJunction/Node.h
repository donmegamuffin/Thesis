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
	double n;			//number of -ve charge carriers
	double p;			//number of +ve charge carriers
	double Nd;			//number of Donor Ions	(+ve)
	double Na;			//number of Acceptor Ions	(-ve)
	double V;			//Potential at that node
	double Ec;			//Conduction Band Energy
	double Ev;			//Valence Band Energy
private:			
};

