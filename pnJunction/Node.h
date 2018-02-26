#pragma once
//Basically just a container for the data 
//needed at each point in the device
class Node
{
public:
	Node(double n, double p, double Nd, double Na);
	Node();				//Default constructor creates 0'd out node
	~Node();
	double n;			//number of -ve charge carriers
	double p;			//number of +ve charge carriers
	double Nd;			//number of Donor Ions	(+ve)
	double Na;			//number of Acceptor Ions	(-ve)
	double V;			//Potential at that node
private:			
};

