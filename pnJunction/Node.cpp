#include "stdafx.h"
#include "Node.h"


Node::Node()
{
	n = 0;
	p = 0;
	Nd = 0;
	Na = 0;
	V = 0;
}


Node::Node(double in_n, double in_p, double in_Nd, double in_Na)
{
	n = in_n;
	p = in_p;
	Nd = in_Nd;
	Na = in_Na;
	V = 0;
}

Node::~Node()
{
}
