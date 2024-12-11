#include "polynomial.h"

class SOSVariable
{
public:
	SOSVariable(class Problem& problem, Polynomial const &halfSupport);
};

class SOSContraint
{
public:
	vector<pair<SOSVariable,Polynomial>>;
	SOSConstrant();
}
class Problem
{
	PolynomialRing r;
	vector<SOSVariables&> sosVariables;
	vector<SOSConstraint> constraints;
public:
	Problem(PolynomialRing const &r_):
		r(r_)
	{
	}
	void addSOSConstraint(SOSConstraint const &constraint)
	{
	}
	void addSOSVariable();
};

