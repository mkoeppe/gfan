#ifndef BUCHBERGER_H_INCLUDED
#define BUCHBERGER_H_INCLUDED

#include "polynomial.h"

/**
 * When autoreducing after lifting a Groebner basis for an initial ideal, the relevant term order is unknown.
 * There are several solutions to this problem:
 * 1) Do autoreduction just with the given markings in a non-systematic way.
 * 2) Represent the relevant term ordering abstractly/symbolically.
 * 3) Compute a weight termordering inducing the marking using Linear Programming.
 * 4) Compute a weight termordering avoiding Linear Programming by providing enough perturbation information.
 *
 * At a high abstraction level, 1) and 2) preferable.
 * At a low implementation level, 3) and 4) are better since monomials can be compared via a few machine integer comparisons.
 * For large
 */
class WeightRecipe
{
	IntegerVector positiveGrading;
	IntegerVector omega;
	IntegerVector perturbationVector;
public:
	WeightRecipe(IntegerVector const &positiveGrading_, IntegerVector const &omega_, IntegerVector const &perturbationVector_):
		positiveGrading(positiveGrading_),
		omega(omega_),
		perturbationVector(perturbationVector_)
	{
	}
	IntegerVector computeWeight(PolynomialSet const &g)const
	{
		bool first=true;
		int64 O1,P1;
		for(auto &f:g)
			for(auto t:f.terms)
			{
				IntegerVector a=f.getMarked().m.exponent-t.first.exponent;
				if(!a.isZero())
				{
					auto O=dotLong(omega,a);
					auto P=dotLong(perturbationVector,a);
					assert(O>=0);
					assert(O>0||P>0);
					if(O>0&&P<0)
					{
						if(first || (O*P1>O1*P)) //possible overflow
						{
							O1=O;
							P1=P;
							first=false;
						}
					}
				}
			}
		if(first)
			return IntegerVector::allOnes(g.getRing().getNumberOfVariables());
		// Epsilon should be in the interval (0,O1/-P1).
		// Let's choose O1/-2P1.
		IntegerVector ret=O1*perturbationVector-2*P1*omega;

		if(!ret.isPositive())
		{
			for(int i=0;i<ret.size();i++)
				if(ret[i]<0)
					ret+=(ret[i]/positiveGrading[i]+1)*positiveGrading;
		}
		assert(ret.isPositive());
		return ret;
	}
};



Polynomial sPolynomial(Polynomial a, Polynomial b);
void buchberger(PolynomialSet *g, TermOrder const &termOrder, bool allowSaturation=false);
void minimize(PolynomialSet *g);
/**
 * The Groebner basis g must be minimal before this function is called.
 */
void autoReduce(PolynomialSet *g, TermOrder const &termOrder, class WeightRecipe const *recipe=0);
static inline void autoReduce_(PolynomialSet *g, TermOrder const &termOrder, class WeightRecipe const *recipe=0){return autoReduce(g,termOrder,recipe);}//<----avoiding scoperules. Should start using namespaces.
bool isMarkedGroebnerBasis(PolynomialSet const &g);

/* For the autoReduction procedure the TermOrder argument is only used as an
argument for the division algorithm. This means that the input to autoReduce
should be marked and that the term order does not have to agree.
*/

/** This routine takes a marked Groebner basis and checks if it is minimal. TODO: there is room for improvement here. */
bool isMinimal(PolynomialSet const &g);
/** This routine takes a marekd Groebner basis and checks if it is autoreduced. This means that it also checks if the basis is minimal.*/
bool isReduced(PolynomialSet const &g);

#endif
