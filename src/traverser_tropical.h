#ifndef TRAVERSER_TROPICAL_H_INCLUDED
#define TRAVERSER_TROPICAL_H_INCLUDED

#include "symmetrictraversal.h"
#include "polynomial.h"

class TropicalTraverser: public ConeTraverser
{
	PolynomialSet coneGroebnerBasis;
	PolynomialSet idealGroebnerBasis;
	PolyhedralCone theCone;
	IntegerVector positiveGrading;

	int n,d;
	void updatePolyhedralCone();
public:
	TropicalTraverser(PolynomialSet const &coneGroebnerBasis_, PolynomialSet const &idealGroebnerBasis_);
	virtual void changeCone(IntegerVector const &ridgeVector, IntegerVector const &rayVector);
	virtual IntegerVectorList link(IntegerVector const &ridgeVector);
	virtual PolyhedralCone & refToPolyhedralCone();
};

#endif
