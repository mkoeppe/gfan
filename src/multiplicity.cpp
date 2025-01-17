#include "multiplicity.h"
#include "lll.h"
#include "polyhedralcone.h"
#include "wallideal.h"
#include "saturation.h"
#include "linalg.h"
#include "log.h"

IntegerVector writeInTermsOf(IntegerVector const &v, IntegerMatrix const &b)
{//write the row vector v as a integer linear combination of the basis elements of the basis b. Asserts if v is not in the lattice generated by b

  /*  AsciiPrinter Q(Stderr);
    fprintf(Stderr,"Dimensions %ix%i\n",b.getHeight(),b.getWidth());
    Q.printVector(v);
    Q.printVectorList(b.getRows());
  */

  int m=b.getHeight();//dimension of lattice
  assert(v.size()==b.getWidth());
  IntegerVectorList equations=b.getRows();
  equations.push_back(v);
  IntegerVectorList inequalities;
  inequalities.push_back(-IntegerVector::standardVector(m+1,m));

  PolyhedralCone p(inequalities,rowsToIntegerMatrix(equations).transposed().getRows(),m+1);
  IntegerVector w=p.getRelativeInteriorPoint();
  //Q.printVector(w);
  assert(w[m]==-1);
  IntegerVector ret=w.subvector(0,m);
  {
    IntegerMatrix ret2(1,m);
    ret2[0]=ret;
    //fprintf(Stderr,"%i\n",(ret2*b)[0].size());
    assert(((ret2*b)[0]-v).isZero());
  }
  /*    Q.printVector(ret);
	fprintf(Stderr,"Returning!!1\n");*/
  return ret;
}

Polynomial notLaurent(Polynomial p)
{
  if(!p.isZero())
    {
      IntegerVector v=p.terms.begin()->first.exponent;
      for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
	v=min(v,i->first.exponent);

      p*=Monomial(p.getRing(),-v);
    }
  return p;
}

Polynomial multiplicativeChangeInv(Polynomial const &p, IntegerMatrix const &lattice, PolynomialRing const &r2)
{
  PolynomialRing theRing=p.getRing();
  Polynomial ret(r2);
  if(!p.isZero())
    {
      IntegerVector rel=p.terms.begin()->first.exponent;
      for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
	ret+=Term(i->second,Monomial(r2,writeInTermsOf(i->first.exponent-rel,lattice)));
    }
  return ret;
}

/*
//Old implementation
PolynomialSet multiplicativeChangeInv(PolynomialSet const &g, IntegerMatrix const &lattice, PolynomialRing const &r2)
{
  PolynomialRing theRing=g.getRing();
  PolynomialSet ret(r2);
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    ret.push_back(multiplicativeChangeInv(*i,lattice,r2));

  return ret;
}
*/

PolynomialSet multiplicativeChangeInv(PolynomialSet const &g, IntegerMatrix const &lattice, PolynomialRing const &r2)
{
  PolynomialRing theRing=g.getRing();
  PolynomialSet ret(r2);
  FieldMatrix bigMatrix=combineLeftRight(integerMatrixToFieldMatrix(lattice,Q),FieldMatrix::identity(Q,lattice.getHeight()));
  bigMatrix.reduce();
  FieldMatrix reducedLatticeBasis=bigMatrix.submatrix(0,0,lattice.getHeight(),lattice.getWidth());
  FieldMatrix inverseMatrix=bigMatrix.submatrix(0,lattice.getWidth(),lattice.getHeight(),lattice.getWidth()+lattice.getHeight());
  for(PolynomialSet::const_iterator i=g.begin();i!=g.end();i++)
    {
      Polynomial q(r2);
      Polynomial const &p=*i;
      if(!p.isZero())
        {
          IntegerVector rel=p.terms.begin()->first.exponent;
          for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
            {
              IntegerVector diff=i->first.exponent-rel;
              FieldVector diff2=integerVectorToFieldVector(diff,Q);
              FieldVector temp(Q,0);
              reducedLatticeBasis.normalForm(diff2,&temp);
              FieldVector temp2=temp*inverseMatrix;
              q+=Term(i->second,Monomial(r2,fieldVectorToIntegerVector(temp2)));
            }
        }
      ret.push_back(q);
    }
  return ret;
}


PolynomialSet notLaurent(PolynomialSet const &s)
{
  PolynomialRing theRing=s.getRing();
  PolynomialSet ret(theRing);
  for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)
    ret.push_back(notLaurent(*i));

  return ret;
}


PolynomialSet idealWithSameMultiplicity(PolynomialSet const &g)
{
  PolynomialRing r1=g.getRing();
  PolyhedralCone H=homogeneitySpace(g).dualCone();

  H.findFacets();

  IntegerVectorList a=H.getEquations();

  //log3  IntegerMatrix latticeBasis=rowsToIntegerMatrix(a);

  IntegerMatrix temp(H.ambientDimension(),0);
  if(a.size()!=0)temp=rowsToIntegerMatrix(a).transposed();
#if 1
  FieldMatrix temp2=integerMatrixToFieldMatrix(temp,Q).transposed();
//  IntegerMatrix latticeBasis=fieldMatrixToIntegerMatrix(temp2.reduceAndComputeKernel());
  IntegerMatrix latticeBasis=fieldMatrixToIntegerMatrix(temp2.basisOfLatticeKernel());
#else
    IntegerMatrix latticeBasis=latticeKernelOfTransposed(temp);
#endif

  mlll(latticeBasis);

//  log3  AsciiPrinter(Stderr).printVectorList(latticeBasis.getRows());

  PolynomialRing r2(r1.getField(),latticeBasis.getHeight());

  return notLaurent(multiplicativeChangeInv(g,latticeBasis,r2));
}

int multiplicity(PolynomialSet const &g)
{
	PolynomialSet g2=idealWithSameMultiplicity(g);
  log3 AsciiPrinter(Stderr).printPolynomialSet(g2);
  PolynomialSet g3=nonHomogeneousSaturation(g2);

  log3  AsciiPrinter(Stderr).printPolynomialSet(g3);

  return numberOfStandardMonomials(g3);
}


bool isStandard(IntegerVector const &v, PolynomialSet const &markedGroebnerBasis)
{
  for(PolynomialSet::const_iterator i=markedGroebnerBasis.begin();i!=markedGroebnerBasis.end();i++)
    if(i->getMarked().m.exponent.divides(v))return false;
  return true;
}


int numberOfStandardMonomials(PolynomialSet const &markedGroebnerBasis)
{
  int ret=0;
  int n=markedGroebnerBasis.numberOfVariablesInRing();

  IntegerVector v(n);

  int i=0;
  while(1)
    {
      //  AsciiPrinter(Stderr).printVector(v);
      if(isStandard(v,markedGroebnerBasis))
	{
	  //AsciiPrinter(Stderr).printVector(v);
	  //fprintf(Stderr,"\n");
	  ret++;
	  i=n-1;
	}
      else
	{
	  v[i]=0;
	  i--;
	}
      if(i==-1)break;
      v[i]++;
    }
  return ret;
}
