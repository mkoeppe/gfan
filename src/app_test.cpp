#include <iostream>
#include <stdlib.h>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "polyhedralcone.h"
#include "gfanapplication.h"
#include "saturation.h"
#include "field_rationals.h"
#include "field_zmodpz.h"
#include "field_rationalfunctions.h"
#include "symmetry.h"
#include "linalg.h"
#include "fieldlp.h"
#include "integer.h"
#include "polynomialgcd.h"
#include "packedmonomial.h"
#include "gfanlib_zcone.h"
#include "gfanlib_tableau.h"
#include "gfanlib_hypersurfaceintersection.h"
#include "gfanlib_circuittableint.h"
#include "gfanlib_mixedvolume.h"
#include "lll.h"
#include "gfanlib_circuittableinteger.h"

using namespace gfan;
/*using gfan::CircuitTableInt32;
using gfan::IntMatrix;
using gfan::Cone;
using namespace gfan::MixedVolumeExamples;
*/
//template class gfan::Vector<gfan::Integer>;
//template class gfan::Vector<gfan::Rational>;
//template class gfan::Matrix<gfan::Integer>;
//template class gfan::Matrix<gfan::Rational>;

class TestApplication : public GFanApplication
{
  StringOption testSuiteFolderOption;
  StringOption executableOption;
  StringOption developerTestOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This runs the test suite and checks against the stored result. If no result exists, it is generated.\n";
  }
  TestApplication():
	testSuiteFolderOption("--suite","Specify the folder which contains the test suite.","testsuite"),
	executableOption("--gfan","Specify name of gfan executable to test.","./gfan"),
	developerTestOption("--dev", "Specify name of particular kind of test.","testsuite")
  {
	  developerTestOption.hide();
	  registerOptions();
  }

  const char *name()
  {
    return "_test";
  }
  int testGCD()
  {
	  PolynomialRing r=StringParser("Q[p126,p125,p124,p123,p027,p026,p025,p024,p023,p017,p016,p015,p014,p013,p012]").parsePolynomialRing();

	  // Crashes (when linked with Singular)
//	  Polynomial p=StringParser("-p126*p027*p026*p015^2+p125*p026*p025*p017*p016-p125*p026^2*p017*p015-p125*p027*p025*p016^2+p125*p027*p026*p016*p015-p126*p025^2*p017*p016+p126*p026*p025*p017*p015+p126*p027*p025*p016*p015").parsePolynomial(r);
//	  Polynomial q=StringParser("p126*p027*p026*p025*p016*p015*p014-p124*p026*p025^2*p017*p016^2+2*p124*p026^2*p025*p017*p016*p015-p124*p026^3*p017*p015^2+p124*p027*p025^2*p016^3-2*p124*p027*p026*p025*p016^2*p015+p124*p027*p026^2*p016*p015^2+p125*p026*p025*p024*p017*p016^2-p125*p026^2*p024*p017*p016*p015-p125*p026^2*p025*p017*p016*p014+p125*p026^3*p017*p015*p014-p125*p027*p025*p024*p016^3+p125*p027*p026*p024*p016^2*p015+p125*p027*p026*p025*p016^2*p014-p125*p027*p026^2*p016*p015*p014-p126*p026*p025*p024*p017*p016*p015+p126*p026*p025^2*p017*p016*p014+p126*p026^2*p024*p017*p015^2-p126*p026^2*p025*p017*p015*p014+p126*p027*p025*p024*p016^2*p015-p126*p027*p025^2*p016^2*p014-p126*p027*p026*p024*p016*p015^2").parsePolynomial(r);

	  // Crashes (even when not linked to Singular)
	  Polynomial p=StringParser("-p126*p027*p026*p015^2*p014+p124*p025^2*p017*p016^2-2*p124*p026*p025*p017*p016*p015+p124*p026^2*p017*p015^2+p125*p026*p025*p017*p016*p014-p125*p026^2*p017*p015*p014-p125*p027*p025*p016^2*p014+p125*p027*p026*p016*p015*p014-p126*p025^2*p017*p016*p014+p126*p026*p025*p017*p015*p014+p126*p027*p025*p016*p015*p014").parsePolynomial(r);
	  Polynomial q=StringParser("p126*p027*p026*p025*p016*p015*p014-p124*p026*p025^2*p017*p016^2+2*p124*p026^2*p025*p017*p016*p015-p124*p026^3*p017*p015^2+p124*p027*p025^2*p016^3-2*p124*p027*p026*p025*p016^2*p015+p124*p027*p026^2*p016*p015^2+p125*p026*p025*p024*p017*p016^2-p125*p026^2*p024*p017*p016*p015-p125*p026^2*p025*p017*p016*p014+p125*p026^3*p017*p015*p014-p125*p027*p025*p024*p016^3+p125*p027*p026*p024*p016^2*p015+p125*p027*p026*p025*p016^2*p014-p125*p027*p026^2*p016*p015*p014-p126*p026*p025*p024*p017*p016*p015+p126*p026*p025^2*p017*p016*p014+p126*p026^2*p024*p017*p015^2-p126*p026^2*p025*p017*p015*p014+p126*p027*p025*p024*p016^2*p015-p126*p027*p025^2*p016^2*p014-p126*p027*p026*p024*p016*p015^2").parsePolynomial(r);


	  // Finds factor
//	  Polynomial p=StringParser("-p027*p026*p016*p015-p026*p025*p017*p016+p026^2*p017*p015+p027*p025*p016^2").parsePolynomial(r);
//	  Polynomial q=StringParser("p126*p027*p025*p016*p014+p124*p026*p025*p017*p016-p124*p026^2*p017*p015-p124*p027*p025*p016^2+p124*p027*p026*p016*p015-p125*p026*p024*p017*p016+p125*p026^2*p017*p014+p125*p027*p024*p016^2-p125*p027*p026*p016*p014+p126*p026*p024*p017*p015-p126*p026*p025*p017*p014-p126*p027*p024*p016*p015").parsePolynomial(r);

	  // Finds only 1 as factor
	  //	  Polynomial p=StringParser("-p026*p015+p025*p016").parsePolynomial(r);
//	  Polynomial q=StringParser("p126*p027*p025*p016*p014+p124*p026*p025*p017*p016-p124*p026^2*p017*p015-p124*p027*p025*p016^2+p124*p027*p026*p016*p015-p125*p026*p024*p017*p016+p125*p026^2*p017*p014+p125*p027*p024*p016^2-p125*p027*p026*p016*p014+p126*p026*p024*p017*p015-p126*p026*p025*p017*p014-p126*p027*p024*p016*p015").parsePolynomial(r);

	  Polynomial g=polynomialGCD(p,q);
	  debug<<g;
	  assert(0);
	  return 0;
  }
int testIntegers()
{
/*	Integer a(1);
	Integer b=1;
	for(int i=0;i<300;i++)
	{


		Integer c=a;
	c+=b;
		b=a;
		a=c;
		cerr<<a<<endl;
	}*/
	return 0;
}
gfan::Matrix<CircuitTableInt32> randomMatrix(int h, int w)
{
	gfan::Matrix<CircuitTableInt32> ret(h,w);
	for(int i=0;i<h;i++)
		for(int j=0;j<w;j++)
			ret[i][j]=rand()%5;
	for(int j=0;j<w;j++)
		ret[ret.getHeight()-1][j]=1;
	return ret;
}
gfan::Matrix<CircuitTableInt32> randomMatrix2(int h, int w)
{
	gfan::Matrix<CircuitTableInt32> ret(h,w);
	for(int i=0;i<h;i++)
		for(int j=0;j<w;j++)
			ret[i][j]=rand()%5-2;
	return ret;
}

template<typename typ>
gfan::Matrix<typ> convertMatrix(IntMatrix const &m)
{
	gfan::Matrix<typ> ret(m.getHeight(),m.getWidth());
	for(int i=0;i<m.getHeight();i++)
		for(int j=0;j<m.getWidth();j++)
			ret[i][j]=typ(m[i][j]);
  return ret;
}

vector<gfan::Matrix<CircuitTableInt32> > convertMatrixVectorT(vector<IntMatrix> const &v)
{
	vector<gfan::Matrix<CircuitTableInt32> > ret;
	for(int i=0;i<v.size();i++)
	{
		ret.push_back(convertMatrix<CircuitTableInt32>(v[i]).transposed());
	}
	return ret;
}

IntMatrix rowsToIntegerMatrix2(IntegerVectorList const &l)//taken from app_mixedvolume.cpp
{
	  assert(l.size());
	  IntMatrix ret(l.size(),l.front().size());
	  IntegerVectorList::const_iterator I=l.begin();
	  for(int i=0;i<ret.getHeight();i++,I++)
		  for(int j=0;j<ret.getWidth();j++)
			  ret[i][j]=(*I)[j];
	  return ret;
}

/*template<typename typ>
int numberOfEdges(gfan::Matrix<typ> const &A)
{
	int ret=0;
	for(int i=0;i<A.getHeight();i++)
		for(int j=0;j<i;j++)
		{
			gfan::Matrix<typ> B(0,A.getWidth());
			for(int k=0;k<A.getHeight();k++)B.appendRow(A[i].toVector()-A[k].toVector());
			for(int k=0;k<A.getHeight();k++)B.appendRow(A[j].toVector()-A[k].toVector());
			auto C=Cone<typ>(B.transposed());
			ret+=C.getDimension()==C.getAmbientDimension()-1;
		}
	std::cerr<<ret<<"\n";
	return ret;
}*/


IntegerVector getInitialExponent(Polynomial p, TermOrder const &order)
{
	assert(!p.isZero());
	p.mark(order);
	return p.getMarked().m.exponent;
}
Polynomial SP(Polynomial a, Polynomial b, TermOrder const &order)
{
	a.mark(order);
	b.mark(order);
	return sPolynomial(a,b);
}

PolynomialSet testGBIdea(PolynomialSet start, TermOrder const &startOrder, TermOrder const &targetOrder)
{
	start.mark(startOrder);
	PolynomialSet ret=start;
	vector<Polynomial> sPolys;
	for(PolynomialSet::const_iterator i=ret.begin();i!=ret.end();i++)
		for(PolynomialSet::const_iterator j=ret.begin();j!=i;j++)
		{
			if(!i->checkMarking(targetOrder)||!j->checkMarking(targetOrder))
			{
				if(!relativelyPrime(getInitialExponent(*j,targetOrder),getInitialExponent(*i,targetOrder)))
					sPolys.push_back(SP(*i,*j,targetOrder));
			}
		}
	ret.mark(targetOrder);
	while(!sPolys.empty())
	{
		AsciiPrinter P(Stderr);
		std::cerr<<"ret:\n";	{P.printPolynomialSet(ret);}
		std::cerr<<"SPO-----------------\n"<<sPolys.size()<<"\n";for(auto &p:sPolys){P.printPolynomial(p);std::cerr<<"\n";}

		if(sPolys.size()>200){std::cerr<<"SKiptest\n";return PolynomialSet(start.getRing());}

		Polynomial p=sPolys.back();sPolys.pop_back();
		auto n=division(p,ret,targetOrder);
		n.mark(targetOrder);
		std::cerr<<"Remainder";P.printPolynomial(n);std::cerr<<"\n";
		if(!n.isZero())
		{
			if(!n.checkMarking(startOrder))
			{
				std::cerr<<"Orders do not agree so we produce all spolys\n";
				//add all spolys
				for(auto &f:ret)
//					if(!relativelyPrime(getInitialExponent(n,targetOrder),getInitialExponent(f,targetOrder)))
					sPolys.push_back(SP(f,n,targetOrder));
			}
			else
			{
				std::cerr<<"Orders agree so we produce only some s-polys \n";
				for(auto &f:ret)
					if(!f.checkMarking(startOrder))
//						if(!relativelyPrime(getInitialExponent(n,targetOrder),getInitialExponent(f,targetOrder)))
						sPolys.push_back(SP(f,n,targetOrder));
				//add some spolys
			}
			ret.push_front(n);
		}
	}
	return ret;
}


int testGfanLib()
	{
	gfan::Integer temp((signed long int)0);
#if 0
	if(1)
		{
	/*		  1  1  1 0
			 -1 -1 -1 0
			  0  1 -1 0
			  0 -1  1 0
			  1  0 -1 1
			  1 -1  0 0
			 -1  1  0 0
			  1  0 -1 1*/

			  gfan::Matrix<CircuitTableInt32> A(9,4);
			  int i=0;
			  A[i][0]= 1;A[i][1]= 1;A[i][2]= 1;A[i][3]= 0;i++;
			  A[i][0]=-1;A[i][1]=-1;A[i][2]=-1;A[i][3]= 0;i++;
			  A[i][0]= 0;A[i][1]= 1;A[i][2]=-1;A[i][3]= 0;i++;
			  A[i][0]= 0;A[i][1]=-1;A[i][2]= 1;A[i][3]= 0;i++;
			  A[i][0]= 1;A[i][1]= 0;A[i][2]=-1;A[i][3]= 1;i++;
			  A[i][0]= 1;A[i][1]=-1;A[i][2]= 0;A[i][3]= 0;i++;
			  A[i][0]=-1;A[i][1]= 1;A[i][2]= 0;A[i][3]= 0;i++;
			  A[i][0]= 1;A[i][1]= 0;A[i][2]=-1;A[i][3]= 1;i++;
			  A[i][0]= 0;A[i][1]= 0;A[i][2]= 0;A[i][3]=-1;i++;
				Cone<CircuitTableInt32> C1(A);
				std::cerr<<C1.toString();
			return 0;
		}
#endif
	if(0)
	{
/*		  1  1  1 0
		 -1 -1 -1 0
		  0  1 -1 0
		  0 -1  1 0
		  1  0 -1 1
		  1 -1  0 0
		 -1  1  0 0
		  1  0 -1 1*/

		  gfan::Matrix<CircuitTableInt32> A(9,4);
		  int i=0;
		  A[i][0]= 1;A[i][1]= 1;A[i][2]= 1;A[i][3]= 0;i++;
		  A[i][0]=-1;A[i][1]=-1;A[i][2]=-1;A[i][3]= 0;i++;
		  A[i][0]= 0;A[i][1]= 1;A[i][2]=-1;A[i][3]= 0;i++;
		  A[i][0]= 0;A[i][1]=-1;A[i][2]= 1;A[i][3]= 0;i++;
		  A[i][0]= 1;A[i][1]= 0;A[i][2]=-1;A[i][3]= 1;i++;
		  A[i][0]= 1;A[i][1]=-1;A[i][2]= 0;A[i][3]= 0;i++;
		  A[i][0]=-1;A[i][1]= 1;A[i][2]= 0;A[i][3]= 0;i++;
		  A[i][0]= 1;A[i][1]= 0;A[i][2]=-1;A[i][3]= 1;i++;
		  A[i][0]= 0;A[i][1]= 0;A[i][2]= 0;A[i][3]=-1;i++;
			Cone<CircuitTableInt32> C1(A.transposed());
			std::cerr<<C1.toString();
		return 0;
	}
	if(0){
		IntegerVectorList temp=StringParser(
		//		  "{(1,1,1,1,0),"
//				  "(-1,-1,-1,-1,0),"
//				  "(0,1,0,-1,0),"
//				  "(0,0,1,-1,0),"
				  "{(1,0,0,-1,0),"
				  "(-1,0,0,1,0),"
//				  "(0,1,0,-1,0),"
				  "(1,0,-1,0,0),"
				  "(-1,0,1,0,0),"
//				  "(0,1,-1,0,0),"
				  "(0,1,0,-1,0),"
				  "(0,-1,0,1,0),"
				  "(-1,1,0,0,1),"
				  "(0,0,0,0,-1)}"
					).parseIntegerVectorList();
		  gfan::Matrix<CircuitTableInt32> A=convertMatrix<CircuitTableInt32>(rowsToIntegerMatrix2(temp));
		  Cone<CircuitTableInt32> C1(A.transposed());
		  std::cerr<<C1.toString();
		  IntegerVectorList empty;
		  std::cerr<<"DDD"<<PolyhedralCone(temp,empty).dimension()<<"\n";
		return 0;
	}
	if(0){
		IntegerVectorList temp=StringParser(
				"{(0,0,1,0,1,0,0),"
				"(0,1,0,0,0,1,0),"
				"(1,0,1,0,0,-1,0),"
				"(0,0,-1,1,0,1,0),"
				"(0,0,0,0,0,0,1)}"
					).parseIntegerVectorList();
		  gfan::Matrix<CircuitTableInt32> A=convertMatrix<CircuitTableInt32>(rowsToIntegerMatrix2(temp));
		  Cone<CircuitTableInt32> C1(A);

		  auto linealitySpace=C1.getLinealitySpace();
		  auto rays=C1.getRays();
		return 0;
	}
	if(0)
	{
		auto v=StringParser("(13,5,4)").parseIntegerVector();
		WeightReverseLexicographicTermOrder T(v);
	    FileParser P(Stdin);
	    auto L=P.parsePolynomialSetListWithRing();
		for(auto &g:L)
			if(g.size()<=3)
		{
			std::cerr<<"Trying:\n";
			AsciiPrinter(Stderr).printPolynomialSet(g,true);
			PolynomialSet h(g.getRing());
			for(auto &f:g)
			{
				auto M=f.getMarked().m;
				Polynomial F(g.getRing());
				for(auto &t:f.terms)
					if(!T(t.first.exponent,M.exponent))
						F=F+Polynomial(Term(t.second,t.first));
				h.push_back(F);
			}
			std::cerr<<"Truncation:\n";
			AsciiPrinter(Stderr).printPolynomialSet(h,true);
			buchberger(&h,T);
			minimize(&h);
			autoReduce(&h,T);
			std::cerr<<h.size()<<"\n";
		}
		return 0;
	}

	if(0)
	{
	    FileParser P(Stdin);
	    auto L=P.parsePolynomialSetListWithRing();

	    for(auto &a:L)
	    	if(a.size()<7)
	    {
	        IntegerVectorList normals=wallInequalities(a);
	        PolyhedralCone C(normals,IntegerVectorList());
	        auto v=intersection(PolyhedralCone::positiveOrthant(C.ambientDimension()),C).getRelativeInteriorPoint();
	        AsciiPrinter P(Stderr);
	        P.printVector(v);

	        WeightTermOrder T(v);
	        for(auto &b:L)
	        	if(b.size()<7)
	        {
		        IntegerVectorList normals=wallInequalities(b);
		        PolyhedralCone C(normals,IntegerVectorList());
		        auto u=intersection(PolyhedralCone::positiveOrthant(C.ambientDimension()),C).getRelativeInteriorPoint();
		        WeightTermOrder S(u);
	        	auto c=testGBIdea(b,S,T);
//	        	std::cerr<<"size"<<c.size()<<"\n";
	        	if(c.size())
	        	{

	        	auto aa=a.markedTermIdeal();
	        	auto cc=c.markedTermIdeal();
	        	minimize(&aa);
	        	minimize(&cc);
	        	if(aa.size()!=cc.size())
	        	{
	        		P.printString("Target:");P.printVector(v);
	        		P.printPolynomialSet(a);
	        		P.printPolynomialSet(aa);
	        		P.printString("Start:");P.printVector(u);
	        		P.printPolynomialSet(b);
	        		P.printString("Computed target monomial ideal:");
	        		P.printPolynomialSet(cc);
	        	//	P.printPolynomialSet(c);



	        	}
	        	assert(aa.size()==cc.size());
	        	}
	        	std::cerr<<"OK\n";
	        }
	    }

		return 0;
	}
	if(0) // 0 means cyclic examples
	{if(0)
		{
			IntegerVectorList temp=StringParser(
			"{(0, 1, 0,0, 1,-1, 1, 1, 1)"
			",(0, 1, 1,0, 0, 0, 1, 0, 1)"
			",(1,-1, 1,0, 0, 0, 0, 0, 1)"
			",(1,-1, 0,0, 0, 0, 0, 0,-1)"
			",(0, 0, 0,0,-1, 0,-1, 0,-1)"
			",(-1,0,-1,0, 0, 0,-1,-1,-1)"
			",(-1,0,-1,0, 0, 0, 0, 0, 0)}").parseIntegerVectorList();
			  gfan::Matrix<CircuitTableInt32> A=convertMatrix<CircuitTableInt32>(rowsToIntegerMatrix2(temp));
			  Cone<CircuitTableInt32> C1(A.submatrixColumns(0,6),A.submatrixColumns(6,9));

//			  std::cerr<<"DDD"<<PolyhedralCone(temp,empty).dimension()<<"\n";
			  std::cerr<<"RAYS"<<C1.getRays().toString();
			  std::cerr<<"The vecotr (-1,-2,-5,-5,-2,-1,-100) is in the cone.\n";
//			  C1.extractFaceComplex();
			  return 0;
		}

	if(1){
		std::stringstream s("((a, b),   ((c)f(d)) ,(e)  )");
		auto v=parseSequence(s);
		for(auto e:v)
		{
			std::cerr<<"\""<<e<<"\"\n\n";
		}

		{
			gfan::Vector<CircuitTableInt32> v;
			std::stringstream s("(1,2,3,43)");
			s>>v;
			std::cerr<<v<<"\n";
		}
		{
			std::stringstream s("((1,2,3),(43,2,2))""((1,2,3),(43,2,2))");
			std::cerr<<gfan::Matrix<CircuitTableInt64>::readMatrix(s,3)<<"\n";
			//std::cerr<<gfan::Matrix<CircuitTableInteger>::readMatrix(s,3)<<"\n";
			//		std::cerr<<gfan::Matrix<CircuitTableInt32>::readMatrix(s,3)<<"\n";
			gfan::GeneratedCone<CircuitTableInt32> C(gfan::Matrix<CircuitTableInt32>::readMatrix(s,3));

			std::stringstream A;
			C.save(A);
			auto As=A.str();
			std::cerr<<"AS"<<As<<"\n";
			std::stringstream B(As);
			auto C2=gfan::GeneratedCone<CircuitTableInt32>::load(B);
			C2.save(std::cerr);
		}
		return 0;
	}

//	typedef gfan::CircuitTableInteger typ;
	typedef gfan::CircuitTableInt32 typ;
//	typedef gfan::CircuitTableInt64 typ;
		IntegerVectorList temp=StringParser(
				  "{(1,1,1),"
				  "(1,1,-1),"
				  "(1,-1,1)}"
					).parseIntegerVectorList();
		  gfan::Matrix<typ> A=convertMatrix<typ>(rowsToIntegerMatrix2(temp));
			IntegerVectorList temp2=StringParser(
					  "{(1,-1,-1)}"
						).parseIntegerVectorList();
			  gfan::Matrix<typ> A2=convertMatrix<typ>(rowsToIntegerMatrix2(temp2));
		  HalfOpenCone<typ> C1(A.transposed(),gfan::Matrix<typ>(3,0),A2.transposed());
		  std::cerr<<C1.toString();
		  IntegerVectorList empty;
		  std::cerr<<"DDD"<<PolyhedralCone(temp,empty).dimension()<<"\n";
		  std::cerr<<"RAYS"<<C1.lifted.getRays().toString();
		  std::cerr<<extractFaceComplex(C1);
		  if(0)
		  {
			  gfan::CircuitTableInt32 A(10);
			  std::cerr<<gfan::CircuitTableInteger(A).toString();
			  gfan::CircuitTableInteger B=gfan::CircuitTableInteger(A);
			  for(int i=0;i<3/*+1*/;i++)B=B*B;
			  auto C=gfan::CircuitTableInt32(B);
			  std::cerr<<C.toString();

			  while(1)
			  {
				  gfan::Matrix<gfan::CircuitTableInteger> M(1,1);
				  M[0][0]=B;
				  std::cerr<<M.toString()<<"\n";
				  auto N=gfan::Matrix<gfan::CircuitTableInt32>(M);
				  std::cerr<<N.toString()<<"\n";
				  auto O=gfan::Matrix<gfan::CircuitTableInteger>(N);
				  std::cerr<<O.toString()<<"\n";
				  B=B*B;

			  }
			  return 1;
		  }
		  if(0)
		  {// This code causes overflow with lex rule, but not Bland rule for CircuitTableInt32?
			  auto A=convertMatrix<gfan::CircuitTableInt64>(rowsToIntegerMatrix2(StringParser("{(0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,-2,-2,-2,-2,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0),"
					  "(0,-2,-2,-2,-2,2,2,2,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1),"
					  "(0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,0,-2,-2,0,2,2,0,0,0,0,1,0,0,-1,-1,0,-1,-1,-1),"
					  "(0,0,-2,0,-2,2,0,2,0,0,0,0,1,0,-1,0,-1,-1,-1,-1),"
					  "(0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,0,0,-2,-2,0,2,2,0,0,0,0,0,1,0,-1,-1,-1,-1,-1),"
					  "(0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0),"
					  "(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1),"
					  "(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1),"
					  "(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0),"
					  "(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0),"
					  "(0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0),"
					  "(0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0),"
					  "(0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0),"
					  "(0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0),"
					  "(0,0,0,2,2,0,-2,-2,0,-1,-1,-1,-1,-1,0,0,0,0,0,1),"
					  "(0,0,2,0,2,-2,0,-2,-1,0,-1,-1,-1,-1,0,0,0,0,1,0),"
					  "(0,0,2,2,0,-2,-2,0,-1,-1,0,-1,-1,-1,0,0,0,1,0,0),"
					  "(0,8,8,8,8,0,0,0,-4,-4,-4,5,8,8,12,12,0,9,0,0),"
					  "(0,8,8,8,8,0,0,0,-4,-4,-4,8,5,8,12,0,12,0,9,0),"
					  "(0,8,8,8,8,0,0,0,-4,-4,-4,8,8,5,0,12,12,0,0,9),"
					  "(0,8,8,8,8,0,0,0,-4,5,8,-4,-4,8,12,9,0,12,0,0),"
					  "(0,8,8,8,8,0,0,0,-4,8,5,-4,-4,8,12,0,9,0,12,0),"
					  "(0,8,8,8,8,0,0,0,-4,8,8,-4,-4,5,0,0,0,12,12,9),"
					  "(0,8,8,8,8,0,0,0,5,-4,8,-4,8,-4,9,12,0,12,0,0),"
					  "(0,8,8,8,8,0,0,0,5,8,-4,8,-4,-4,9,0,12,0,12,0),"
					  "(0,8,8,8,8,0,0,0,8,-4,5,-4,8,-4,0,12,9,0,0,12),"
					  "(0,8,8,8,8,0,0,0,8,-4,8,-4,5,-4,0,0,0,12,9,12),"
					  "(0,8,8,8,8,0,0,0,8,5,-4,8,-4,-4,0,9,12,0,0,12),"
					  "(0,8,8,8,8,0,0,0,8,8,-4,5,-4,-4,0,0,0,9,12,12)}").parseIntegerVectorList()));
			  GeneratedCone<gfan::CircuitTableInt64> C1(A.transposed(),gfan::Matrix<gfan::CircuitTableInt64>(A.getWidth(),0));
			  std::cerr<<"DIM"<<C1.getDimension()<<"LINEALITY"<<C1.getLinealitySpace()<<"\n";
		  }


		  {//fan example
			  auto g=StringParser(
#if 0
					  "Q[a,b,c,d,e,f,g,h,i,j]"
					  "{a+b+c+d+e+f+g+h+i+j,"
					  "ab+bc+cd+de+ef+fg+gh+hi+ij+ja,"
					  "abc+bcd+cde+def+efg+fgh+ghi+hij+ija+jab,"
					  "abcd+bcde+cdef+defg+efgh+fghi+ghij+hija+ijab+jabc,"
					  "abcde+bcdef+cdefg+defgh+efghi+fghij+ghija+hijab+ijabc+jabcd,"
					  "abcdef+bcdefg+cdefgh+defghi+efghij+fghija+ghijab+hijabc+ijabcd+jabcde,"
					  "abcdefg+bcdefgh+cdefghi+defghij+efghija+fghijab+ghijabc+hijabcd+ijabcde+jabcdef,"
					  "abcdefgh+bcdefghi+cdefghij+defghija+efghijab+fghijabc+ghijabcd+hijabcde+ijabcdef"
					  "+jabcdefg,"
					  "abcdefghi+abcdefghj+abcdefgij+abcdefhij+abcdeghij+abcdfghij+abcefghij+abdefghij+"
					  "acdefghij+bcdefghij,"
					  "abcdefghij-1}"
#endif
					  "Q[x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20]"
					        "{1-x02*x09^2+x03*x09^2,"
					        "1-x02*x10^2+x04*x10^2,"
					        "1-x02*x11^2+x05*x11^2,"
					        "1-x03*x12^2+x04*x12^2,"
					        "1-x03*x13^2+x05*x13^2,"
					        "1-x04*x14^2+x05*x14^2,"
					        "1-x15^2+x06*x15^2,"
					        "1-x16^2+x07*x16^2,"
					        "1-x17^2+x08*x17^2,"
					        "1-x06*x18^2+x07*x18^2,"
					        "1-x06*x19^2+x08*x19^2,"
					        "1-x07*x20^2+x08*x20^2,"
					        "x02-x01^3*x09*x15^3-x01^5*x10*x16^3-x01^7*x11*x17^3,"
					        "x03+x01^2*x09*x15^3-x01^5*x12*x18^3-x01^7*x13*x19^3,"
					        "x04+x01^2*x10*x16^3+x01^3*x12*x18^3-x01^7*x14*x20^3,"
					        "x05+x01^2*x11*x17^3+x01^3*x13*x19^3+x01^5*x14*x20^3,"
					        "1-x01^3*x09^3*x15-x01^5*x10^3*x16-x01^7*x11^3*x17,"
					        "x06+x01^2*x09^3*x15-x01^5*x12^3*x18-x01^7*x13^3*x19,"
					        "x07+x01^2*x10^3*x16+x01^3*x12^3*x18-x01^7*x14^3*x20,"
					        "x08+x01^2*x11^3*x17+x01^3*x13^3*x19+x01^5*x14^3*x20},"////////
					        "x01^2*x02+x01^3*x03+x01^5*x04+x01^7*x05,"
					        "x01^2+x01^3*x06+x01^5*x07+x01^7*x08"        "},"
					        "-x21+x01^2*x02+x01^3*x03*x06+x01^5*x04*x07+x01^7*x05*x08,"
					        "x01^5*x10^2*x11^2*x12^2*x13^2*x14^2*x16^2*x17^2*x18^2*x19^2*x20^2+x01^7*x09^2*x11^2*x12^2*x13^2*x14^2*x15^2*x17^2*x18^2*x19^2*x20^2-x01^2*x09^2*x10^2*x11^2*x12^2*x13^2*x14^2*x15^2*x16^2*x17^2*x18^2*x19^2*x20^2*x21+x01^8*x09^2*x10^2*x11^2*x13^2*x14^2*x15^2*x16^2*x17^2*x19^2*x20^2-x01^3*x09^2*x10^2*x11^2*x12^2*x13^2*x14^2*x15^2*x16^2*x17^2*x18^2*x19^2*x20^2*x21+x01^9*x09^2*x10^2*x12^2*x13^2*x14^2*x15^2*x16^2*x18^2*x19^2*x20^2+x01^10*x09^2*x10^2*x11^2*x12^2*x14^2*x15^2*x16^2*x17^2*x18^2*x20^2-x01^5*x09^2*x10^2*x11^2*x12^2*x13^2*x14^2*x15^2*x16^2*x17^2*x18^2*x19^2*x20^2*x21+x01^12*x09^2*x10^2*x11^2*x12^2*x13^2*x15^2*x16^2*x17^2*x18^2*x19^2-x01^7*x09^2*x10^2*x11^2*x12^2*x13^2*x14^2*x15^2*x16^2*x17^2*x18^2*x19^2*x20^2*x21,"
					        "-x21+x01^5*x09*x15+x01^7*x10*x16+x01^8*x12*x18+x01^9*x11*x17+x01^10*x13*x19+x01^12*x14*x20}"
					  "Q[x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13]"
      "{1-x02*x07^2+x03*x07^2,"
      "1-x02*x08^2+x04*x08^2,"
      "1-x03*x09^2+x04*x09^2,"
      "1-x10^2+x05*x10^2,"
      "1-x11^2+x06*x11^2,"
      "1-x05*x12^2+x06*x12^2,"
      "x02-x01^3*x07*x10^3-x01^5*x08*x11^3,"
      "x03+x01^2*x07*x10^3-x01^5*x09*x12^3,"
      "x04+x01^2*x08*x11^3+x01^3*x09*x12^3,"
      "1-x01^3*x07^3*x10-x01^5*x08^3*x11,"
      "x05+x01^2*x07^3*x10-x01^5*x09^3*x12,"
      "x06+x01^2*x08^3*x11+x01^3*x09^3*x12,"
      "x01^2*x02+x01^3*x03+x01^5*x04,"
      "x01^2+x01^3*x05+x01^5*x06,"
      "-x13+x01^2*x02+x01^3*x03*x05+x01^5*x04*x06,"
      "x01^5*x08^2*x09^2*x11^2*x12^2+x01^7*x07^2*x09^2*x10^2*x12^2-x01^2*x07^2*x08^2*x09^2*x10^2*x11^2*x12^2*x13+x01^8*x07^2*x08^2*x10^2*x11^2-x01^3*x07^2*x08^2*x09^2*x10^2*x11^2*x12^2*x13-x01^5*x07^2*x08^2*x09^2*x10^2*x11^2*x12^2*x13,"
      "-x13+x01^5*x07*x10+x01^7*x08*x11+x01^8*x09*x12}"


					  "Q[x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22]"
					  "{1-x02*x10^2+x03*x10^2,"
					  "1-x02*x11^2+x04*x11^2,"
					  "1-x02*x12^2+x05*x12^2,"
					  "1-x03*x13^2+x04*x13^2,"
					  "1-x03*x14^2+x05*x14^2,"
					  "1-x04*x15^2+x05*x15^2,"
					  "1-x06*x16^2+x07*x16^2,"
					  "1-x06*x17^2+x08*x17^2,"
					  "1-x06*x18^2+x09*x18^2,"
					  "1-x07*x19^2+x08*x19^2,"
					  "1-x07*x20^2+x09*x20^2,"
					  "1-x08*x21^2+x09*x21^2,"
					  "x02+x01^3*x10*x16^3+x01^5*x11*x17^3+x01^7*x12*x18^3,"
					  "x03+x01^2*x10*x16^3+x01^5*x13*x19^3+x01^7*x14*x20^3,"
					  "x04+x01^2*x11*x17^3+x01^3*x13*x19^3+x01^7*x15*x21^3,"
					  "x05+x01^2*x12*x18^3+x01^3*x14*x20^3+x01^5*x15*x21^3,"
					  "x06+x01^3*x10^3*x16+x01^5*x11^3*x17+x01^7*x12^3*x18,"
					  "x07+x01^2*x10^3*x16+x01^5*x13^3*x19+x01^7*x14^3*x20,"
					  "x08+x01^2*x11^3*x17+x01^3*x13^3*x19+x01^7*x15^3*x21,"
					  "x09+x01^2*x12^3*x18+x01^3*x14^3*x20+x01^5*x15^3*x21,"
					  "x01^2*x02+x01^3*x03+x01^5*x04+x01^7*x05,"
					  "x01^2*x06+x01^3*x07+x01^5*x08+x01^7*x09,"
					  "x01^5*x11^2*x12^2*x13^2*x14^2*x15^2*x17^2*x18^2*x19^2*x20^2*x21^2-x10^2*x11^2*x12^2*x13^2*x14^2*x15^2*x16^2*x17^2*x18^2*x19^2*x20^2*x21^2*x22+x01^7*x10^2*x12^2*x13^2*x14^2*x15^2*x16^2*x18^2*x19^2*x20^2*x21"
					  "^2+x01^8*x10^2*x11^2*x12^2*x14^2*x15^2*x16^2*x17^2*x18^2*x20^2*x21^2+x01^9*x10^2*x11^2*x13^2*x14^2*x15^2*x16^2*x17^2*x19^2*x20^2*x21^2+x01^10*x10^2*x11^2*x12^2*x13^2*x15^2*x16^2*x17^2*x18^2*x19^2*x21^2+x01"
					  "^12*x10^2*x11^2*x12^2*x13^2*x14^2*x16^2*x17^2*x18^2*x19^2*x20^2,"
					  "-x22+x01^7*x10*x16+x01^9*x11*x17+x01^10*x13*x19+x01^11*x12*x18+x01^12*x14*x20+x01^14*x15*x21}"

			  "Q[z1,z2,z3,z4,w1,w2,w3,w4,a12,a13,a14,a23,a24,a34,b12,b13,b14,b23,b24,b34,t]{"
"a12^2-z2+z1,"
"b12^2-w2+w1,"
"a13^2-z3+z1,"
"b13^2-w3+w1,"
"a23^2-z3+z2,"
"b23^2-w3+w2,"
"a14^2-z4+z1,"
"b14^2-w4+w1,"
//"a24^2-z4+z2,"
//"b24^2-w4+w2,"
//"a34^2-z4+z3,"
//"b34^2-w4+w3,"
//"a12^(-1)*b12^(-3)*t^3+a13^(-1)*b13^(-3)*t^5+a14^(-1)*b14^(-3)*t^7-z1,"
//"a12^(-3)*b12^(-1)*t^3+a13^(-3)*b13^(-1)*t^5+a14^(-3)*b14^(-1)*t^7-w1,"
//"a12^(-1)*b12^(-3)*t^2+a23^(-1)*b23^(-3)*t^5+a24^(-1)*b24^(-3)*t^7-z2,"
//"a12^(-3)*b12^(-1)*t^2+a23^(-3)*b23^(-1)*t^5+a24^(-3)*b24^(-1)*t^7-w2,"
//"a13^(-1)*b13^(-3)*t^2+a23^(-1)*b23^(-3)*t^3+a34^(-1)*b34^(-3)*t^7-z3,"
//"a13^(-3)*b13^(-1)*t^2+a23^(-3)*b23^(-1)*t^3+a34^(-3)*b34^(-1)*t^7-w3,"
//"a14^(-1)*b14^(-3)*t^2+a24^(-1)*b24^(-3)*t^3+a34^(-1)*b34^(-3)*t^5-z4,"
"a14^(-3)*b14^(-1)*t^2+a24^(-3)*b24^(-1)*t^3+a34^(-3)*b34^(-1)*t^5-w4}"
			//  "Q[t,x,y]{x+y+t}"
			).parsePolynomialSetWithRing();
			  vector<gfan::Matrix<typ> > configurations;
			  for(auto &p:g)
				  configurations.push_back(convertMatrix<typ>(rowsToIntegerMatrix2(p.exponents())));
			  vector<HalfOpenCone<typ> > f;
			  //int n=3;
//			  vector<gfan::Matrix<typ> > configurations=convertMatrixVectorT(MixedVolumeExamples::cyclic(7));
//			  configurations.pop_back();configurations.pop_back();configurations.pop_back();//configurations.pop_back();
			  //We pretend that first coordinate is speciale (that it is the t parameter)
				int n=configurations[0].getWidth();
				int nonEmptyIntersections=0;
//				vector<int> edgeCounts;
//				for(int i=0;i<configurations.size();i++)edgeCounts.push_back(numberOfEdges(configurations[i]));
				vector<int> used(configurations.size());
				PolytopeIntersectionData<typ> data(configurations);
				gfan::Matrix<typ> strictInequalities(n,1);
				strictInequalities[0][0]=typ(-1);
				std::cerr<<"ETSEEETET\n";
				ProgressCounter progress;
				CommonStatistics statistics;
				commonRefinement(						gfan::HalfOpenCone<typ>(
						gfan::Matrix<typ>(n,0),
						gfan::Matrix<typ>(n,0),
						strictInequalities),used,configurations,f,nonEmptyIntersections/*,&edgeCounts*/,data,RelationTable(data.layout),progress,&statistics,8/*8*/);
				assert(!f.empty());
//				std::cerr<<"ABD"<<f.begin()->getDimension();//<<f.begin()->closure().getLinealitySpace();
				auto linealitySpace=f.begin()->closure().getLinealitySpace();
				RayCollector<typ> collector(f.begin()->closure().getLinealitySpace());

				if(0){
					vector<HalfOpenCone<typ> > F;
					for(auto &c:f)
						if(!c.representedPolyhedronIsBounded())
							if(!c.representedPolyhedronsRecessionConeIsContainedInNegativeHalfSpace())
								F.push_back(c);
					f=F;
				}
				int numberOfNonEmptyHalfOpenCones=0;
				for(auto &c:f)
				{
					if(!c.isEmpty())
					{
						numberOfNonEmptyHalfOpenCones++;
//						std::cerr<<"++++++++++++++++++++++++++++++++++++++\n";
						auto rays=c.closure().getRays();
						for(int i=0;i<rays.getHeight();i++)
							collector.lookup(normalize(rays[i].toVector()));
					}
				}
//				std::cerr<<"THE RAYS\n"<<collector.getRays()<<"\n";
				std::cerr<<"N_RAYS\n"<<collector.getRays().getHeight()<<"\n";
/*				for(auto &c:f)
				{
					std::cerr<<"HALFOPENCONE:\n";
					c.extractFaceComplex();
				}*/
///*				std::cerr<<*/ extractFaceComplex(f,collector,linealitySpace);

				std::cerr<<fvector(f).toString()<<"\n";

				std::cerr<<"Number of computed half open cones:"<<f.size()<<"\n";
				std::cerr<<"Number of non-empty half open cones:"<<numberOfNonEmptyHalfOpenCones<<"\n";
				std::cerr<<"NonEmptyIntersections:"<<statistics.numberOfIntermediateVertices<<"\n";
		  }

			static_assert(std::is_move_constructible<gfan::Matrix<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_nothrow_move_constructible<gfan::Matrix<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_move_constructible<gfan::HalfOpenCone<gfan::CircuitTableInt32>>::value,"");
//			static_assert(std::is_trivially_move_constructible<gfan::HalfOpenCone<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_nothrow_move_constructible<gfan::Tableau<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_nothrow_move_constructible<gfan::GeneratedCone<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_nothrow_move_constructible<gfan::Cone<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_move_assignable<gfan::Cone<gfan::CircuitTableInt32>>::value,"");
			static_assert(std::is_nothrow_move_constructible<gfan::HalfOpenCone<gfan::CircuitTableInt32>>::value,"");

		  return 0;
	}
	{
//		typedef gfan::CircuitTableInt64 typ;
//		typedef gfan::CircuitTableInteger typ;
		typedef gfan::CircuitTableInt128 typ;
		std::cerr<<"Digits:"<<std::numeric_limits<char>::digits<<"\n";
		std::cerr<<"Digits:"<<std::numeric_limits<short>::digits<<"\n";
		std::cerr<<"Digits:"<<std::numeric_limits<int>::digits<<"\n";
		std::cerr<<"Digits:"<<std::numeric_limits<long int>::digits<<"\n";
		std::cerr<<"Digits:"<<std::numeric_limits<__int128_t>::digits<<"\n";


		{
			  __int128 t=0;
			  __int128 A=2;
			  __int128 s=-4;
			  __int128 B=1;
			  B=(B<<127)+3;
			  std::cerr<<"---B\n";
			  std::cerr<<toStr(t)<<"\n";
			  std::cerr<<toStr(A)<<"\n";
			  std::cerr<<toStr(extMul(t,A))<<"\n";
			  std::cerr<<toStr(s)<<"\n";
			  std::cerr<<toStr(B)<<"\n";
			  std::cerr<<toStr(extMul(s,B))<<"\n";
			  std::cerr<<"---B\n";
			  std::cerr<<toStr(unsignedProd128(s,B));

		  }


		FileParser F(stdin);
		std::cerr<<"___\n";
		auto L=F.parseIntegerVectorList();
		std::cerr<<"___\n";
		auto rays=convertMatrix<typ>(rowsToIntegerMatrix2(L)).transposed();
		auto lines=gfan::Matrix<typ>(rays.getHeight(),0);
		GeneratedCone<typ> C(rays,lines);
		std::cerr<<"Dimension of Lineality Space:"<<C.getDimensionOfLinealitySpace()<<"\n";

		return 0;
	}
		{
			if(MixedVolumeExamples::cyclic(3).size()!=3)std::cerr<<"This mistake could be caused by linking libs without/without bound checking.\n";
		gfan::IntMatrix ex=MixedVolumeExamples::cyclic(3)[1];
		std::cerr<<"1c\n";
		gfan::Matrix<CircuitTableInt32>  A=convertMatrix<CircuitTableInt32>(ex);
/*		gfan::Matrix<CircuitTableInt32>  A(3,2);
		A[1][0]=1;
		A[2][1]=1;*/
		/*=
				restrictedTropicalHypersurface(
						gfan::HalfOpenCone<CircuitTableInt32>(
								gfan::Matrix<CircuitTableInt32>(0,3),
								gfan::Matrix<CircuitTableInt32>(0,3)),
								A
				);*/

		initializeCddlibIfRequired();//SEEMS TO BE REQUIRED HERE. WHERE DOES THE INITIALISATION USUALLY HAPPEN?
		int largest=10+2;
		vector<int> halfOpenConesStatistics(largest+1);
		vector<int> nonEmptyIntersectionsStatistics(largest+1);
		for(int k=2/*largest/*2*/;k<=largest;k++)
		{
			vector<HalfOpenCone<CircuitTableInt32> > f;
			vector<gfan::Matrix<CircuitTableInt32> > configurations=convertMatrixVectorT(MixedVolumeExamples::cyclic(k));
if(0)			{
				auto temp=configurations;
				temp[0]=configurations[7];
				temp[1]=configurations[3];
				temp[2]=configurations[1];
				temp[3]=configurations[5];
				temp[4]=configurations[6];
				temp[5]=configurations[0];
				temp[6]=configurations[4];
				temp[7]=configurations[2];
				configurations=temp;
			}
if(0)			{
				auto temp=configurations;
				temp[0]=configurations[9];
				temp[1]=configurations[4];
				temp[2]=configurations[1];
				temp[3]=configurations[7];
				temp[4]=configurations[3];
				temp[5]=configurations[5];
				temp[6]=configurations[0];
				temp[7]=configurations[8];
				temp[8]=configurations[2];
				temp[9]=configurations[6];
				configurations=temp;
			}
//	while(configurations.size()>3)configurations.pop_back();
			int n=configurations[0].getWidth();
			std::cerr<<"n:"<<n<<"\n";
			int nonEmptyIntersections=0;
//			vector<int> edgeCounts;
//			for(int i=0;i<configurations.size();i++)edgeCounts.push_back(numberOfEdges(configurations[i]));
			vector<int> used(configurations.size());
			PolytopeIntersectionData<CircuitTableInt32> data(configurations);
			//std::cerr<<data.toString();
			ProgressCounter progress;
			CommonStatistics statistics;
			commonRefinement(						gfan::HalfOpenCone<CircuitTableInt32>(
					gfan::Matrix<CircuitTableInt32>(n,0),
					gfan::Matrix<CircuitTableInt32>(n,0),
					gfan::Matrix<CircuitTableInt32>(n,0)),used,configurations,f,nonEmptyIntersections,/*&edgeCounts,*/data,RelationTable(data.layout),progress,&statistics,8);
			if(1)
			{
				Vector<int> sum(n+1);
				for(auto &a:f)sum=sum+a.fVector();
				std::cerr<<"F_VECTOR:"<<sum<<"\n";
			}

/*			{
				for(int i=0;i<f.size();i++)
					for(int j=0;j<i;j++)
					{
						assert(intersection2(f[i],f[j]).isEmpty());
					}
			}
*/

//			std::cerr<<"RESULT:"<<gfan::toString(f,n);

//			std::cerr<<"Refinementsize:"<<f.size()<<"\n";
			halfOpenConesStatistics[k]=f.size();
			nonEmptyIntersectionsStatistics[k]=statistics.numberOfIntermediateVertices;//nonEmptyIntersections;


			if(1)
			{
				assert(!f.begin()->isEmpty());
				auto linealitySpace=f.begin()->closure().getLinealitySpace();

				RayCollector<CircuitTableInt32> collector(linealitySpace);
				for(auto &c:f)
				{
					if(!c.isEmpty())
					{
						auto rays=c.closure().getRays();
						for(int i=0;i<rays.getHeight();i++)
							collector.lookup(normalize(rays[i].toVector()));
					}
				}
				std::cerr<<"N_RAYS\n"<<collector.getRays().getHeight()<<"\n";

				//std::cerr<<extractFaceComplex(f,collector,linealitySpace,0).toString(FPF_cones|FPF_maximalCones);
			}

			//			print(f);
		}
		std::cerr<<"Example\tResult\tIntermediate\n";
		for(int i=2;i<=largest;i++)
			std::cerr<<"Cyclic-"<<i<<"\t:"<<halfOpenConesStatistics[i]<<"\t:"<<nonEmptyIntersectionsStatistics[i]<<"\n";

		return 0;

		for(int i=0;i<10;i++)
		{
			gfan::Matrix<CircuitTableInt32> A=randomMatrix2(4,3);
			Cone<CircuitTableInt32> C(A.transposed());
			std::cerr<<"------------------------------------\n"<<matrixToString(A)<<"(rows give inequalities)\n";
			std::cerr<<C.toString();
		}
		return 0;

//			gfan::TableauSolver<CircuitTableInt32> M2(randomMatrix(4,15),true);
		for(int i=0;i<100;i++)
		{
			gfan::Matrix<CircuitTableInt32> A=randomMatrix(3,3);
			debug<<A.toString();
			gfan::GeneratedCone<CircuitTableInt32> M2(A);
			debug<<M2.toString();
			debug<<M2.getRays().toString();
			debug<<"orthogonal complement:"<<M2.getOrthogonalComplement().toString();
			debug<<"rays:"<<M2.getRays().toString();
		}
			return 0;
		}
		//		CircuitTableInt32
		{
			gfan::Matrix<CircuitTableInt32> M(4,8);

			M[0][0]=1;M[0][1]=-1;M[0][2]=1;M[0][3]=-1;M[0][4]=1;M[0][5]=-1;M[0][6]=1;M[0][7]=-1;
			M[1][0]=1;M[1][1]=1;M[1][2]=-1;M[1][3]=-1;M[1][4]=1;M[1][5]=1;M[1][6]=-1;M[1][7]=-1;
			M[2][0]=1;M[2][1]=1;M[2][2]=1;M[2][3]=1;M[2][4]=-1;M[2][5]=-1;M[2][6]=-1;M[2][7]=-1;
			M[3][0]=1;M[3][1]=1;M[3][2]=1;M[3][3]=1;M[3][4]=1;M[3][5]=1;M[3][6]=1;M[3][7]=1;
			gfan::GeneratedCone<CircuitTableInt32> M2(M);

			debug<<M2.getRays().toString();
//			M2.coneInfo();
			return 0;
		}

		{
			gfan::Matrix<CircuitTableInt32> M(3,5);

			M[0][0]=0;M[0][1]=0;M[0][2]=2;M[0][3]=2;M[0][4]=1;
			M[1][0]=1;M[1][1]=0;M[1][2]=3;M[1][3]=2;M[1][4]=1;
			M[2][0]=2;M[2][1]=0;M[2][2]=2;M[2][3]=0;M[2][4]=0;
			gfan::GeneratedCone<CircuitTableInt32> M2(M);

			debug<<M2.getRays().toString();
//			M2.coneInfo();
			return 0;
		}


		{
			gfan::Matrix<CircuitTableInt32> M(2,2);
			M[0][0]=2;
			M[1][0]=1;
			M[1][1]=2;
			M[0][1]=2;
			gfan::Tableau<CircuitTableInt32> M2(M,true);
			debug<<M2.toString();
			M2.exchange(1,2);
			debug<<M2.toString();
		}
		gfan::Matrix<CircuitTableInt32> M(3,6);
		M[0][0]=1;
		M[1][0]=0;
		M[2][0]=1;
		M[0][1]=0;
		M[1][1]=1;
		M[2][1]=1;
		M[0][2]=2;
		M[1][2]=2;
		M[2][2]=1;
		M[0][3]=5;
		M[1][3]=1;
		M[2][3]=1;
		M[0][4]=4;
		M[1][4]=0;
		M[2][4]=1;
		M[0][5]=3;
		M[1][5]=3;
		M[2][5]=3;
		gfan::GeneratedCone<CircuitTableInt32> M2(M);
		debug<<M2.toString();
		M2.exchange(0,3);
		debug<<M2.toString();
		M2.exchange(1,4);
		debug<<M2.toString();
		M2.exchange(2,5);
		debug<<M2.toString();

//		M2.solve();

		return 0;

		debug<<"ETST3213312\n";
		ZMatrix A=ZMatrix(3,3);
		A[0][0]=1;
		A[1][0]=1;

		A[1][1]=1;
		ZMatrix B=ZMatrix(0,3);
		for(int i=0;i<4000;i++)
		{		ZCone C(A,B);}
		return 0;
	}

	int testLLL()
	{
		IntegerMatrix A=rowsToIntegerMatrix(StringParser(
				"{"
/*				"(1,0,1,1, 0,0, 0, 0,0),"
				"(1,1,0, 0,1,0, 0, 0,0),"
				"(1,0,0, 0, 0,1,1, 0,0),"
				"(0,1,0, 0, 0,1, 0,1,0),"
				"(0,1,0, 0, 0,0, 0, 0,1),"
				"(0,0,0, 0, 0,1, 0, 0,1)}"
*/
				"( 1, 0,1,1, 0, 0, 0, 0,0),"
				"( 1, 1,0,0, 1, 0, 0, 0,0),"
				"( 1, 0,0,0, 0, 1, 1, 0,0),"
				"(-1, 0,1,1,-1, 0,-1, 1,0),"
				"( 0, 0,0,0, 0,-1, 0,-1,1),"
				"( 1,-1,0,0, 1, 0, 0,-1,0)}"
				).parseIntegerVectorList());
		A=A.submatrixColumnSubsetBoolean(IntegerVector(StringParser("( 1,1,1, 1, 1,1, 1, 1,1)").parseIntegerVector()));

		{
			FieldMatrix B=integerMatrixToFieldMatrix(A,Q);
			debug <<B.getHeight()<<"x"<<B.getWidth()<<"rank"<<B.reduceAndComputeRank();
		}

		debug<<"Doing LLL....\n"<<A;
		mlll(A);
		debug<<"LLL done\n";
		return 0;
	}

	int testPolynomialGCD()
	{
		PolynomialSet A=FileParser(stdin).parsePolynomialSetWithRing();
		Polynomial g=NonMonomialPolynomialGCDForZ(A);
		pout<<g<<"\n";
//		debug<<"nterms:"<<(int)g.terms.size()<<"\n"<<g<<"\n";
		return 0;
		{
			FILE *f=fopen("gcdexamples","r");
			FileParser F(f);
//			while(!feof(f))
			for(int i=0;i<20;i++)
			{
				PolynomialSet A=F.parsePolynomialSetWithRing();
				int nterms=F.parseInt();
//				Polynomial g=NonMonomialPolynomialGCDForZModP(A);
				Polynomial g=NonMonomialPolynomialGCDForZ(A);
				debug<<"nterms:"<<(int)g.terms.size()<<"\n"<<g<<"\n";
			}
			assert(0);
		}
//		PolynomialRing r=StringParser("Z/9931Z[a,b]").parsePolynomialRing();
		PolynomialRing r=StringParser("Z/9851Z[a,b]").parsePolynomialRing();
#if 1
		if(1)
		{
			PolynomialSet A=StringParser(
//					"Z/9851Z[x0,x1,x2,x3,x4,x5]"
#if 0
					"Z/31081Z[x1,x3,x4]{"
					"10480*x3^4+23521*x3^4*x4+17445*x3^4*x4^2+14816*x3^4*x4^3+21952*x3^4*x4^4+15402*x3^5+30564*x3^5*x4+7415*x3^5*x4^2+16201*x3^5*x4^3+27099*x3^5*x4^4+19917*x3^5*x4^5+19612*x3^6+15302*x3^6*x4+3612*x3^6*x4^2+12310*x3^6*x4^3+16477*x3^6*x4^4+7051*x3^6*x4^5+20564*x3^7+19791*x3^7*x4+21991*x3^7*x4^2+12637*x3^7*x4^3+6057*x3^7*x4^4+23231*x3^7*x4^5+14535*x3^8+26782*x3^8*x4+20681*x3^8*x4^2+6649*x3^8*x4^3+23550*x3^8*x4^4+14190*x3^9+1471*x3^9*x4+2489*x3^9*x4^2+7531*x3^9*x4^3+22060*x3^10+8835*x3^10*x4+7850*x3^10*x4^2+30722*x1*x3^3+19791*x1*x3^3*x4+1353*x1*x3^3*x4^2+24071*x1*x3^3*x4^3+29424*x1*x3^3*x4^4+21613*x1*x3^3*x4^5+8339*x1*x3^4+7465*x1*x3^4*x4+13651*x1*x3^4*x4^2+5486*x1*x3^4*x4^3"
					"+17869*x1*x3^4*x4^4+4951*x1*x3^4*x4^5+23731*x1*x3^4*x4^6+15655*x1*x3^5+28546*x1*x3^5*x4+11909*x1*x3^5*x4^2+5517*x1*x3^5*x4^3+27556*x1*x3^5*x4^4+26370*x1*x3^5*x4^5+10806*x1*x3^5*x4^6+2870*x1*x3^6+9972*x1*x3^6*x4+18303*x1*x3^6*x4^2+5649*x1*x3^6*x4^3+23575*x1*x3^6*x4^4+29235*x1*x3^6*x4^5-x1*x3^6*x4^6+11084*x1*x3^7+26810*x1*x3^7*x4+13588*x1*x3^7*x4^2+21113*x1*x3^7*x4^3+17736*x1*x3^7*x4^4+3*x1*x3^7*x4^5+23398*x1*x3^8+14883*x1*x3^8*x4+13836*x1*x3^8*x4^2+20085*x1*x3^8*x4^3+31078*x1*x3^8*x4^4+5003*x1*x3^9+15699*x1*x3^9*x4+15381*x1*x3^9*x4^2+x1*x3^9*x4^3+359*x1^2*x3^2+17068*x1^2*x3^2*x4+15595*x1^2*x3^2*x4^2+25377*x1^2*x3^2*x4^3+28664*x1^2*x3^2*x4^4+28404*x1^2*x3^2*x4^5+6895*x1^2*x3^3+10231*x1^2*x3^3*x4+19516*x1^2*x3^3*x4^2+19906*x1^2*x3^3*x4^3+8900*x1^2*x3^3*x4^4+14974*x1^2*x3^3*x4^5+29400*x1^2*x3^3*x4^6+21849*x1^2*x3^4+8606*x1^2*x3^4*x4+22936*x1^2*x3^4*x4^2+22769*x1^2*x3^4*x4^3+15774*x1^2*x3^4*x4^4+15899*x1^2*x3^4*x4^5+18938*x1^2*x3^4*x4^6+23411*x1^2*x3^5"
					"+14860*x1^2*x3^5*x4+22744*x1^2*x3^5*x4^2+21410*x1^2*x3^5*x4^3+15407*x1^2*x3^5*x4^4+23722*x1^2*x3^5*x4^5+4*x1^2*x3^5*x4^6+3862*x1^2*x3^6+24483*x1^2*x3^6*x4+14308*x1^2*x3^6*x4^2+20672*x1^2*x3^6*x4^3+4366*x1^2*x3^6*x4^4+31069*x1^2*x3^6*x4^5+71*x1^2*x3^7+11950*x1^2*x3^7*x4+7606*x1^2*x3^7*x4^2+30836*x1^2*x3^7*x4^3+12*x1^2*x3^7*x4^4+8036*x1^2*x3^8+4259*x1^2*x3^8*x4+15381*x1^2*x3^8*x4^2+31077*x1^2*x3^8*x4^3+20601*x1^3*x3+29974*x1^3*x3*x4+19295*x1^3*x3*x4^2+4255*x1^3*x3*x4^3+15240*x1^3*x3*x4^4+2677*x1^3*x3*x4^5+23909*x1^3*x3^2+22930*x1^3*x3^2*x4+2963*x1^3*x3^2*x4^2+25085*x1^3*x3^2*x4^3+30230*x1^3*x3^2*x4^4+5931*x1^3*x3^2*x4^5+18062*x1^3*x3^2*x4^6+17375*x1^3*x3^3+22858*x1^3*x3^3*x4+10104*x1^3*x3^3*x4^2+24139*x1^3*x3^3*x4^3+26625*x1^3*x3^3*x4^4+3302*x1^3*x3^3*x4^5+2674*x1^3*x3^3*x4^6+17717*x1^3*x3^4+8488*x1^3*x3^4*x4+6093*x1^3*x3^4*x4^2+24154*x1^3*x3^4*x4^3+407*x1^3*x3^4*x4^4+18410*x1^3*x3^4*x4^5+31075*x1^3*x3^4*x4^6+25208*x1^3*x3^5+26248*x1^3*x3^5*x4"
					"+10629*x1^3*x3^5*x4^2+24545*x1^3*x3^5*x4^3+17958*x1^3*x3^5*x4^4+18*x1^3*x3^5*x4^5+28309*x1^3*x3^6+15547*x1^3*x3^6*x4+24681*x1^3*x3^6*x4^2+22482*x1^3*x3^6*x4^3+31063*x1^3*x3^6*x4^4+5003*x1^3*x3^7+22246*x1^3*x3^7*x4+638*x1^3*x3^7*x4^2+6*x1^3*x3^7*x4^3+2889*x1^4*x4+8474*x1^4*x4^2+24724*x1^4*x4^3+29044*x1^4*x4^4+9468*x1^4*x4^5+7617*x1^4*x3+20496*x1^4*x3*x4+30333*x1^4*x3*x4^2+4730*x1^4*x3*x4^3+16679*x1^4*x3*x4^4+10313*x1^4*x3*x4^5+29400*x1^4*x3*x4^6+18752*x1^4*x3^2+17501*x1^4*x3^2*x4+693*x1^4*x3^2*x4^2+17364*x1^4*x3^2*x4^3+17629*x1^4*x3^2*x4^4+22531*x1^4*x3^2*x4^5+18938*x1^4*x3^2*x4^6+28681*x1^4*x3^3+17168*x1^4*x3^3*x4+16144*x1^4*x3^3*x4^2+19248*x1^4*x3^3*x4^3+14866*x1^4*x3^3*x4^4+810*x1^4*x3^3*x4^5+4*x1^4*x3^3*x4^6+7473*x1^4*x3^4+5877*x1^4*x3^4*x4+26709*x1^4*x3^4*x4^2+19845*x1^4*x3^4*x4^3+10940*x1^4*x3^4*x4^4+31069*x1^4*x3^4*x4^5+27275*x1^4*x3^5+21574*x1^4*x3^5*x4"
					"+15436*x1^4*x3^5*x4^2+24262*x1^4*x3^5*x4^3+12*x1^4*x3^5*x4^4+22060*x1^4*x3^6+22246*x1^4*x3^6*x4+7212*x1^4*x3^6*x4^2+31077*x1^4*x3^6*x4^3+1557*x1^5*x4+19365*x1^5*x4^2+21835*x1^5*x4^3+23547*x1^5*x4^4+6076*x1^5*x4^5+23731*x1^5*x4^6+430*x1^5*x3*x4+12908*x1^5*x3*x4^2+11144*x1^5*x3*x4^3+20263*x1^5*x3*x4^4+18090*x1^5*x3*x4^5+10806*x1^5*x3*x4^6+22964*x1^5*x3^2*x4+7968*x1^5*x3^2*x4^2+10145*x1^5*x3^2*x4^3+1850*x1^5*x3^2*x4^4+28916*x1^5*x3^2*x4^5-x1^5*x3^2*x4^6+14124*x1^5*x3^3*x4+7328*x1^5*x3^3*x4^2+419*x1^5*x3^3*x4^3+18693*x1^5*x3^3*x4^4+3*x1^5*x3^3*x4^5+27818*x1^5*x3^4*x4+29195*x1^5*x3^4*x4^2+19128*x1^5*x3^4*x4^3+31078*x1^5*x3^4*x4^4+19958*x1^5*x3^5*x4+15700*x1^5*x3^5*x4^2+x1^5*x3^5*x4^3,"
					"2889*x3^5+13894*x3^5*x4+26401*x3^5*x4^2+9129*x3^5*x4^3+27088*x3^6+24707*x3^6*x4+15089*x3^6*x4^2+15002*x3^6*x4^3+11164*x3^6*x4^4+15051*x3^7+3809*x3^7*x4+15009*x3^7*x4^2+9693*x3^7*x4^3+24030*x3^7*x4^4+10070*x3^8+30804*x3^8*x4+2943*x3^8*x4^2+5066*x3^8*x4^3+7850*x3^8*x4^4+23299*x3^9+16033*x3^9*x4+22144*x3^9*x4^2+7531*x3^9*x4^3+22329*x3^10+30880*x3^10*x4+23550*x3^10*x4^2+11123*x3^11+23231*x3^11*x4+19525*x1*x3^4+1166*x1*x3^4*x4+17043*x1*x3^4*x4^2+14480*x1*x3^4*x4^3+9468*x1*x3^4*x4^4+17529*x1*x3^5+28033*x1*x3^5*x4+16621*x1*x3^5*x4^2+1290*x1*x3^5*x4^3+3802*x1*x3^5*x4^4+7350*x1*x3^5*x4^5+2388*x1*x3^6+7014*x1*x3^6*x4+23679*x1*x3^6*x4^2+3783*x1*x3^6*x4^3+18813*x1*x3^6*x4^4+20275*x1*x3^6*x4^5+13765*x1*x3^7+11545*x1*x3^7*x4+19388*x1*x3^7*x4^2+15044*x1*x3^7*x4^3+17227*x1*x3^7*x4^4+x1*x3^7*x4^5+14171*x1*x3^8+25057*x1*x3^8*x4+5913*x1*x3^8*x4^2+29364*x1*x3^8*x4^3+31078*x1*x3^8*x4^4+664*x1*x3^9+8495*x1*x3^9*x4+26058*x1*x3^9*x4^2+3*x1*x3^9*x4^3+6547*x1*x3^10+319*x1*x3^10*x4-x1*x3^10*x4^2+17334*x1^2*x3^3+11801*x1^2*x3^3*x4+9709*x1^2*x3^3*x4^2+6195*x1^2*x3^3*x4^3"
					"+24290*x1^2*x3^3*x4^4+895*x1^2*x3^4+8041*x1^2*x3^4*x4+26035*x1^2*x3^4*x4^2+16446*x1^2*x3^4*x4^3+14716*x1^2*x3^4*x4^4+25412*x1^2*x3^4*x4^5+26424*x1^2*x3^5+3816*x1^2*x3^5*x4+6147*x1^2*x3^5*x4^2+30203*x1^2*x3^5*x4^3+6297*x1^2*x3^5*x4^4+22949*x1^2*x3^5*x4^5+30726*x1^2*x3^6+16799*x1^2*x3^6*x4+4313*x1^2*x3^6*x4^2+20404*x1^2*x3^6*x4^3+13682*x1^2*x3^6*x4^4+31076*x1^2*x3^6*x4^5+21136*x1^2*x3^7+6487*x1^2*x3^7*x4+26119*x1^2*x3^7*x4^2+19944*x1^2*x3^7*x4^3+15*x1^2*x3^7*x4^4+22702*x1^2*x3^8+28105*x1^2*x3^8*x4+13756*x1^2*x3^8*x4^2+31066*x1^2*x3^8*x4^3+17987*x1^2*x3^9+22912*x1^2*x3^9*x4+5*x1^2*x3^9*x4^2+19525*x1^3*x3^2+5147*x1^3*x3^2*x4+8658*x1^3*x3^2*x4^2+20812*x1^3*x3^2*x4^3+25727*x1^3*x3^2*x4^4+25314*x1^3*x3^3+1472*x1^3*x3^3*x4+18409*x1^3*x3^3*x4^2+25329*x1^3*x3^3*x4^3+21727*x1^3*x3^3*x4^4+11338*x1^3*x3^3*x4^5+4538*x1^3*x3^4+24902*x1^3*x3^4*x4+28981*x1^3*x3^4*x4^2+20470*x1^3*x3^4*x4^3+20945*x1^3*x3^4*x4^4+16264*x1^3*x3^4*x4^5+4261*x1^3*x3^5+9380*x1^3*x3^5*x4+9737*x1^3*x3^5*x4^2+8515*x1^3*x3^5*x4^3+20055*x1^3*x3^5*x4^4+10*x1^3*x3^5*x4^5+22629*x1^3*x3^6+24442*x1^3*x3^6*x4+15134*x1^3*x3^6*x4^2+4341*x1^3*x3^6*x4^3"
					"+31051*x1^3*x3^6*x4^4+15430*x1^3*x3^7+24217*x1^3*x3^7*x4+21502*x1^3*x3^7*x4^2+30*x1^3*x3^7*x4^3+13094*x1^3*x3^8+31071*x1^3*x3^8*x4^2+2889*x1^4*x3+4493*x1^4*x3*x4+2028*x1^4*x3*x4^2+22712*x1^4*x3*x4^3+24290*x1^4*x3*x4^4+20860*x1^4*x3^2+11266*x1^4*x3^2*x4+17638*x1^4*x3^2*x4^2+22541*x1^4*x3^2*x4^3+14015*x1^4*x3^2*x4^4+19743*x1^4*x3^2*x4^5+13331*x1^4*x3^3+23771*x1^4*x3^3*x4+4789*x1^4*x3^3*x4^2+9281*x1^4*x3^3*x4^3+3504*x1^4*x3^3*x4^4+14817*x1^4*x3^3*x4^5+11457*x1^4*x3^4+8419*x1^4*x3^4*x4+17696*x1^4*x3^4*x4^2+29654*x1^4*x3^4*x4^3+2857*x1^4*x3^4*x4^4+31071*x1^4*x3^4*x4^5+28965*x1^4*x3^5+18307*x1^4*x3^5*x4+28214*x1^4*x3^5*x4^2+20166*x1^4*x3^5*x4^3+30*x1^4*x3^5*x4^4+4300*x1^4*x3^6+18675*x1^4*x3^6*x4+16153*x1^4*x3^6*x4^2+31051*x1^4*x3^6*x4^3+24534*x1^4*x3^7+8169*x1^4*x3^7*x4+10*x1^4*x3^7*x4^2+25661*x1^5*x4+29404*x1^5*x4^2+19915*x1^5*x4^3+9468*x1^5*x4^4+1557*x1^5*x3+25453*x1^5*x3*x4+4690*x1^5*x3*x4^2+27495*x1^5*x3*x4^3+10579*x1^5*x3*x4^4+5669*x1^5*x3*x4^5+430*x1^5*x3^2+22131*x1^5*x3^2*x4+12333*x1^5*x3^2*x4^2+1872*x1^5*x3^2*x4^3"
					"+8615*x1^5*x3^2*x4^4+8132*x1^5*x3^2*x4^5+22964*x1^5*x3^3+18249*x1^5*x3^3*x4+21114*x1^5*x3^3*x4^2+7644*x1^5*x3^3*x4^3+25887*x1^5*x3^3*x4^4+5*x1^5*x3^3*x4^5+14124*x1^5*x3^4+15601*x1^5*x3^4*x4+4237*x1^5*x3^4*x4^2+16754*x1^5*x3^4*x4^3+31066*x1^5*x3^4*x4^4+27818*x1^5*x3^5+16039*x1^5*x3^5*x4+11708*x1^5*x3^5*x4^2+15*x1^5*x3^5*x4^3+19958*x1^5*x3^6+30762*x1^5*x3^6*x4+31076*x1^5*x3^6*x4^2+25352*x1^6*x4+25842*x1^6*x4^2+16221*x1^6*x4^3+17240*x1^6*x4^4+23731*x1^6*x4^5+7800*x1^6*x3*x4+2305*x1^6*x3*x4^2+17941*x1^6*x3*x4^3+11039*x1^6*x3*x4^4+10806*x1^6*x3*x4^5+29128*x1^6*x3^2*x4+18052*x1^6*x3^2*x4^2+6916*x1^6*x3^2*x4^3+5685*x1^6*x3^2*x4^4-x1^6*x3^2*x4^5+18397*x1^6*x3^3*x4+22563*x1^6*x3^3*x4^2+26224*x1^6*x3^3*x4^3+3*x1^6*x3^3*x4^4+28994*x1^6*x3^4*x4+11597*x1^6*x3^4*x4^2+31078*x1^6*x3^4*x4^3+7850*x1^6*x3^5*x4+x1^6*x3^5*x4^2}"
#elif 1
					"Z/31081Z[x0,x1,x2,x3,x4]"
					"{x0*x1^3*x2*x4^2+x2^2*x3^3*x4^2-x2^2*x3^4*x4-x2^3*x3^3*x4+x2^3*x3^4-2*x1*x2*x3^2*x4^3+2*x1*x2*x3^3*x4^2+x1*x2^2*x3^2*x4^2-x1*x2^2*x3^3*x4+x1*x2^3*x3^2*x4-x1*x2^3*x3^3+x1^2*x3*x4^4-x1^2*x3^2*x4^3+x1^2*x2*x3*x4^3-x1^2*x2*x3^2*x4^2-2*x1^2*x2^2*x3*x4^2+2*x1^2*x2^2*x3^2*x4-x1^3*x4^4+x1^3*x3*x4^3+x1^3*x2*x4^3-x1^3*x2*x3*x4^2+x0*x2^2*x3^3*x4-x0*x2^3*x3^3-2*x0*x1*x2*x3^2*x4^2+x0*x1*x2^2*x3^2*x4+x0*x1*x2^3*x3^2+x0*x1^2*x3*x4^3+x0*x1^2*x2*x3*x4^2-2*x0*x1^2*x2^2*x3*x4-x0*x1^3*x4^3,"
					"x0*x1^6*x4^4+x2^4*x3^6*x4-x2^4*x3^7-4*x1*x2^3*x3^5*x4^2+4*x1*x2^3*x3^6*x4-2*x1*x2^4*x3^5*x4+2*x1*x2^4*x3^6+6*x1^2*x2^2*x3^4*x4^3-6*x1^2*x2^2*x3^5*x4^2+8*x1^2*x2^3*x3^4*x4^2-8*x1^2*x2^3*x3^5*x4+x1^2*x2^4*x3^4*x4-x1^2*x2^4*x3^5-4*x1^3*x2*x3^3*x4^4+4*x1^3*x2*x3^4*x4^3-12*x1^3*x2^2*x3^3*x4^3+12*x1^3*x2^2*x3^4*x4^2-4*x1^3*x2^3*x3^3*x4^2+4*x1^3*x2^3*x3^4*x4+x1^4*x3^2*x4^5-x1^4*x3^3*x4^4+8*x1^4*x2*x3^2*x4^4-8*x1^4*x2*x3^3*x4^3+6*x1^4*x2^2*x3^2*x4^3-6*x1^4*x2^2*x3^3*x4^2-2*x1^5*x3*x4^5+2*x1^5*x3^2*x4^4-4*x1^5*x2*x3*x4^4+4*x1^5*x2*x3^2*x4^3+x1^6*x4^5-x1^6*x3*x4^4+x0*x2^4*x3^6-4*x0*x1*x2^3*x3^5*x4-2*x0*x1*x2^4*x3^5+6*x0*x1^2*x2^2*x3^4*x4^2+8*x0*x1^2*x2^3*x3^4*x4+x0*x1^2*x2^4*x3^4-4*x0*x1^3*x2*x3^3*x4^3-12*x0*x1^3*x2^2*x3^3*x4^2-4*x0*x1^3*x2^3*x3^3*x4+x0*x1^4*x3^2*x4^4+8*x0*x1^4*x2*x3^2*x4^3+6*x0*x1^4*x2^2*x3^2*x4^2-2*x0*x1^5*x3*x4^4-4*x0*x1^5*x2*x3*x4^3}"
//					"-x2^2*x3^3*x4+x2^2*x3^4+2*x1*x2*x3^2*x4^2-2*x1*x2*x3^3*x4+x1*x2^2*x3^2*x4-x1*x2^2*x3^3-x1^2*x3*x4^3+x1^2*x3^2*x4^2-2*x1^2*x2*x3*x4^2+2*x1^2*x2*x3^2*x4+x1^3*x4^3-x1^3*x3*x4^2-x0*x2^2*x3^3+2*x0*x1*x2*x3^2*x4+x0*x1*x2^2*x3^2-x0*x1^2*x3*x4^2-2*x0*x1^2*x2*x3*x4+x0*x1^3*x4^2"

#elif 0
					"Z/31081Z[x0,x1,x2,x3,x4]"
					"{x2*x3^7*x4^5-3*x2*x3^8*x4^4+3*x2*x3^9*x4^3-x2*x3^10*x4^2-2*x2^2*x3^7*x4^4+6*x2^2*x3^8*x4^3-6*x2^2*x3^9*x4^2+2*x2^2*x3^10*x4+x2^3*x3^7*x4^3-3*x2^3*x3^8*x4^2+3*x2^3*x3^9*x4-x2^3*x3^10-x1*x3^6*x4^6+3*x1*x3^7*x4^5-3*x1*x3^8*x4^4+x1*x3^9*x4^3-2*x1*x2*x3^6*x4^5+6*x1*x2*x3^7*x4^4-6*x1*x2*x3^8*x4^3+2*x1*x2*x3^9*x4^2+7*x1*x2^2*x3^6*x4^4-21*x1*x2^2*x3^7*x4^3+21*x1*x2^2*x3^8*x4^2-7*x1*x2^2*x3^9*x4-4*x1*x2^3*x3^6*x4^3+12*x1*x2^3*x3^7*x4^2-12*x1*x2^3*x3^8*x4+4*x1*x2^3*x3^9+4*x1^2*x3^5*x4^6-12*x1^2*x3^6*x4^5+12*x1^2*x3^7*x4^4-4*x1^2*x3^8*x4^3-2*x1^2*x2*x3^5*x4^5+6*x1^2*x2*x3^6*x4^4-6*x1^2*x2*x3^7*x4^3+2*x1^2*x2*x3^8*x4^2-8*x1^2*x2^2*x3^5*x4^4+24*x1^2*x2^2*x3^6*x4^3-24*x1^2*x2^2*x3^7*x4^2+8*x1^2*x2^2*x3^8*x4+6*x1^2*x2^3*x3^5*x4^3-18*x1^2*x2^3*x3^6*x4^2+18*x1^2*x2^3*x3^7*x4-6*x1^2*x2^3*x3^8-6*x1^3*x3^4*x4^6+18*x1^3*x3^5*x4^5-18*x1^3*x3^6*x4^4+6*x1^3*x3^7*x4^3+8*x1^3*x2*x3^4*x4^5-24*x1^3*x2*x3^5*x4^4"
					"+24*x1^3*x2*x3^6*x4^3-8*x1^3*x2*x3^7*x4^2+2*x1^3*x2^2*x3^4*x4^4-6*x1^3*x2^2*x3^5*x4^3+6*x1^3*x2^2*x3^6*x4^2-2*x1^3*x2^2*x3^7*x4-4*x1^3*x2^3*x3^4*x4^3+12*x1^3*x2^3*x3^5*x4^2-12*x1^3*x2^3*x3^6*x4+4*x1^3*x2^3*x3^7+4*x1^4*x3^3*x4^6-12*x1^4*x3^4*x4^5+12*x1^4*x3^5*x4^4-4*x1^4*x3^6*x4^3-7*x1^4*x2*x3^3*x4^5+21*x1^4*x2*x3^4*x4^4-21*x1^4*x2*x3^5*x4^3+7*x1^4*x2*x3^6*x4^2+2*x1^4*x2^2*x3^3*x4^4-6*x1^4*x2^2*x3^4*x4^3+6*x1^4*x2^2*x3^5*x4^2-2*x1^4*x2^2*x3^6*x4+x1^4*x2^3*x3^3*x4^3-3*x1^4*x2^3*x3^4*x4^2+3*x1^4*x2^3*x3^5*x4-x1^4*x2^3*x3^6-x1^5*x3^2*x4^6+3*x1^5*x3^3*x4^5-3*x1^5*x3^4*x4^4+x1^5*x3^5*x4^3+2*x1^5*x2*x3^2*x4^5-6*x1^5*x2*x3^3*x4^4+6*x1^5*x2*x3^4*x4^3-2*x1^5*x2*x3^5*x4^2-x1^5*x2^2*x3^2*x4^4+3*x1^5*x2^2*x3^3*x4^3-3*x1^5*x2^2*x3^4*x4^2+x1^5*x2^2*x3^5*x4-x0*x2*x3^6*x4^4-2*x0*x2*x3^6*x4^5+3*x0*x2*x3^7*x4^3+9*x0*x2*x3^7*x4^4-3*x0*x2*x3^8*x4^2-12*x0*x2*x3^8*x4^3+x0*x2*x3^9*x4+5*x0*x2*x3^9*x4^2+x0*x2^2*x3^6*x4^3+4*x0*x2^2*x3^6*x4^4-3*x0*x2^2*x3^7*x4^2-18*x0*x2^2*x3^7*x4^3+3*x0*x2^2*x3^8*x4+24*x0*x2^2*x3^8*x4^2-x0*x2^2*x3^9-10*x0*x2^2*x3^9*x4-2*x0*x2^3*x3^6*x4^3+9*x0*x2^3*x3^7*x4^2-12*x0*x2^3*x3^8*x4+5*x0*x2^3*x3^9+x0*x1*x3^5*x4^5+2*x0*x1*x3^5*x4^6-3*x0*x1*x3^6*x4^4-9*x0*x1*x3^6*x4^5+3*x0*x1*x3^7*x4^3+12*x0*x1*x3^7*x4^4-x0*x1*x3^8*x4^2-5*x0*x1*x3^8*x4^3+2*x0*x1*x2*x3^5*x4^4+4*x0*x1*x2*x3^5*x4^5-6*x0*x1*x2*x3^6*x4^3-18*x0*x1*x2*x3^6*x4^4+6*x0*x1*x2*x3^7*x4^2+24*x0*x1*x2*x3^7*x4^3-2*x0*x1*x2*x3^8*x4-10*x0*x1*x2*x3^8*x4^2-3*x0*x1*x2^2*x3^5*x4^3-14*x0*x1*x2^2*x3^5*x4^4+9*x0*x1*x2^2*x3^6*x4^2+63*x0*x1*x2^2*x3^6*x4^3-9*x0*x1*x2^2*x3^7*x4-84*x0*x1*x2^2*x3^7*x4^2+3*x0*x1*x2^2*x3^8+35*x0*x1*x2^2*x3^8*x4"
					"+8*x0*x1*x2^3*x3^5*x4^3-36*x0*x1*x2^3*x3^6*x4^2+48*x0*x1*x2^3*x3^7*x4-20*x0*x1*x2^3*x3^8-3*x0*x1^2*x3^4*x4^5-8*x0*x1^2*x3^4*x4^6+9*x0*x1^2*x3^5*x4^4+36*x0*x1^2*x3^5*x4^5-9*x0*x1^2*x3^6*x4^3-48*x0*x1^2*x3^6*x4^4+3*x0*x1^2*x3^7*x4^2+20*x0*x1^2*x3^7*x4^3+4*x0*x1^2*x2*x3^4*x4^5-18*x0*x1^2*x2*x3^5*x4^4+24*x0*x1^2*x2*x3^6*x4^3-10*x0*x1^2*x2*x3^7*x4^2+3*x0*x1^2*x2^2*x3^4*x4^3+16*x0*x1^2*x2^2*x3^4*x4^4-9*x0*x1^2*x2^2*x3^5*x4^2-72*x0*x1^2*x2^2*x3^5*x4^3+9*x0*x1^2*x2^2*x3^6*x4+96*x0*x1^2*x2^2*x3^6*x4^2-3*x0*x1^2*x2^2*x3^7-40*x0*x1^2*x2^2*x3^7*x4-12*x0*x1^2*x2^3*x3^4*x4^3+54*x0*x1^2*x2^3*x3^5*x4^2-72*x0*x1^2*x2^3*x3^6*x4+30*x0*x1^2*x2^3*x3^7+3*x0*x1^3*x3^3*x4^5+12*x0*x1^3*x3^3*x4^6-9*x0*x1^3*x3^4*x4^4-54*x0*x1^3*x3^4*x4^5+9*x0*x1^3*x3^5*x4^3+72*x0*x1^3*x3^5*x4^4-3*x0*x1^3*x3^6*x4^2-30*x0*x1^3*x3^6*x4^3-2*x0*x1^3*x2*x3^3*x4^4-16*x0*x1^3*x2*x3^3*x4^5+6*x0*x1^3*x2*x3^4*x4^3+72*x0*x1^3*x2*x3^4*x4^4-6*x0*x1^3*x2*x3^5*x4^2-96*x0*x1^3*x2*x3^5*x4^3+2*x0*x1^3*x2*x3^6*x4+40*x0*x1^3*x2*x3^6*x4^2-x0*x1^3*x2^2*x3^3*x4^3-4*x0*x1^3*x2^2*x3^3*x4^4+3*x0*x1^3*x2^2*x3^4*x4^2+18*x0*x1^3*x2^2*x3^4*x4^3-3*x0*x1^3*x2^2*x3^5*x4-24*x0*x1^3*x2^2*x3^5*x4^2+x0*x1^3*x2^2*x3^6+10*x0*x1^3*x2^2*x3^6*x4+8*x0*x1^3*x2^3*x3^3*x4^3-36*x0*x1^3*x2^3*x3^4*x4^2+48*x0*x1^3*x2^3*x3^5*x4-20*x0*x1^3*x2^3*x3^6-x0*x1^4*x3^2*x4^5-8*x0*x1^4*x3^2*x4^6+3*x0*x1^4*x3^3*x4^4+36*x0*x1^4*x3^3*x4^5-3*x0*x1^4*x3^4*x4^3-48*x0*x1^4*x3^4*x4^4+x0*x1^4*x3^5*x4^2+20*x0*x1^4*x3^5*x4^3+x0*x1^4*x2*x3^2*x4^4+14*x0*x1^4*x2*x3^2*x4^5-3*x0*x1^4*x2*x3^3*x4^3-63*x0*x1^4*x2*x3^3*x4^4+3*x0*x1^4*x2*x3^4*x4^2+84*x0*x1^4*x2*x3^4*x4^3-x0*x1^4*x2*x3^5*x4-35*x0*x1^4*x2*x3^5*x4^2-4*x0*x1^4*x2^2*x3^2*x4^4+18*x0*x1^4*x2^2*x3^3*x4^3-24*x0*x1^4*x2^2*x3^4*x4^2+10*x0*x1^4*x2^2*x3^5*x4-2*x0*x1^4*x2^3*x3^2*x4^3+9*x0*x1^4*x2^3*x3^3*x4^2-12*x0*x1^4*x2^3*x3^4*x4+5*x0*x1^4*x2^3*x3^5+2*x0*x1^5*x3*x4^6-9*x0*x1^5*x3^2*x4^5+12*x0*x1^5*x3^3*x4^4-5*x0*x1^5*x3^4*x4^3-4*x0*x1^5*x2*x3*x4^5+18*x0*x1^5*x2*x3^2*x4^4-24*x0*x1^5*x2*x3^3*x4^3+10*x0*x1^5*x2*x3^4*x4^2+2*x0*x1^5*x2^2*x3*x4^4-9*x0*x1^5*x2^2*x3^2*x4^3+12*x0*x1^5*x2^2*x3^3*x4^2-5*x0*x1^5*x2^2*x3^4*x4+2*x0^2*x2*x3^5*x4^4+x0^2*x2*x3^5*x4^5-9*x0^2*x2*x3^6*x4^3-9*x0^2*x2*x3^6*x4^4+12*x0^2*x2*x3^7*x4^2+18*x0^2*x2*x3^7*x4^3-5*x0^2*x2*x3^8*x4-10*x0^2*x2*x3^8*x4^2-2*x0^2*x2^2*x3^5*x4^3-2*x0^2*x2^2*x3^5*x4^4+9*x0^2*x2^2*x3^6*x4^2+18*x0^2*x2^2*x3^6*x4^3-12*x0^2*x2^2*x3^7*x4-36*x0^2*x2^2*x3^7*x4^2+5*x0^2*x2^2*x3^8+20*x0^2*x2^2*x3^8*x4+x0^2*x2^3*x3^5*x4^3-9*x0^2*x2^3*x3^6*x4^2+18*x0^2*x2^3*x3^7*x4-10*x0^2*x2^3*x3^8-2*x0^2*x1*x3^4*x4^5-x0^2*x1*x3^4*x4^6+9*x0^2*x1*x3^5*x4^4+9*x0^2*x1*x3^5*x4^5-12*x0^2*x1*x3^6*x4^3-18*x0^2*x1*x3^6*x4^4+5*x0^2*x1*x3^7*x4^2+10*x0^2*x1*x3^7*x4^3-4*x0^2*x1*x2*x3^4*x4^4-2*x0^2*x1*x2*x3^4*x4^5+18*x0^2*x1*x2*x3^5*x4^3+18*x0^2*x1*x2*x3^5*x4^4-24*x0^2*x1*x2*x3^6*x4^2-36*x0^2*x1*x2*x3^6*x4^3+10*x0^2*x1*x2*x3^7*x4+20*x0^2*x1*x2*x3^7*x4^2+6*x0^2*x1*x2^2*x3^4*x4^3+7*x0^2*x1*x2^2*x3^4*x4^4-27*x0^2*x1*x2^2*x3^5*x4^2-63*x0^2*x1*x2^2*x3^5*x4^3+36*x0^2*x1*x2^2*x3^6*x4+126*x0^2*x1*x2^2*x3^6*x4^2-15*x0^2*x1*x2^2*x3^7-70*x0^2*x1*x2^2*x3^7*x4-4*x0^2*x1*x2^3*x3^4*x4^3+36*x0^2*x1*x2^3*x3^5*x4^2-72*x0^2*x1*x2^3*x3^6*x4+40*x0^2*x1*x2^3*x3^7+6*x0^2*x1^2*x3^3*x4^5+4*x0^2*x1^2*x3^3*x4^6-27*x0^2*x1^2*x3^4*x4^4-36*x0^2*x1^2*x3^4*x4^5+36*x0^2*x1^2*x3^5*x4^3+72*x0^2*x1^2*x3^5*x4^4-15*x0^2*x1^2*x3^6*x4^2-40*x0^2*x1^2*x3^6*x4^3-2*x0^2*x1^2*x2*x3^3*x4^5+18*x0^2*x1^2*x2*x3^4*x4^4-36*x0^2*x1^2*x2*x3^5*x4^3+20*x0^2*x1^2*x2*x3^6*x4^2-6*x0^2*x1^2*x2^2*x3^3*x4^3-8*x0^2*x1^2*x2^2*x3^3*x4^4+27*x0^2*x1^2*x2^2*x3^4*x4^2+72*x0^2*x1^2*x2^2*x3^4*x4^3-36*x0^2*x1^2*x2^2*x3^5*x4-144*x0^2*x1^2*x2^2*x3^5*x4^2+15*x0^2*x1^2*x2^2*x3^6+80*x0^2*x1^2*x2^2*x3^6*x4+6*x0^2*x1^2*x2^3*x3^3*x4^3-54*x0^2*x1^2*x2^3*x3^4*x4^2+108*x0^2*x1^2*x2^3*x3^5*x4-60*x0^2*x1^2*x2^3*x3^6-6*x0^2*x1^3*x3^2*x4^5-6*x0^2*x1^3*x3^2*x4^6+27*x0^2*x1^3*x3^3*x4^4+54*x0^2*x1^3*x3^3*x4^5-36*x0^2*x1^3*x3^4*x4^3-108*x0^2*x1^3*x3^4*x4^4+15*x0^2*x1^3*x3^5*x4^2+60*x0^2*x1^3*x3^5*x4^3+4*x0^2*x1^3*x2*x3^2*x4^4+8*x0^2*x1^3*x2*x3^2*x4^5-18*x0^2*x1^3*x2*x3^3*x4^3-72*x0^2*x1^3*x2*x3^3*x4^4+24*x0^2*x1^3*x2*x3^4*x4^2+144*x0^2*x1^3*x2*x3^4*x4^3-10*x0^2*x1^3*x2*x3^5*x4"
					"-80*x0^2*x1^3*x2*x3^5*x4^2+2*x0^2*x1^3*x2^2*x3^2*x4^3+2*x0^2*x1^3*x2^2*x3^2*x4^4-9*x0^2*x1^3*x2^2*x3^3*x4^2-18*x0^2*x1^3*x2^2*x3^3*x4^3+12*x0^2*x1^3*x2^2*x3^4*x4+36*x0^2*x1^3*x2^2*x3^4*x4^2-5*x0^2*x1^3*x2^2*x3^5-20*x0^2*x1^3*x2^2*x3^5*x4-4*x0^2*x1^3*x2^3*x3^2*x4^3+36*x0^2*x1^3*x2^3*x3^3*x4^2-72*x0^2*x1^3*x2^3*x3^4*x4+40*x0^2*x1^3*x2^3*x3^5+2*x0^2*x1^4*x3*x4^5+4*x0^2*x1^4*x3*x4^6-9*x0^2*x1^4*x3^2*x4^4-36*x0^2*x1^4*x3^2*x4^5+12*x0^2*x1^4*x3^3*x4^3+72*x0^2*x1^4*x3^3*x4^4-5*x0^2*x1^4*x3^4*x4^2-40*x0^2*x1^4*x3^4*x4^3-2*x0^2*x1^4*x2*x3*x4^4-7*x0^2*x1^4*x2*x3*x4^5+9*x0^2*x1^4*x2*x3^2*x4^3+63*x0^2*x1^4*x2*x3^2*x4^4-12*x0^2*x1^4*x2*x3^3*x4^2-126*x0^2*x1^4*x2*x3^3*x4^3+5*x0^2*x1^4*x2*x3^4*x4+70*x0^2*x1^4*x2*x3^4*x4^2+2*x0^2*x1^4*x2^2*x3*x4^4-18*x0^2*x1^4*x2^2*x3^2*x4^3+36*x0^2*x1^4*x2^2*x3^3*x4^2-20*x0^2*x1^4*x2^2*x3^4*x4+x0^2*x1^4*x2^3*x3*x4^3-9*x0^2*x1^4*x2^3*x3^2*x4^2+18*x0^2*x1^4*x2^3*x3^3*x4-10*x0^2*x1^4*x2^3*x3^4-x0^2*x1^5*x4^6+9*x0^2*x1^5*x3*x4^5-18*x0^2*x1^5*x3^2*x4^4+10*x0^2*x1^5*x3^3*x4^3+2*x0^2*x1^5*x2*x4^5-18*x0^2*x1^5*x2*x3*x4^4+36*x0^2*x1^5*x2*x3^2*x4^3-20*x0^2*x1^5*x2*x3^3*x4^2-x0^2*x1^5*x2^2*x4^4+9*x0^2*x1^5*x2^2*x3*x4^3-18*x0^2*x1^5*x2^2*x3^2*x4^2+10*x0^2*x1^5*x2^2*x3^3*x4-x0^3*x2*x3^4*x4^4+9*x0^3*x2*x3^5*x4^3+3*x0^3*x2*x3^5*x4^4-18*x0^3*x2*x3^6*x4^2-12*x0^3*x2*x3^6*x4^3+10*x0^3*x2*x3^7*x4+10*x0^3*x2*x3^7*x4^2+x0^3*x2^2*x3^4*x4^3-9*x0^3*x2^2*x3^5*x4^2-6*x0^3*x2^2*x3^5*x4^3+18*x0^3*x2^2*x3^6*x4+24*x0^3*x2^2*x3^6*x4^2-10*x0^3*x2^2*x3^7-20*x0^3*x2^2*x3^7*x4+3*x0^3*x2^3*x3^5*x4^2-12*x0^3*x2^3*x3^6*x4+10*x0^3*x2^3*x3^7+x0^3*x1*x3^3*x4^5-9*x0^3*x1*x3^4*x4^4-3*x0^3*x1*x3^4*x4^5+18*x0^3*x1*x3^5*x4^3+12*x0^3*x1*x3^5*x4^4-10*x0^3*x1*x3^6*x4^2-10*x0^3*x1*x3^6*x4^3+2*x0^3*x1*x2*x3^3*x4^4-18*x0^3*x1*x2*x3^4*x4^3-6*x0^3*x1*x2*x3^4*x4^4+36*x0^3*x1*x2*x3^5*x4^2+24*x0^3*x1*x2*x3^5*x4^3-20*x0^3*x1*x2*x3^6*x4-20*x0^3*x1*x2*x3^6*x4^2-3*x0^3*x1*x2^2*x3^3*x4^3+27*x0^3*x1*x2^2*x3^4*x4^2+21*x0^3*x1*x2^2*x3^4*x4^3-54*x0^3*x1*x2^2*x3^5*x4-84*x0^3*x1*x2^2*x3^5*x4^2+30*x0^3*x1*x2^2*x3^6+70*x0^3*x1*x2^2*x3^6*x4-12*x0^3*x1*x2^3*x3^4*x4^2+48*x0^3*x1*x2^3*x3^5*x4-40*x0^3*x1*x2^3*x3^6-3*x0^3*x1^2*x3^2*x4^5+27*x0^3*x1^2*x3^3*x4^4+12*x0^3*x1^2*x3^3*x4^5-54*x0^3*x1^2*x3^4*x4^3-48*x0^3*x1^2*x3^4*x4^4+30*x0^3*x1^2*x3^5*x4^2+40*x0^3*x1^2*x3^5*x4^3-6*x0^3*x1^2*x2*x3^3*x4^4+24*x0^3*x1^2*x2*x3^4*x4^3-20*x0^3*x1^2*x2*x3^5*x4^2+3*x0^3*x1^2*x2^2*x3^2*x4^3-27*x0^3*x1^2*x2^2*x3^3*x4^2-24*x0^3*x1^2*x2^2*x3^3*x4^3+54*x0^3*x1^2*x2^2*x3^4*x4+96*x0^3*x1^2*x2^2*x3^4*x4^2-30*x0^3*x1^2*x2^2*x3^5-80*x0^3*x1^2*x2^2*x3^5*x4+18*x0^3*x1^2*x2^3*x3^3*x4^2-72*x0^3*x1^2*x2^3*x3^4*x4+60*x0^3*x1^2*x2^3*x3^5+3*x0^3*x1^3*x3*x4^5-27*x0^3*x1^3*x3^2*x4^4-18*x0^3*x1^3*x3^2*x4^5+54*x0^3*x1^3*x3^3*x4^3+72*x0^3*x1^3*x3^3*x4^4-30*x0^3*x1^3*x3^4*x4^2-60*x0^3*x1^3*x3^4*x4^3-2*x0^3*x1^3*x2*x3*x4^4+18*x0^3*x1^3*x2*x3^2*x4^3+24*x0^3*x1^3*x2*x3^2*x4^4-36*x0^3*x1^3*x2*x3^3*x4^2-96*x0^3*x1^3*x2*x3^3*x4^3+20*x0^3*x1^3*x2*x3^4*x4+80*x0^3*x1^3*x2*x3^4*x4^2-x0^3*x1^3*x2^2*x3*x4^3+9*x0^3*x1^3*x2^2*x3^2*x4^2+6*x0^3*x1^3*x2^2*x3^2*x4^3-18*x0^3*x1^3*x2^2*x3^3*x4-24*x0^3*x1^3*x2^2*x3^3*x4^2+10*x0^3*x1^3*x2^2*x3^4+20*x0^3*x1^3*x2^2*x3^4*x4-12*x0^3*x1^3*x2^3*x3^2*x4^2+48*x0^3*x1^3*x2^3*x3^3*x4-40*x0^3*x1^3*x2^3*x3^4-x0^3*x1^4*x4^5+9*x0^3*x1^4*x3*x4^4+12*x0^3*x1^4*x3*x4^5-18*x0^3*x1^4*x3^2*x4^3-48*x0^3*x1^4*x3^2*x4^4+10*x0^3*x1^4*x3^3*x4^2+40*x0^3*x1^4*x3^3*x4^3+x0^3*x1^4*x2*x4^4-9*x0^3*x1^4*x2*x3*x4^3-21*x0^3*x1^4*x2*x3*x4^4+18*x0^3*x1^4*x2*x3^2*x4^2+84*x0^3*x1^4*x2*x3^2*x4^3-10*x0^3*x1^4*x2*x3^3*x4-70*x0^3*x1^4*x2*x3^3*x4^2+6*x0^3*x1^4*x2^2*x3*x4^3-24*x0^3*x1^4*x2^2*x3^2*x4^2+20*x0^3*x1^4*x2^2*x3^3*x4+3*x0^3*x1^4*x2^3*x3*x4^2-12*x0^3*x1^4*x2^3*x3^2*x4+10*x0^3*x1^4*x2^3*x3^3-3*x0^3*x1^5*x4^5+12*x0^3*x1^5*x3*x4^4-10*x0^3*x1^5*x3^2*x4^3+6*x0^3*x1^5*x2*x4^4-24*x0^3*x1^5*x2*x3*x4^3+20*x0^3*x1^5*x2*x3^2*x4^2-3*x0^3*x1^5*x2^2*x4^3+12*x0^3*x1^5*x2^2*x3*x4^2-10*x0^3*x1^5*x2^2*x3^2*x4-3*x0^4*x2*x3^4*x4^3+12*x0^4*x2*x3^5*x4^2+3*x0^4*x2*x3^5*x4^3-10*x0^4*x2*x3^6*x4"
					"-5*x0^4*x2*x3^6*x4^2+3*x0^4*x2^2*x3^4*x4^2-12*x0^4*x2^2*x3^5*x4-6*x0^4*x2^2*x3^5*x4^2+10*x0^4*x2^2*x3^6+10*x0^4*x2^2*x3^6*x4+3*x0^4*x2^3*x3^5*x4-5*x0^4*x2^3*x3^6+3*x0^4*x1*x3^3*x4^4-12*x0^4*x1*x3^4*x4^3-3*x0^4*x1*x3^4*x4^4+10*x0^4*x1*x3^5*x4^2+5*x0^4*x1*x3^5*x4^3+6*x0^4*x1*x2*x3^3*x4^3-24*x0^4*x1*x2*x3^4*x4^2-6*x0^4*x1*x2*x3^4*x4^3+20*x0^4*x1*x2*x3^5*x4+10*x0^4*x1*x2*x3^5*x4^2-9*x0^4*x1*x2^2*x3^3*x4^2+36*x0^4*x1*x2^2*x3^4*x4+21*x0^4*x1*x2^2*x3^4*x4^2-30*x0^4*x1*x2^2*x3^5-35*x0^4*x1*x2^2*x3^5*x4-12*x0^4*x1*x2^3*x3^4*x4+20*x0^4*x1*x2^3*x3^5-9*x0^4*x1^2*x3^2*x4^4+36*x0^4*x1^2*x3^3*x4^3+12*x0^4*x1^2*x3^3*x4^4-30*x0^4*x1^2*x3^4*x4^2-20*x0^4*x1^2*x3^4*x4^3-6*x0^4*x1^2*x2*x3^3*x4^3+10*x0^4*x1^2*x2*x3^4*x4^2+9*x0^4*x1^2*x2^2*x3^2*x4^2-36*x0^4*x1^2*x2^2*x3^3*x4-24*x0^4*x1^2*x2^2*x3^3*x4^2+30*x0^4*x1^2*x2^2*x3^4+40*x0^4*x1^2*x2^2*x3^4*x4+18*x0^4*x1^2*x2^3*x3^3*x4-30*x0^4*x1^2*x2^3*x3^4+9*x0^4*x1^3*x3*x4^4-36*x0^4*x1^3*x3^2*x4^3-18*x0^4*x1^3*x3^2*x4^4+30*x0^4*x1^3*x3^3*x4^2+30*x0^4*x1^3*x3^3*x4^3-6*x0^4*x1^3*x2*x3*x4^3+24*x0^4*x1^3*x2*x3^2*x4^2+24*x0^4*x1^3*x2*x3^2*x4^3-20*x0^4*x1^3*x2*x3^3*x4-40*x0^4*x1^3*x2*x3^3*x4^2-3*x0^4*x1^3*x2^2*x3*x4^2+12*x0^4*x1^3*x2^2*x3^2*x4+6*x0^4*x1^3*x2^2*x3^2*x4^2-10*x0^4*x1^3*x2^2*x3^3-10*x0^4*x1^3*x2^2*x3^3*x4-12*x0^4*x1^3*x2^3*x3^2*x4+20*x0^4*x1^3*x2^3*x3^3-3*x0^4*x1^4*x4^4+12*x0^4*x1^4*x3*x4^3+12*x0^4*x1^4*x3*x4^4-10*x0^4*x1^4*x3^2*x4^2-20*x0^4*x1^4*x3^2*x4^3+3*x0^4*x1^4*x2*x4^3-12*x0^4*x1^4*x2*x3*x4^2-21*x0^4*x1^4*x2*x3*x4^3+10*x0^4*x1^4*x2*x3^2*x4+35*x0^4*x1^4*x2*x3^2*x4^2+6*x0^4*x1^4*x2^2*x3*x4^2-10*x0^4*x1^4*x2^2*x3^2*x4+3*x0^4*x1^4*x2^3*x3*x4-5*x0^4*x1^4*x2^3*x3^2-3*x0^4*x1^5*x4^4+5*x0^4*x1^5*x3*x4^3+6*x0^4*x1^5*x2*x4^3-10*x0^4*x1^5*x2*x3*x4^2-3*x0^4*x1^5*x2^2*x4^2+5*x0^4*x1^5*x2^2*x3*x4-3*x0^5*x2*x3^4*x4^2+5*x0^5*x2*x3^5*x4+x0^5*x2*x3^5*x4^2+3*x0^5*x2^2*x3^4*x4-5*x0^5*x2^2*x3^5-2*x0^5*x2^2*x3^5*x4+x0^5*x2^3*x3^5+3*x0^5*x1*x3^3*x4^3-5*x0^5*x1*x3^4*x4^2-x0^5*x1*x3^4*x4^3+6*x0^5*x1*x2*x3^3*x4^2-10*x0^5*x1*x2*x3^4*x4-2*x0^5*x1*x2*x3^4*x4^2-9*x0^5*x1*x2^2*x3^3*x4+15*x0^5*x1*x2^2*x3^4+7*x0^5*x1*x2^2*x3^4*x4-4*x0^5*x1*x2^3*x3^4-9*x0^5*x1^2*x3^2*x4^3+15*x0^5*x1^2*x3^3*x4^2+4*x0^5*x1^2*x3^3*x4^3-2*x0^5*x1^2*x2*x3^3*x4^2+9*x0^5*x1^2*x2^2*x3^2*x4-15*x0^5*x1^2*x2^2*x3^3-8*x0^5*x1^2*x2^2*x3^3*x4+6*x0^5*x1^2*x2^3*x3^3+9*x0^5*x1^3*x3*x4^3-15*x0^5*x1^3*x3^2*x4^2-6*x0^5*x1^3*x3^2*x4^3-6*x0^5*x1^3*x2*x3*x4^2+10*x0^5*x1^3*x2*x3^2*x4+8*x0^5*x1^3*x2*x3^2*x4^2-3*x0^5*x1^3*x2^2*x3*x4+5*x0^5*x1^3*x2^2*x3^2+2*x0^5*x1^3*x2^2*x3^2*x4-4*x0^5*x1^3*x2^3*x3^2-3*x0^5*x1^4*x4^3+5*x0^5*x1^4*x3*x4^2+4*x0^5*x1^4*x3*x4^3+3*x0^5*x1^4*x2*x4^2-5*x0^5*x1^4*x2*x3*x4-7*x0^5*x1^4*x2*x3*x4^2+2*x0^5*x1^4*x2^2*x3*x4+x0^5*x1^4*x2^3*x3-x0^5*x1^5*x4^3+2*x0^5*x1^5*x2*x4^2-x0^5*x1^5*x2^2*x4-x0^6*x2*x3^4*x4+x0^6*x2^2*x3^4+x0^6*x1*x3^3*x4^2+2*x0^6*x1*x2*x3^3*x4-3*x0^6*x1*x2^2*x3^3-3*x0^6*x1^2*x3^2*x4^2+3*x0^6*x1^2*x2^2*x3^2+3*x0^6*x1^3*x3*x4^2-2*x0^6*x1^3*x2*x3*x4-x0^6*x1^3*x2^2*x3-x0^6*x1^4*x4^2+x0^6*x1^4*x2*x4,"
					"-x2*x3^8*x4^4+3*x2*x3^9*x4^3-3*x2*x3^10*x4^2+x2*x3^11*x4+x2^2*x3^8*x4^3-3*x2^2*x3^9*x4^2+3*x2^2*x3^10*x4-x2^2*x3^11+x1*x3^7*x4^5-3*x1*x3^8*x4^4+3*x1*x3^9*x4^3-x1*x3^10*x4^2+4*x1*x2*x3^7*x4^4-12*x1*x2*x3^8*x4^3+12*x1*x2*x3^9*x4^2-4*x1*x2*x3^10*x4-5*x1*x2^2*x3^7*x4^3+15*x1*x2^2*x3^8*x4^2-15*x1*x2^2*x3^9*x4+5*x1*x2^2*x3^10-5*x1^2*x3^6*x4^5+15*x1^2*x3^7*x4^4-15*x1^2*x3^8*x4^3+5*x1^2*x3^9*x4^2-5*x1^2*x2*x3^6*x4^4+15*x1^2*x2*x3^7*x4^3-15*x1^2*x2*x3^8*x4^2+5*x1^2*x2*x3^9*x4+10*x1^2*x2^2*x3^6*x4^3-30*x1^2*x2^2*x3^7*x4^2+30*x1^2*x2^2*x3^8*x4-10*x1^2*x2^2*x3^9+10*x1^3*x3^5*x4^5"
					"-30*x1^3*x3^6*x4^4+30*x1^3*x3^7*x4^3-10*x1^3*x3^8*x4^2-10*x1^3*x2^2*x3^5*x4^3+30*x1^3*x2^2*x3^6*x4^2-30*x1^3*x2^2*x3^7*x4+10*x1^3*x2^2*x3^8-10*x1^4*x3^4*x4^5+30*x1^4*x3^5*x4^4-30*x1^4*x3^6*x4^3+10*x1^4*x3^7*x4^2+5*x1^4*x2*x3^4*x4^4-15*x1^4*x2*x3^5*x4^3+15*x1^4*x2*x3^6*x4^2-5*x1^4*x2*x3^7*x4+5*x1^4*x2^2*x3^4*x4^3-15*x1^4*x2^2*x3^5*x4^2+15*x1^4*x2^2*x3^6*x4-5*x1^4*x2^2*x3^7+5*x1^5*x3^3*x4^5-15*x1^5*x3^4*x4^4+15*x1^5*x3^5*x4^3-5*x1^5*x3^6*x4^2-4*x1^5*x2*x3^3*x4^4+12*x1^5*x2*x3^4*x4^3-12*x1^5*x2*x3^5*x4^2+4*x1^5*x2*x3^6*x4-x1^5*x2^2*x3^3*x4^3+3*x1^5*x2^2*x3^4*x4^2-3*x1^5*x2^2*x3^5*x4+x1^5*x2^2*x3^6-x1^6*x3^2*x4^5+3*x1^6*x3^3*x4^4-3*x1^6*x3^4*x4^3+x1^6*x3^5*x4^2+x1^6*x2*x3^2*x4^4-3*x1^6*x2*x3^3*x4^3+3*x1^6*x2*x3^4*x4^2-x1^6*x2*x3^5*x4+x0*x2*x3^7*x4^3+2*x0*x2*x3^7*x4^4-3*x0*x2*x3^8*x4^2-9*x0*x2*x3^8*x4^3+3*x0*x2*x3^9*x4+12*x0*x2*x3^9*x4^2-x0*x2*x3^10-5*x0*x2*x3^10*x4-2*x0*x2^2*x3^7*x4^3+9*x0*x2^2*x3^8*x4^2-12*x0*x2^2*x3^9*x4+5*x0*x2^2*x3^10-x0*x1*x3^6*x4^4-2*x0*x1*x3^6*x4^5+3*x0*x1*x3^7*x4^3+9*x0*x1*x3^7*x4^4-3*x0*x1*x3^8*x4^2-12*x0*x1*x3^8*x4^3+x0*x1*x3^9*x4+5*x0*x1*x3^9*x4^2-4*x0*x1*x2*x3^6*x4^3-8*x0*x1*x2*x3^6*x4^4+12*x0*x1*x2*x3^7*x4^2+36*x0*x1*x2*x3^7*x4^3-12*x0*x1*x2*x3^8*x4-48*x0*x1*x2*x3^8*x4^2+4*x0*x1*x2*x3^9+20*x0*x1*x2*x3^9*x4+10*x0*x1*x2^2*x3^6*x4^3-45*x0*x1*x2^2*x3^7*x4^2+60*x0*x1*x2^2*x3^8*x4-25*x0*x1*x2^2*x3^9+4*x0*x1^2*x3^5*x4^4+10*x0*x1^2*x3^5*x4^5-12*x0*x1^2*x3^6*x4^3-45*x0*x1^2*x3^6*x4^4+12*x0*x1^2*x3^7*x4^2+60*x0*x1^2*x3^7*x4^3-4*x0*x1^2*x3^8*x4-25*x0*x1^2*x3^8*x4^2+6*x0*x1^2*x2*x3^5*x4^3+10*x0*x1^2*x2*x3^5*x4^4-18*x0*x1^2*x2*x3^6*x4^2-45*x0*x1^2*x2*x3^6*x4^3+18*x0*x1^2*x2*x3^7*x4+60*x0*x1^2*x2*x3^7*x4^2-6*x0*x1^2*x2*x3^8-25*x0*x1^2*x2*x3^8*x4-20*x0*x1^2*x2^2*x3^5*x4^3+90*x0*x1^2*x2^2*x3^6*x4^2-120*x0*x1^2*x2^2*x3^7*x4+50*x0*x1^2*x2^2*x3^8-6*x0*x1^3*x3^4*x4^4-20*x0*x1^3*x3^4*x4^5+18*x0*x1^3*x3^5*x4^3+90*x0*x1^3*x3^5*x4^4-18*x0*x1^3*x3^6*x4^2-120*x0*x1^3*x3^6*x4^3+6*x0*x1^3*x3^7*x4+50*x0*x1^3*x3^7*x4^2-4*x0*x1^3*x2*x3^4*x4^3+12*x0*x1^3*x2*x3^5*x4^2-12*x0*x1^3*x2*x3^6*x4+4*x0*x1^3*x2*x3^7+20*x0*x1^3*x2^2*x3^4*x4^3-90*x0*x1^3*x2^2*x3^5*x4^2+120*x0*x1^3*x2^2*x3^6*x4-50*x0*x1^3*x2^2*x3^7+4*x0*x1^4*x3^3*x4^4+20*x0*x1^4*x3^3*x4^5-12*x0*x1^4*x3^4*x4^3-90*x0*x1^4*x3^4*x4^4+12*x0*x1^4*x3^5*x4^2+120*x0*x1^4*x3^5*x4^3-4*x0*x1^4*x3^6*x4-50*x0*x1^4*x3^6*x4^2+x0*x1^4*x2*x3^3*x4^3-10*x0*x1^4*x2*x3^3*x4^4-3*x0*x1^4*x2*x3^4*x4^2+45*x0*x1^4*x2*x3^4*x4^3+3*x0*x1^4*x2*x3^5*x4-60*x0*x1^4*x2*x3^5*x4^2-x0*x1^4*x2*x3^6+25*x0*x1^4*x2*x3^6*x4-10*x0*x1^4*x2^2*x3^3*x4^3+45*x0*x1^4*x2^2*x3^4*x4^2-60*x0*x1^4*x2^2*x3^5*x4+25*x0*x1^4*x2^2*x3^6-x0*x1^5*x3^2*x4^4-10*x0*x1^5*x3^2*x4^5+3*x0*x1^5*x3^3*x4^3+45*x0*x1^5*x3^3*x4^4-3*x0*x1^5*x3^4*x4^2-60*x0*x1^5*x3^4*x4^3+x0*x1^5*x3^5*x4+25*x0*x1^5*x3^5*x4^2+8*x0*x1^5*x2*x3^2*x4^4-36*x0*x1^5*x2*x3^3*x4^3+48*x0*x1^5*x2*x3^4*x4^2-20*x0*x1^5*x2*x3^5*x4+2*x0*x1^5*x2^2*x3^2*x4^3-9*x0*x1^5*x2^2*x3^3*x4^2+12*x0*x1^5*x2^2*x3^4*x4-5*x0*x1^5*x2^2*x3^5+2*x0*x1^6*x3*x4^5-9*x0*x1^6*x3^2*x4^4+12*x0*x1^6*x3^3*x4^3-5*x0*x1^6*x3^4*x4^2-2*x0*x1^6*x2*x3*x4^4+9*x0*x1^6*x2*x3^2*x4^3-12*x0*x1^6*x2*x3^3*x4^2+5*x0*x1^6*x2*x3^4*x4-2*x0^2*x2*x3^6*x4^3-x0^2*x2*x3^6*x4^4+9*x0^2*x2*x3^7*x4^2+9*x0^2*x2*x3^7*x4^3-12*x0^2*x2*x3^8*x4-18*x0^2*x2*x3^8*x4^2+5*x0^2*x2*x3^9+10*x0^2*x2*x3^9*x4+x0^2*x2^2*x3^6*x4^3-9*x0^2*x2^2*x3^7*x4^2+18*x0^2*x2^2*x3^8*x4-10*x0^2*x2^2*x3^9+2*x0^2*x1*x3^5*x4^4+x0^2*x1*x3^5*x4^5-9*x0^2*x1*x3^6*x4^3-9*x0^2*x1*x3^6*x4^4+12*x0^2*x1*x3^7*x4^2+18*x0^2*x1*x3^7*x4^3-5*x0^2*x1*x3^8*x4-10*x0^2*x1*x3^8*x4^2+8*x0^2*x1*x2*x3^5*x4^3+4*x0^2*x1*x2*x3^5*x4^4-36*x0^2*x1*x2*x3^6*x4^2-36*x0^2*x1*x2*x3^6*x4^3+48*x0^2*x1*x2*x3^7*x4+72*x0^2*x1*x2*x3^7*x4^2-20*x0^2*x1*x2*x3^8-40*x0^2*x1*x2*x3^8*x4-5*x0^2*x1*x2^2*x3^5*x4^3+45*x0^2*x1*x2^2*x3^6*x4^2-90*x0^2*x1*x2^2*x3^7*x4+50*x0^2*x1*x2^2*x3^8-8*x0^2*x1^2*x3^4*x4^4-5*x0^2*x1^2*x3^4*x4^5+36*x0^2*x1^2*x3^5*x4^3+45*x0^2*x1^2*x3^5*x4^4-48*x0^2*x1^2*x3^6*x4^2-90*x0^2*x1^2*x3^6*x4^3+20*x0^2*x1^2*x3^7*x4+50*x0^2*x1^2*x3^7*x4^2-12*x0^2*x1^2*x2*x3^4*x4^3-5*x0^2*x1^2*x2*x3^4*x4^4+54*x0^2*x1^2*x2*x3^5*x4^2+45*x0^2*x1^2*x2*x3^5*x4^3-72*x0^2*x1^2*x2*x3^6*x4-90*x0^2*x1^2*x2*x3^6*x4^2+30*x0^2*x1^2*x2*x3^7+50*x0^2*x1^2*x2*x3^7*x4+10*x0^2*x1^2*x2^2*x3^4*x4^3-90*x0^2*x1^2*x2^2*x3^5*x4^2+180*x0^2*x1^2*x2^2*x3^6*x4-100*x0^2*x1^2*x2^2*x3^7+12*x0^2*x1^3*x3^3*x4^4+10*x0^2*x1^3*x3^3*x4^5-54*x0^2*x1^3*x3^4*x4^3-90*x0^2*x1^3*x3^4*x4^4+72*x0^2*x1^3*x3^5*x4^2+180*x0^2*x1^3*x3^5*x4^3-30*x0^2*x1^3*x3^6*x4-100*x0^2*x1^3*x3^6*x4^2+8*x0^2*x1^3*x2*x3^3*x4^3-36*x0^2*x1^3*x2*x3^4*x4^2+48*x0^2*x1^3*x2*x3^5*x4-20*x0^2*x1^3*x2*x3^6-10*x0^2*x1^3*x2^2*x3^3*x4^3+90*x0^2*x1^3*x2^2*x3^4*x4^2-180*x0^2*x1^3*x2^2*x3^5*x4+100*x0^2*x1^3*x2^2*x3^6-8*x0^2*x1^4*x3^2*x4^4-10*x0^2*x1^4*x3^2*x4^5+36*x0^2*x1^4*x3^3*x4^3+90*x0^2*x1^4*x3^3*x4^4-48*x0^2*x1^4*x3^4*x4^2"
					"-180*x0^2*x1^4*x3^4*x4^3+20*x0^2*x1^4*x3^5*x4+100*x0^2*x1^4*x3^5*x4^2-2*x0^2*x1^4*x2*x3^2*x4^3+5*x0^2*x1^4*x2*x3^2*x4^4+9*x0^2*x1^4*x2*x3^3*x4^2-45*x0^2*x1^4*x2*x3^3*x4^3-12*x0^2*x1^4*x2*x3^4*x4+90*x0^2*x1^4*x2*x3^4*x4^2+5*x0^2*x1^4*x2*x3^5-50*x0^2*x1^4*x2*x3^5*x4+5*x0^2*x1^4*x2^2*x3^2*x4^3-45*x0^2*x1^4*x2^2*x3^3*x4^2+90*x0^2*x1^4*x2^2*x3^4*x4-50*x0^2*x1^4*x2^2*x3^5+2*x0^2*x1^5*x3*x4^4+5*x0^2*x1^5*x3*x4^5-9*x0^2*x1^5*x3^2*x4^3-45*x0^2*x1^5*x3^2*x4^4+12*x0^2*x1^5*x3^3*x4^2+90*x0^2*x1^5*x3^3*x4^3-5*x0^2*x1^5*x3^4*x4-50*x0^2*x1^5*x3^4*x4^2-4*x0^2*x1^5*x2*x3*x4^4+36*x0^2*x1^5*x2*x3^2*x4^3-72*x0^2*x1^5*x2*x3^3*x4^2+40*x0^2*x1^5*x2*x3^4*x4-x0^2*x1^5*x2^2*x3*x4^3+9*x0^2*x1^5*x2^2*x3^2*x4^2-18*x0^2*x1^5*x2^2*x3^3*x4+10*x0^2*x1^5*x2^2*x3^4-x0^2*x1^6*x4^5+9*x0^2*x1^6*x3*x4^4-18*x0^2*x1^6*x3^2*x4^3+10*x0^2*x1^6*x3^3*x4^2+x0^2*x1^6*x2*x4^4-9*x0^2*x1^6*x2*x3*x4^3+18*x0^2*x1^6*x2*x3^2*x4^2-10*x0^2*x1^6*x2*x3^3*x4+x0^3*x2*x3^5*x4^3-9*x0^3*x2*x3^6*x4^2-3*x0^3*x2*x3^6*x4^3+18*x0^3*x2*x3^7*x4+12*x0^3*x2*x3^7*x4^2-10*x0^3*x2*x3^8-10*x0^3*x2*x3^8*x4+3*x0^3*x2^2*x3^6*x4^2-12*x0^3*x2^2*x3^7*x4+10*x0^3*x2^2*x3^8-x0^3*x1*x3^4*x4^4+9*x0^3*x1*x3^5*x4^3+3*x0^3*x1*x3^5*x4^4-18*x0^3*x1*x3^6*x4^2-12*x0^3*x1*x3^6*x4^3+10*x0^3*x1*x3^7*x4+10*x0^3*x1*x3^7*x4^2-4*x0^3*x1*x2*x3^4*x4^3+36*x0^3*x1*x2*x3^5*x4^2+12*x0^3*x1*x2*x3^5*x4^3-72*x0^3*x1*x2*x3^6*x4-48*x0^3*x1*x2*x3^6*x4^2+40*x0^3*x1*x2*x3^7+40*x0^3*x1*x2*x3^7*x4-15*x0^3*x1*x2^2*x3^5*x4^2+60*x0^3*x1*x2^2*x3^6*x4-50*x0^3*x1*x2^2*x3^7+4*x0^3*x1^2*x3^3*x4^4-36*x0^3*x1^2*x3^4*x4^3-15*x0^3*x1^2*x3^4*x4^4+72*x0^3*x1^2*x3^5*x4^2+60*x0^3*x1^2*x3^5*x4^3-40*x0^3*x1^2*x3^6*x4-50*x0^3*x1^2*x3^6*x4^2+6*x0^3*x1^2*x2*x3^3*x4^3-54*x0^3*x1^2*x2*x3^4*x4^2-15*x0^3*x1^2*x2*x3^4*x4^3+108*x0^3*x1^2*x2*x3^5*x4+60*x0^3*x1^2*x2*x3^5*x4^2-60*x0^3*x1^2*x2*x3^6-50*x0^3*x1^2*x2*x3^6*x4+30*x0^3*x1^2*x2^2*x3^4*x4^2-120*x0^3*x1^2*x2^2*x3^5*x4+100*x0^3*x1^2*x2^2*x3^6-6*x0^3*x1^3*x3^2*x4^4+54*x0^3*x1^3*x3^3*x4^3+30*x0^3*x1^3*x3^3*x4^4-108*x0^3*x1^3*x3^4*x4^2-120*x0^3*x1^3*x3^4*x4^3+60*x0^3*x1^3*x3^5*x4+100*x0^3*x1^3*x3^5*x4^2-4*x0^3*x1^3*x2*x3^2*x4^3"
					"+36*x0^3*x1^3*x2*x3^3*x4^2-72*x0^3*x1^3*x2*x3^4*x4+40*x0^3*x1^3*x2*x3^5-30*x0^3*x1^3*x2^2*x3^3*x4^2+120*x0^3*x1^3*x2^2*x3^4*x4-100*x0^3*x1^3*x2^2*x3^5+4*x0^3*x1^4*x3*x4^4-36*x0^3*x1^4*x3^2*x4^3-30*x0^3*x1^4*x3^2*x4^4+72*x0^3*x1^4*x3^3*x4^2+120*x0^3*x1^4*x3^3*x4^3-40*x0^3*x1^4*x3^4*x4-100*x0^3*x1^4*x3^4*x4^2+x0^3*x1^4*x2*x3*x4^3-9*x0^3*x1^4*x2*x3^2*x4^2+15*x0^3*x1^4*x2*x3^2*x4^3+18*x0^3*x1^4*x2*x3^3*x4-60*x0^3*x1^4*x2*x3^3*x4^2-10*x0^3*x1^4*x2*x3^4+50*x0^3*x1^4*x2*x3^4*x4+15*x0^3*x1^4*x2^2*x3^2*x4^2-60*x0^3*x1^4*x2^2*x3^3*x4+50*x0^3*x1^4*x2^2*x3^4-x0^3*x1^5*x4^4+9*x0^3*x1^5*x3*x4^3+15*x0^3*x1^5*x3*x4^4-18*x0^3*x1^5*x3^2*x4^2-60*x0^3*x1^5*x3^2*x4^3+10*x0^3*x1^5*x3^3*x4+50*x0^3*x1^5*x3^3*x4^2-12*x0^3*x1^5*x2*x3*x4^3+48*x0^3*x1^5*x2*x3^2*x4^2-40*x0^3*x1^5*x2*x3^3*x4-3*x0^3*x1^5*x2^2*x3*x4^2+12*x0^3*x1^5*x2^2*x3^2*x4-10*x0^3*x1^5*x2^2*x3^3-3*x0^3*x1^6*x4^4+12*x0^3*x1^6*x3*x4^3-10*x0^3*x1^6*x3^2*x4^2+3*x0^3*x1^6*x2*x4^3-12*x0^3*x1^6*x2*x3*x4^2+10*x0^3*x1^6*x2*x3^2*x4+3*x0^4*x2*x3^5*x4^2-12*x0^4*x2*x3^6*x4-3*x0^4*x2*x3^6*x4^2+10*x0^4*x2*x3^7+5*x0^4*x2*x3^7*x4+3*x0^4*x2^2*x3^6*x4-5*x0^4*x2^2*x3^7-3*x0^4*x1*x3^4*x4^3+12*x0^4*x1*x3^5*x4^2+3*x0^4*x1*x3^5*x4^3-10*x0^4*x1*x3^6*x4-5*x0^4*x1*x3^6*x4^2-12*x0^4*x1*x2*x3^4*x4^2+48*x0^4*x1*x2*x3^5*x4+12*x0^4*x1*x2*x3^5*x4^2-40*x0^4*x1*x2*x3^6-20*x0^4*x1*x2*x3^6*x4-15*x0^4*x1*x2^2*x3^5*x4+25*x0^4*x1*x2^2*x3^6+12*x0^4*x1^2*x3^3*x4^3-48*x0^4*x1^2*x3^4*x4^2-15*x0^4*x1^2*x3^4*x4^3+40*x0^4*x1^2*x3^5*x4+25*x0^4*x1^2*x3^5*x4^2+18*x0^4*x1^2*x2*x3^3*x4^2-72*x0^4*x1^2*x2*x3^4*x4-15*x0^4*x1^2*x2*x3^4*x4^2+60*x0^4*x1^2*x2*x3^5+25*x0^4*x1^2*x2*x3^5*x4+30*x0^4*x1^2*x2^2*x3^4*x4-50*x0^4*x1^2*x2^2*x3^5-18*x0^4*x1^3*x3^2*x4^3+72*x0^4*x1^3*x3^3*x4^2+30*x0^4*x1^3*x3^3*x4^3-60*x0^4*x1^3*x3^4*x4-50*x0^4*x1^3*x3^4*x4^2-12*x0^4*x1^3*x2*x3^2*x4^2+48*x0^4*x1^3*x2*x3^3*x4-40*x0^4*x1^3*x2*x3^4-30*x0^4*x1^3*x2^2*x3^3*x4+50*x0^4*x1^3*x2^2*x3^4+12*x0^4*x1^4*x3*x4^3-48*x0^4*x1^4*x3^2*x4^2-30*x0^4*x1^4*x3^2*x4^3+40*x0^4*x1^4*x3^3*x4+50*x0^4*x1^4*x3^3*x4^2+3*x0^4*x1^4*x2*x3*x4^2-12*x0^4*x1^4*x2*x3^2*x4+15*x0^4*x1^4*x2*x3^2*x4^2+10*x0^4*x1^4*x2*x3^3-25*x0^4*x1^4*x2*x3^3*x4+15*x0^4*x1^4*x2^2*x3^2*x4-25*x0^4*x1^4*x2^2*x3^3-3*x0^4*x1^5*x4^3+12*x0^4*x1^5*x3*x4^2+15*x0^4*x1^5*x3*x4^3-10*x0^4*x1^5*x3^2*x4-25*x0^4*x1^5*x3^2*x4^2-12*x0^4*x1^5*x2*x3*x4^2+20*x0^4*x1^5*x2*x3^2*x4-3*x0^4*x1^5*x2^2*x3*x4+5*x0^4*x1^5*x2^2*x3^2-3*x0^4*x1^6*x4^3+5*x0^4*x1^6*x3*x4^2+3*x0^4*x1^6*x2*x4^2-5*x0^4*x1^6*x2*x3*x4+3*x0^5*x2*x3^5*x4-5*x0^5*x2*x3^6-x0^5*x2*x3^6*x4+x0^5*x2^2*x3^6-3*x0^5*x1*x3^4*x4^2+5*x0^5*x1*x3^5*x4+x0^5*x1*x3^5*x4^2-12*x0^5*x1*x2*x3^4*x4+20*x0^5*x1*x2*x3^5+4*x0^5*x1*x2*x3^5*x4-5*x0^5*x1*x2^2*x3^5+12*x0^5*x1^2*x3^3*x4^2-20*x0^5*x1^2*x3^4*x4-5*x0^5*x1^2*x3^4*x4^2+18*x0^5*x1^2*x2*x3^3*x4-30*x0^5*x1^2*x2*x3^4-5*x0^5*x1^2*x2*x3^4*x4+10*x0^5*x1^2*x2^2*x3^4-18*x0^5*x1^3*x3^2*x4^2+30*x0^5*x1^3*x3^3*x4+10*x0^5*x1^3*x3^3*x4^2-12*x0^5*x1^3*x2*x3^2*x4+20*x0^5*x1^3*x2*x3^3-10*x0^5*x1^3*x2^2*x3^3+12*x0^5*x1^4*x3*x4^2-20*x0^5*x1^4*x3^2*x4-10*x0^5*x1^4*x3^2*x4^2+3*x0^5*x1^4*x2*x3*x4-5*x0^5*x1^4*x2*x3^2+5*x0^5*x1^4*x2*x3^2*x4+5*x0^5*x1^4*x2^2*x3^2-3*x0^5*x1^5*x4^2+5*x0^5*x1^5*x3*x4+5*x0^5*x1^5*x3*x4^2-4*x0^5*x1^5*x2*x3*x4-x0^5*x1^5*x2^2*x3-x0^5*x1^6*x4^2+x0^5*x1^6*x2*x4+x0^6*x2*x3^5-x0^6*x1*x3^4*x4-4*x0^6*x1*x2*x3^4+4*x0^6*x1^2*x3^3*x4+6*x0^6*x1^2*x2*x3^3-6*x0^6*x1^3*x3^2*x4-4*x0^6*x1^3*x2*x3^2+4*x0^6*x1^4*x3*x4+x0^6*x1^4*x2*x3-x0^6*x1^5*x4}"
#endif
//					"Z/9851Z[x0,x1,x2]"
//					"{-x2+2*x2^2-x2^3+x1+9848*x1*x2+2*x1*x2^2+x1^2-x1^2*x2-x0+3*x0*x2+9849*x0*x2^2+9849*x0*x1+2*x0*x1*x2+x0^2-x0^2*x2,"
//					"x2+9849*x2^2+x2^3+x1*x2-x1*x2^2+x0+9849*x0*x2+x0*x2^2}"

//					"Z/9851Z[x0,x2]"
//					"{8644+9315*x2+1200*x2^2+8445*x2^3+7996*x0+2194*x0*x2+1415*x0*x2^2-x0*x2^3+994*x0^2+7048*x0^2*x2+9848*x0^2*x2^2+4227*x0^3+9848*x0^3*x2-x0^4,"
//					"1319+206*x2+4218*x2^2+x2^3+206*x0+8436*x0*x2+3*x0*x2^2+4218*x0^2+3*x0^2*x2+x0^3}"
//					"{3700*x1+9598*x1^2-x1^3+6025*x0+505*x0*x1+3*x0*x1^2+9599*x0^2+9848*x0^2*x1+x0^3,"
//					"9228+1627*x1+9473*x1^2-x1^3+8224*x0+756*x0*x1+3*x0*x1^2+9473*x0^2+9848*x0^2*x1+x0^3}"
					).parsePolynomialSetWithRing();
			A.mark(LexicographicTermOrder());
			debug<<homogeneitySpace(A);
			debug<<NonMonomialPolynomialGCDForZModP(A);
			A.pop_front();
			debug<<homogeneitySpace(A);
		}
		else
		if(0)
		{
			PolynomialSet p=StringParser("{1+a+b+ab,1+a+bb+abb}").parsePolynomialSet(r);

			debug<<p<<"\n";

			debug<<NonMonomialPolynomialGCDForZModP(p);
		}
		else if(1)
		{
			debug<<"A1\n";
			PolynomialRing r=StringParser("Q[a,b]").parsePolynomialRing();
			Polynomial p=StringParser("1+a+3b+2abb").parsePolynomial(r);
			Polynomial q=StringParser("2a+b+2aab").parsePolynomial(r);
			Polynomial a=StringParser("1/2+2/2a+3/2b+4/2ab").parsePolynomial(r);
			PolynomialSet A(r);
			A.push_back(p*a);
			A.push_back(q*a);
			debug<<A<<"\n";

			debug<<NonMonomialPolynomialGCDForZ(A);

		}
		else
		{
			Polynomial p=StringParser("1+a+3b+2abb").parsePolynomial(r);
			Polynomial q=StringParser("2a+b+2aab").parsePolynomial(r);
			Polynomial a=StringParser("1+2a+3b+4ab").parsePolynomial(r);
//			Polynomial p=StringParser("1+3b").parsePolynomial(r);
//			Polynomial q=StringParser("3a+b").parsePolynomial(r);
//			Polynomial a=StringParser("1+a+b").parsePolynomial(r);
			PolynomialSet A(r);
			A.push_back(p*a);
			A.push_back(q*a);
			debug<<A<<"\n";

			debug<<NonMonomialPolynomialGCDForZModP(A);

		}
#endif
		return 0;
	}

  void lpRationalFunctionTest()
  {
    //int n=3;
	//     IntegerVectorList L=StringParser("{(1,0,2),(1,2,0)}").parseIntegerVectorList();
    //  IntegerVectorList L=StringParser("{(1,2,0)}").parseIntegerVectorList();

    //     int n=4;
    //  IntegerVectorList L=StringParser("{(1,0,2,3),(1,2,3,0)}").parseIntegerVectorList();

    /*     int n=16;
     IntegerVectorList L=StringParser("{(1,0,2,3,5,4,6,7,9,8,10,11,13,12,14,15),"
				      "(3,0,1,2,7,4,5,6,11,8,9,10,15,12,13,14),"
				      "(0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15)}").parseIntegerVectorList();
    */
    int n=9;
     IntegerVectorList L=StringParser("{(1,0,2,4,3,5,7,6,8),"
				      "(1,2,0,4,5,3,7,8,6),"
				      "(3,4,5,0,1,2,6,7,8),"
				      "(3,4,5,6,7,8,0,1,2)"
				      "}").parseIntegerVectorList();
     /*
    int n=9;
     IntegerVectorList L=StringParser("{(1,0,2,4,3,5,7,6,8),"
				      "(1,2,0,4,5,3,7,8,6),"
				      "(0,3,6,1,4,7,2,5,8)}").parseIntegerVectorList();
     */

    ::SymmetryGroup s(n);
    s.computeClosure(L);

    FieldRationalFunctions F(Q,"t");
    FieldMatrix M(F,s.elements.size(),n);
    int I=0;
    for(::SymmetryGroup::ElementContainer::const_iterator i=s.elements.begin();i!=s.elements.end();i++,I++)
      {
	for(int j=0;j<n;j++)
	  {
	    M[I][j]=M[I][j]+F.exponent(j);
	    M[I][j]=M[I][j]-F.exponent(((::SymmetryGroup::inverse(*i))[j]));
	  }
      }
    AsciiPrinter P(Stderr);
    //    M.printMatrix(P);


    FieldMatrix M2(Q,M.getHeight(),M.getWidth());
    for(int i=0;i<M.getHeight();i++)
      for(int j=0;j<M.getWidth();j++)
	M2[i][j]=F.substitute(M[i][j],Q.zHomomorphism(10000).inverse());

        Field &myField=F;
	{
      M=M2;
      myField=Q;//HERE
      }
    IntegerVectorList extreme;

    M.printMatrix(P);


    I=0;
    FieldElement minusOne=myField.zHomomorphism(-1);//HERE
    for(::SymmetryGroup::ElementContainer::const_iterator i=s.elements.begin();i!=s.elements.end();i++,I++)
      //    for(int i=0;i<s.elements.size();i++)
      {
	fprintf(Stderr,"IIII:%i\n",I);
	M[I]=minusOne*M[I];
	FieldVector b(myField,s.elements.size());//HERE
	b[I]=minusOne;
	FieldLP lp(M,b);
	FieldLP lp2=lp.withNoLineality();

	AsciiPrinter P(Stderr);
	//lp2.print(P);

	//	cerr <<"Result" << lp2.findFeasibleBasis()<<endl;

	if(lp2.findFeasibleBasis())
	  extreme.push_back(*i);

	M[I]=minusOne*M[I];
      }
    fprintf(Stdout,"Extreme permutations (%i):\n",(int)extreme.size());
    P.printVectorList(extreme);
  }


  class TestCase
  {
  public:
	  string folder;
	  TestCase(string const &folder_):
		  folder(folder_)
		  {

		  }
	  void fail()
	  {
		  cerr<<"Test failed:"<<folder<<endl;
		  assert(0);
	  }
	  bool compare(string a, string b)
	  {
		  FILE *A=fopen(a.c_str(),"r");
		  FILE *B=fopen(b.c_str(),"r");
		  assert(A);
		  assert(B);
		  while((!feof(A))&&(!feof(B)))
		  {
			  if(fgetc(A)!=fgetc(B))return false;
		  }
		  if(feof(A)!=feof(B))return false;
		  return true;
	  }
	  bool fileExists(string name)
	  {
		  FILE *f=fopen(name.c_str(),"r");
		  if(f)fclose(f);
		  return f;
	  }
	  void replaceStrings(char *dest, char *s, const char* exe, const char *examplePath)
	  {
		  while(*s)
		  {
			  if(*s=='%')			  {
				  if(s[1]=='s')
				  {
					  strcpy(dest,exe);dest+=strlen(exe);
				  }
				  else if(s[1]=='p')
				  {
					  strcpy(dest,examplePath);dest+=strlen(examplePath);
				  }
				  else
					  assert(0);
				  s+=2;			  }
			  else
				  *(dest++)=*(s++);
		  }
		  dest[0]=0;
	  }
	  /*
	   * Returns true if test was successful.
	   * Returns false if test was performed for the first time.
	   * Returns false if test fails // old:Asserts if test fails
	   */
	  bool perform(const char *exe)
	  {
		  string fileName=folder+"/command";
		  FILE *f=fopen(fileName.c_str(),"r");
		  if(!f)
		  {
			  cerr<<"Could not open file:\""<<fileName<<"\""<<endl;
			  assert(f);
		  }
		  char command[4096];
		  char *temp=fgets(command,4095,f);
		  fclose(f);
		  assert(temp);
		  for(int i=0;i<4096 && command[i];i++)if(command[i]=='\n'){command[i]=0;}
		  char command2[4096];
		  string input=folder+"/input";
//		  sprintf(command2,command,exe,exe,exe,exe);
		  replaceStrings(command2,command,exe,folder.c_str());
		  bool outputExists=fileExists(folder+"/output");
		  string outputName=folder+"/output";
		  if(outputExists)
		  {
			  outputName=outputName+"New";
		  }
		  int err=0;
		  {
			  string t="rm -f "+outputName;
			  err|=system(t.c_str());
		  }
		  string command3="cat <"+input+"|"+string(command2)+">"+outputName;
		  cerr<<"Running command:\""<<command3<<"\""<<endl;
		  err|=system(command3.c_str());
		  //assert(err==0);
		  if(err)return false;
		  if(outputExists)
			  outputExists=compare(folder+"/output",folder+"/outputNew");
		  return outputExists;
	  }
  };

  list<string> subFolderNames()
  {
#define tempName "GfAnTeMpTeStS"
	  char command[256];
	  int err=system("rm " tempName);
	  err=0;//Having err!=0 above is probably not at mistake. Rather the file did not exist.
	  sprintf(command,"ls %s>" tempName ,testSuiteFolderOption.getValue());
	  err|=system(command);
	  assert(err==0);

	  list<string> ret;
	  FILE *f=fopen(tempName,"r");
	  assert(f);
	  char name[256];
	  while(fgets(name,255,f))
	  {
		  for(int i=0;i<255 && name[i];i++)if(name[i]=='\n'){name[i]=0;}
		  if(name[0]>='0' && name[0]<='9')ret.push_back(string(testSuiteFolderOption.getValue())+"/"+string(name));
	  }
	  fclose(f);
	  return ret;
  }

  int main()
  {
	  {

//		  gfan::Matrix<CircuitTableInteger>::readMatrix(cin,4);
	  }

	  if(developerTestOption.getValue()!="testsuite")
	  {
//    lpRationalFunctionTest();
//    testRationalFunctionField();
//	  packedTest();return 0;
//	  return testIntegers();
//	  return testGCD();
	  return testGfanLib();
//	  return testPolynomialGCD();
//	  return testLLL();
	  }
	  else
	  {
	  list<string> testFolders=subFolderNames();
	  list<TestCase> testList;
	  for(list<string>::const_iterator i=testFolders.begin();i!=testFolders.end();i++)
	  {
		  testList.push_back(TestCase(*i));
	  }

	  cout<<"Number of tests to perform "<<testList.size()<<endl;

	  list<string> failed;

	  int good=0;
	  int bad=0;
	  for(list<TestCase>::iterator i=testList.begin();i!=testList.end();i++)
		  if(i->perform(executableOption.getValue()))
			  good++;
		  else
		  {
			  bad++;
			  failed.push_back(i->folder);
		  }
	  cout<<"\n";
	  if(!failed.empty())
	  {
		  cout<<"Failed tests:\n-------------\n";
		  for(list<string>::iterator i=failed.begin();i!=failed.end();i++)
		  {
			  cout<<*i<<" FAILED!\n";
		  }
	  }

	  cout<<"Number of succesful tests "<<good<<endl;
	  cout<<"Number of failed tests "<<bad<<endl;

	  }
	  return 0;
  }
};

static TestApplication theApplication;
