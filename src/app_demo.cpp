/*
 * app_demo.cpp
 *
 *  Created on: Jan 18, 2022
 *      Author: anders
 */




#include <iostream>
#include <fstream>
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

class TropicalPrevariety : public GFanApplication
{
	StringOption saveAs;
	StringOption loadFrom;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the tropical prevariety of a set of polynomials given in a polynomialring with first variable being regarded as a field element with valuation 1.\n";
  }
  TropicalPrevarietyComponentsApplication():
	  saveAs("--saveas","Specify filename for saving intermediate result."),
	  loadFrom("--loadfrom","Specify filename for loading intermediate result.")
  {
	  registerOptions();
  }

  const char *name()
  {
    return "_tropicalprevariety";
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

template<typename typ>
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
}


  int main()
  {
    //	typedef gfan::CircuitTableInt128 typ; // type for integers
    typedef gfan::CircuitTableInt64 typ; // type for integers


    FileParser S(Stdin); // allows to read and parse polynomials, etc.
    const PolynomialSet g = S.parsePolynomialSetWithRing(); // system of polynomials
    std::cerr<<"DONE reading input\n";

          vector<gfan::Matrix<typ> > configurations;
	  for(auto &p:g)
		  configurations.push_back(convertMatrix<typ>(rowsToIntegerMatrix2(p.exponents())));

	  vector<HalfOpenCone<typ> > f;

	  //We pretend that first coordinate is special (that it is the t parameter)
	  int n=configurations[0].getWidth();
	  int nonEmptyIntersections=0;
	  vector<int> edgeCounts;
	  for(int i=0;i<configurations.size();i++)edgeCounts.push_back(numberOfEdges(configurations[i]));
	  vector<int> used(configurations.size());
	  PolytopeIntersectionData<typ> data(configurations);

	  gfan::Matrix<typ> strictInequalities(n,1);
	  strictInequalities[0][0]=typ(-1);

	  ProgressCounter progress;
	  CommonStatistics statistics;

	  if(loadFrom.getValue())
	  {
		  string loadName(loadFrom.getValue());
		  std::cerr<<"Loading vector of HalfOpenCones from file \""<<loadName<<"\"\n";
		  std::ifstream F;
		  F.open(loadName.c_str(),std::ifstream::in);
		  f=loadHalfOpenConeVector<typ>(F);
		  F.close();
	  }
	  else
	  {
	  std::cerr<<"ABOUT TO INTERSECT\n";
	  commonRefinement(
			  gfan::HalfOpenCone<typ>(
					  gfan::Matrix<typ>(n,0),
					  gfan::Matrix<typ>(n,0),
					  strictInequalities),
			  used,configurations,f,nonEmptyIntersections,&edgeCounts,data,RelationTable(data.layout),progress,&statistics,32/*8*/);
	  }
	  if(saveAs.getValue())
	  {
		  string saveName(saveAs.getValue());
		  std::cerr<<"Saving vector of HalfOpenCones in file \""<<saveName<<"\"\n";
		  std::ofstream F;
		  F.open(saveName.c_str(),std::ofstream::out);
		  save(f,F);
		  F.close();
	  }

	  assert(!f.empty());

	  std::cerr<<"INTERSECTION COMPLETED\n";


      vector<HalfOpenCone<CircuitTableInteger> > f2;
      vector<gfan::Matrix<CircuitTableInteger> > f2Rays;
      f2.reserve(f.size());
      f2Rays.reserve(f.size());
      for(auto &a:f)
    	  if(!a.isEmpty())
    	  {
    		  auto temp=HalfOpenCone<CircuitTableInteger>(a);
    		  f2Rays.push_back(temp.closure().getRays());
    		  f2.push_back(temp);
    	  }
		  else
		  {
			  std::cerr<<"Empty set computed\n"<<a.toString()<<"\nLifted:\n"<<a.lifted.toString()<<"\n";
		  }

      std::cerr<<"DONE CONVERTING FROM INT64 AND COMPUTING RAYS\n";

	  auto linealitySpace=f2.begin()->closure().getLinealitySpace();
	  RayCollector<CircuitTableInteger> collector(linealitySpace);

#if 0

//				std::cerr<<"THE RAYS\n"<<collector.getRays()<<"\n";
	  std::cerr<<"N_RAYS\n"<<collector.getRays().getHeight()<<"\n";
/*				for(auto &c:f)
		{
			std::cerr<<"HALFOPENCONE:\n";
			c.extractFaceComplex();
		}*/
///*				std::cerr<<*/ extractFaceComplex(f,collector,linealitySpace);

	  std::cerr<<fvector(f).toString()<<"\n";
#endif

//	  std::cerr<<"Number of computed half open cones:"<<f2.size()<<"\n";
//	  std::cerr<<"Number of non-empty half open cones:"<<numberOfNonEmptyHalfOpenCones<<"\n";
//	  std::cerr<<"NonEmptyIntersections:"<<statistics.numberOfIntermediateVertices<<"\n";

	  return 0;
  }
};

static TropicalPrevarietyApplication theApplication;

