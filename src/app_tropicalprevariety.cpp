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
#include "gfanlib_tropicalhomotopy.h"
#include "log.h"

using namespace gfan;

typedef gfan::CircuitTableInt32 mvtyp;

class TropicalPrevarietyApplication : public GFanApplication
{
	StringOption saveAs;
	StringOption loadFrom;
	IntegerOption optionNThreads;
	SimpleOption optionEarlyExit;
	IntegerOption optionBitsInIntegerType;
	SimpleOption optionUseValuation;
	SimpleOption optionMinT;
	SimpleOption optionMinX;

public:
  bool includeInDefaultInstallation()
  {
    return true;
  }
  const char *helpText()
  {
    return "This program computes the tropical prevariety of a set of polynomials. "
    		"The program works in two different modes. "
    		"By default the coefficient field is regarded as having trivial valuation and the output is a polyhedral fan. "
    		"With the option --usevaluation the program uses the valuation of the given coefficient field when defining the tropical prevariety. "
    		"As a consequence the output in this case is a polyhedral complex represented by the half-open polyhedral fan over it. \n"
    		"It is often useful to use the field Q(t). "
    		"However, because of the limitations of the current parser, "
    		"the way to specify elements of for example Q(t)[x1,x2] is to write them as polynomials in Q[t,x1,x2].\n\n"
    		"When --usevaluation is used the options --mint and --minx can switch the program from max to min convention.\n"
    		"In particular, using both --mint and --minx switches conventions to the [Speyer,Sturmfels]/[Maclagan,Sturmfels] convention.\n"
    		"Example:\n"
    		"./gfan _tropicalprevariety -j32 --usevaluation --bits64\n"
    		"with input:\n"
    		" Q(t)[x,y]{1+t*x+x+y}\n"
    		"produces a polyhedral complex."
    		;
  }
  TropicalPrevarietyApplication():
	  optionNThreads("-j","Number of threads",8),
	  optionEarlyExit("--earlyexit","Terminate after doing the intersection. Do not extract complex. This is useful for timing the code."),
	  saveAs("--saveas","Specify filename for saving intermediate result."),
	  loadFrom("--loadfrom","Specify filename for loading intermediate result."),
	  optionBitsInIntegerType("--bits","Number of bits used in intermediate integer type. Allowed values are 64, 32 and 0 for arbitrary precision. In the postprocessing step the computed prevariety is converted in arbitrary precission.",0),
	  optionUseValuation("--usevaluation","Use valuations in the defintion of the tropical hypersurfaces. The used valuation comes from the ring specified. For example Q(t) will come with its usual valuation. In order to represent the output an additional 0th coordinate is introduced."),
	  optionMinT("--mint","In a field with parameter t this option decides if maximal or minimal degrees of t determine the valuation."),
	  optionMinX("--minx","Switch the programme to min convention for the variables.")

  {
	  optionEarlyExit.hide();
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

template<class mvtyp, class mvtypDouble, class mvtypDivisor>
class SpecializedRTraverserStable: public SpecializedRTraverser<mvtyp,mvtypDouble,mvtypDivisor>
{
	std::vector<gfan::Matrix<mvtyp> > const tupleCopy; // This could be avoided, but I do not understand C++.
	gfan::Vector<mvtyp> const targetCopy;
public:
	SpecializedRTraverserStable(std::vector<gfan::Matrix<mvtyp> > const &tuple_, gfan::Vector<mvtyp> const target):
		SpecializedRTraverser<mvtyp,mvtypDouble,mvtypDivisor>(tuple_,&target),
		tupleCopy(tuple_),
		targetCopy(target)
	{
	}
	void process()
	{
		typedef SpecializedRTraverser<mvtyp,mvtypDouble,mvtypDivisor> P;
		std::cerr<<"CELL\n";
		gfan::Matrix<mvtyp> M(0,tupleCopy.size());
		gfan::Matrix<mvtyp> M2(P::T.traversers[P::T.level].choices.size(),1);
		int i=0;
		for(auto &e:P::T.traversers[P::T.level].choices)
		{
			std::cerr<<"CELL\n";
			M.appendRow(tupleCopy[i].column(e.first)-tupleCopy[i].column(e.second));
			std::cerr<<"CELL\n";
			M2[i][0]=(targetCopy)[e.first]-(targetCopy)[e.second];
			std::cerr<<"CELL\n";
			std::cerr<<e.first<<e.second<<"\n";
			i++;
		}
		auto N=combineLeftRight(M2,M);
		auto N2=gfan::Matrix<CircuitTableInteger>(N);
		std::cerr<<N2.toString();
		gfan::Cone<CircuitTableInteger> C(
				gfan::Matrix<CircuitTableInteger>::rowVectorMatrix(
				gfan::Vector<CircuitTableInteger>::standardVector(N2.getWidth(),0)).transposed(),
				N2.transposed());
		std::cerr<<C.getRays().toString();
//		std::cerr<<		N2.reduceAndComputeKernel().toString();
//		mixedVolume.addWithOverflowCheck(SpecializedRTraverser<mvtyp,mvtypDouble,mvtypDivisor>::T.traversers[SpecializedRTraverser<mvtyp,mvtypDouble,mvtypDivisor>::T.level].inequalityTable.getVolume().extend());
	}
};

	template<typename typ>
	vector<HalfOpenCone<typ> > mainPart1(bool useValuation, bool maxt, bool switchtsign, bool minconvention=false) //ideally maxt should be known to the field
		{

	    FileParser S(Stdin); // allows to read and parse polynomials, etc.
	//StringParser S("Q[t,x,y]{t+x+y,t+x+y}");
	    const PolynomialSet g = S.parsePolynomialSetWithRing(); // system of polynomials
	    int n=g.getRing().getNumberOfVariables();
	    log1 std::cerr<<"DONE reading input\n";

	          vector<gfan::Matrix<typ> > configurations;
		  for(auto &p:g)
		  {
			  if(useValuation)
			  {
				  IntegerVectorList temp;
				  for(auto &t:p.terms)
					  temp.push_back(concatenation(IntegerVector::singleEntry(t.second.valuation(maxt)*((switchtsign)?-1:1)),t.first.exponent));
				  configurations.push_back(convertMatrix<typ>(rowsToIntegerMatrix2(temp)));
			  }
			  else
			  {
				  IntegerVectorList temp;
				  for(auto &t:p.terms)
					  temp.push_back(t.first.exponent);
				  configurations.push_back(convertMatrix<typ>(rowsToIntegerMatrix2(temp)));
//				  configurations.push_back(convertMatrix<typ>(rowsToIntegerMatrix2(p.exponents())));
			  }
		  }

		  if(minconvention)
			  for(auto &A:configurations)A=-A;

#if 0
		  if(0)
		  {// Stable intersection -- assuming square system
	          vector<gfan::Matrix<mvtyp> > configurations;
		  for(auto &p:g)
			  configurations.push_back(convertMatrix<mvtyp>(rowsToIntegerMatrix2(p.exponents())));

			  vector<gfan::Matrix<mvtyp> > configurations2;
			  gfan::Vector<mvtyp> target;
			  int newD=configurations.size();
			  for(auto &A:configurations)
			  {
				  target=concatenation(target,A[newD].toVector());
				  configurations2.push_back(A.submatrix(0,0,newD,A.getWidth()));
			  }
			  SpecializedRTraverserStable<mvtyp,mvtyp::Double,mvtyp::Divisor> traverser(configurations2,target);
			  traverse_simple(&traverser);
			  return 0;
		  }
#endif


		  vector<HalfOpenCone<typ> > f;

		  //We pretend that first coordinate is special (that it is the t parameter)
//		  int n=configurations[0].getWidth();
		  int nonEmptyIntersections=0;
	//	  vector<int> edgeCounts;
	//	  for(int i=0;i<configurations.size();i++)edgeCounts.push_back(numberOfEdges(configurations[i]));
		  vector<int> used(configurations.size());
		  PolytopeIntersectionData<typ> data(configurations);
		  log1 debug<<"Check the line above. Does it try to figure out intersection data for a part of space?\n";
		  /*
		   * The answer is no. The PolytopeData does not know in which part of the space we wish to do the
		   * intersection. There may be room for improvement. Notice however, that PolytopeData's relation
		   * table are only used for _guiding_ the dynamic enumeration.
		   */

		  gfan::Matrix<typ> strictInequalities(n+useValuation,useValuation);
		  if(useValuation)strictInequalities[0][0]=typ(+1);//--------

		  ProgressCounter progress;
		  CommonStatistics statistics;

		  if(loadFrom.getValue())
		  {
			  string loadName(loadFrom.getValue());
			  log1 std::cerr<<"Loading vector of HalfOpenCones from file \""<<loadName<<"\"\n";
			  std::ifstream F;
			  F.open(loadName.c_str(),std::ifstream::in);
			  f=loadHalfOpenConeVector<typ>(F);
			  F.close();
		  }
		  else
		  {
			  log1 std::cerr<<"ABOUT TO INTERSECT\n";
		  commonRefinement(
				  gfan::HalfOpenCone<typ>(
						  gfan::Matrix<typ>(n+useValuation,0),
						  gfan::Matrix<typ>(n+useValuation,0),
						  strictInequalities),
				  used,configurations,f,nonEmptyIntersections/*,&edgeCounts*/,data,RelationTable(data.layout),progress,&statistics,optionNThreads.getValue()/*8*/);
		  }
		  if(saveAs.getValue())
		  {
			  string saveName(saveAs.getValue());
			  log1 std::cerr<<"Saving vector of HalfOpenCones in file \""<<saveName<<"\"\n";
			  std::ofstream F;
			  F.open(saveName.c_str(),std::ofstream::out);
			  save(f,F);
			  F.close();
		  }

		  log1 std::cerr<<"Number of halfopen cones in output:"<<f.size()<<"\n";
//		  std::cerr<<"NonEmptyIntersections:"<<statistics.numberOfIntermediateVertices<<"\n";

		  assert(!f.empty());

		  log1 std::cerr<<"INTERSECTION COMPLETED\n";
		  return f;
		}
	template<typename typ>
	vector<HalfOpenCone<CircuitTableInteger> > mainPart2(vector<HalfOpenCone<typ> > &f)
	{
    vector<HalfOpenCone<CircuitTableInteger> > f2;
    f2.reserve(f.size());
    for(auto &a:f)
  	  if(!a.isEmpty())
  	  {
  		  auto temp=HalfOpenCone<CircuitTableInteger>(a);
  		  f2.push_back(temp);
  	  }
		  else
		  {
			  std::cerr<<"Empty set computed\n"<<a.toString()<<"\nLifted:\n"<<a.lifted.toString()<<"\n";
		  }

    return f2;
 }

	int main()
  {
		vector<HalfOpenCone<CircuitTableInteger> >f2;
		vector<gfan::Matrix<CircuitTableInteger> > f2Rays;

	  //	typedef gfan::CircuitTableInt128 typ; // type for integers
//    typedef gfan::CircuitTableInt64 typ; // type for integers

		bool maxt=!optionMinT.getValue();
		bool minconvention=optionMinX.getValue();

		bool switchtsign=(!maxt)^minconvention;

		if(optionBitsInIntegerType.getValue()!=0)
			cerr<<"Warning: Using machine integer arithmetic. Overflows may not all be caught!\n";
		if(optionBitsInIntegerType.getValue()==32)
    {
    	using typ=gfan::CircuitTableInt32;
    	auto f=mainPart1<typ>(optionUseValuation.getValue(),maxt,switchtsign,minconvention);
    	if(optionEarlyExit.getValue())return 0;
    	f2=mainPart2<typ>(f);
    	f2Rays.reserve(f.size()); //strange size!
    	for(auto &a:f2)f2Rays.push_back(a.closure().getRays());
    }
    else if(optionBitsInIntegerType.getValue()==64)
    {
    	using typ=gfan::CircuitTableInt32;
    	auto f=mainPart1<typ>(optionUseValuation.getValue(),maxt,switchtsign,minconvention);
    	if(optionEarlyExit.getValue())return 0;
    	f2=mainPart2<typ>(f);
    	f2Rays.reserve(f.size()); //strange size!
    	for(auto &a:f2)f2Rays.push_back(a.closure().getRays());
    }
    else if(optionBitsInIntegerType.getValue()==0)
    {
    	using typ=gfan::CircuitTableInteger;
    	auto f=mainPart1<typ>(optionUseValuation.getValue(),maxt,switchtsign,minconvention);
    	if(optionEarlyExit.getValue())return 0;
    	f2=mainPart2<typ>(f);
    	f2Rays.reserve(f.size()); //strange size!
    	for(auto &a:f2)f2Rays.push_back(a.closure().getRays());
    }


    if(optionBitsInIntegerType.getValue()!=0)
    	std::cerr<<"DONE CONVERTING FROM MACHINE INTEGERS\n";
    log1 std::cerr<<"DONE COMPUTING RAYS\n";

 	  auto linealitySpace=f2.begin()->closure().getLinealitySpace();
 	  RayCollector<CircuitTableInteger> collector(linealitySpace);

 	  for(int i=0;i<f2Rays.size();i++)
 		  {
 			  auto rays=f2Rays[i];
 			  for(int i=0;i<rays.getHeight();i++)
 			  {
 				  int index=collector.lookup(normalize(rays[i].toVector()));
 			  }
 		  }
 #if 1

 //				std::cerr<<"THE RAYS\n"<<collector.getRays()<<"\n";
 	  log1 std::cerr<<"N_RAYS\n"<<collector.getRays().getHeight()<<"\n";
 /*				for(auto &c:f)
 		{
 			std::cerr<<"HALFOPENCONE:\n";
 			c.extractFaceComplex();
 		}*/
 	  std::cout<<extractFaceComplex(f2,collector,linealitySpace,12,0,0,logLevel).toString(2+4+8+2048);

 //	  std::cerr<<fvector(f).toString()<<"\n";
 #endif

 //	  std::cerr<<"Number of computed half open cones:"<<f2.size()<<"\n";
 //	  std::cerr<<"Number of non-empty half open cones:"<<numberOfNonEmptyHalfOpenCones<<"\n";
 //	  std::cerr<<"NonEmptyIntersections:"<<statistics.numberOfIntermediateVertices<<"\n";
	  return 0;
  }
};

static TropicalPrevarietyApplication theApplication;

