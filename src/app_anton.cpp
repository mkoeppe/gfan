/*
 * app_anton.cpp
 *
 *  Created on: 12 Mar 2020
 *      Authors: anders, anton
 */

#include "app_anton.h"

using namespace gfan;

class AntonApplication : public GFanApplication
{
	  SimpleOption fVectorOption;
	  SimpleOption filterOption;
	  SimpleOption doubleEdgeOption;
	  IntegerOption nThreadsOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This application reads in a polynomial ring and a set of polynomials and computes the intersection of the tropical hypersurfaces of the polynomials. The first variable is treated specially. Namely, it is the parameter of the field and has valuation 1.\n";
  }
  AntonApplication():
	fVectorOption("--f-vector","This option enables the computation of f-vector."),
	filterOption("--filter","This option enables the filter."),
	doubleEdgeOption("--double-edge","This option uses a specialization of valuations assuming vertices 1 and 2 in Albouy-Kaloshin diagrams are connected with a double edge."),
	nThreadsOption("-j","Sets the number of cpu threads to be used.",8)
  {
    registerOptions();
  }

  const char *name()
  {
    return "_anton";
  }


  // A few methods that should go somewhere else:
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

/*  template<typename typ>
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
  template<typename typ>
  static bool testIfUnboundedInCoordinateC(HalfOpenCone<typ> &c, int specialCoordinate=0, int cCoordinate=0)
  {
    auto closed = c.closure();
    int n = closed.getAmbientDimension();
    //auto hyperplane_t_0=Cone<typ>::hyperPlane(Vector<typ>::standardVector(n,specialCoordinate));
    auto recessionCone=c.representedPolyhedronsRecessionConeInHigherDim();
    auto hyperplane_c_0=Cone<typ>::hyperPlane(Vector<typ>::standardVector(n,cCoordinate));
    auto temp=intersection(recessionCone,hyperplane_c_0);
    return temp.getDimension()>0;
  }
  
  template<typename typ>
  static bool testRecessionConeAgainstAllOnesVector(HalfOpenCone<typ> &c)
  {
    //if(c.isEmpty())return true;
    auto recessionCone=c.representedPolyhedronsRecessionConeInHigherDim();
    int n=recessionCone.getAmbientDimension()-1;
    auto H2=Cone<typ>::halfSpace(Vector<typ>::allOnes(n+1));
    return intersection(H2,recessionCone).getDimension()>0;
  }
  template<typename mvtyp>
  class L: public std::function<bool(gfan::HalfOpenCone<mvtyp>&)>
		{
	  	  bool &noFilter;
	  	  int const &numberOfVariables;
		public:
	  	  L(bool &noFilter_,int const &numberOfVariables_):
	  		  noFilter(noFilter_),
			  numberOfVariables(numberOfVariables_)
		{
		}
	  	  bool oprator(gfan::HalfOpenCone<mvtyp>& c)
	  	  {
	  	      return noFilter or // to enable: ./gfan _anton --filter
	  	      // keep cones that...
	  	      ( testIfUnboundedInCoordinateC(c,0,numberOfVariables-1/*c coordinate*/)  // ... have infinite fibers for c-projection
	  	        and testRecessionConeAgainstAllOnesVector(c)
	  	        ); // ... and are "positive"
	  	  }
		};

  int main() 
  {
//	  typedef gfan::CircuitTableInteger typ;
    typedef gfan::CircuitTableInt64 typ; // type for integers
    FileParser S(Stdin); // allows to read and parse polynomials, etc.
    const PolynomialSet g = S.parsePolynomialSetWithRing(); // system of polynomials
    const int numberOfVariables = g.getRing().getNumberOfVariables();
    vector< gfan::Matrix<typ> > configurations; // vectors of exponents (of the monomials in polynomials)
    for(auto &p : g)
      configurations.push_back(convertMatrix<typ>(rowsToIntegerMatrix2(p.exponents())));

    vector< HalfOpenCone<typ> > remainingCones; // empty (in the beginning) 

    int n = configurations[0].getWidth(); // dimension of the ambient space
    int m = configurations.size(); // # of input polynomials
//    vector<int> edgeCounts; // # edges in Newton polytope of input polynomials
//    for(int i=0; i<m; i++)
//      edgeCounts.push_back(numberOfEdges(configurations[i]));
    gfan::Matrix<typ> strictInequalities(n,1);
    strictInequalities[0][0]=typ(-1); // -(first variable) > 0
    auto nonstrictInequalities = gfan::Matrix<typ>(n,0);
    if (doubleEdgeOption.getValue()) { 
      // hardcoded for examples/n-5--p5--IIUc-BIG-not-normalized
      StringParser S("{"
		     "t^-25*a12^4-z1-t^-125*a13*b13^3-t^-625*a14*b14^3-t^-3125*a15*b15^3,"
		     "t^-5*a12^4-z2+t^-125*a23*b23^3-t^-625*a24*b24^3-t^-3125*a25*b25^3,"
		     "t^-25*a12^4-w1-t^-125*a13^3*b13-t^-625*a14^3*b14-t^-3125*a15^3*b15,"
		     "t^-5*a12^4-w2+t^-125*a23^3*b23-t^-625*a24^3*b24-t^-3125*a25^3*b25,"
		     "t^-15625*a12^4 + a34^3*b34 + b34^3*a23 + a35^3*b35 + b35^3*a35 + a45^3*b45 + b45^3*b45 + z3+z4+z5 + w3+w4+w5"  
		     "}");
      const PolynomialSet G = S.parsePolynomialSet(g.getRing()); // system of polynomials
      nonstrictInequalities = convertMatrix<typ>(rowsToIntegerMatrix2(wallInequalities(G))).transposed();
    }
    const auto HOC = gfan::HalfOpenCone<typ>(nonstrictInequalities, // inequalities (...>=0)
					     gfan::Matrix<typ>(n,0), // equations (...=0)
					     strictInequalities);    // strict inequalities (...>0)
    bool noFilter = not filterOption.getValue(); // custom filter parameter used in the lambda-function below 
    const auto filterFunction = [noFilter,numberOfVariables](HalfOpenCone<typ>& c) -> bool {
      return noFilter or // to enable: ./gfan _anton --filter
      // keep cones that...
      ( testIfUnboundedInCoordinateC(c,0,numberOfVariables-1/*c coordinate*/)  // ... have infinite fibers for c-projection 
        // and testRecessionConeAgainstAllOnesVector(c)
        ); // ... and are "positive"
    }; 
    // THINGS THAT ARE NOT USED (beyond commonRefinement call)
    CommonStatistics statistics;
    int nonEmptyIntersections = 0; // some statistic...
    vector<int> used(m);   
    PolytopeIntersectionData<typ> data(configurations);
    ProgressCounter progress;
    std::function<bool (HalfOpenCone<typ>&)> filterF2=filterFunction; // check if this is necessary!!!
    // (end) THINGS...
    commonRefinement(HOC,used,configurations,remainingCones,
		     nonEmptyIntersections,/*&edgeCounts,*/data,
		     RelationTable(data.layout),progress,
		     &statistics,nThreadsOption.getValue(),filterF2);
    if (fVectorOption.getValue()) 
      std::cerr<<"FVECTOR:"<<fvector(remainingCones)<<"\n";
    std::cerr<<"remainingCones.size():"<<remainingCones.size()<<"\n"<<std::flush;
 
    if (remainingCones.empty()) std::cerr<<"No remaining cones!!!\n";
    else {
      auto linealitySpace=remainingCones.begin()->closure().getLinealitySpace();
      std::cerr<<"LinealitySpace:\n"<<linealitySpace.toString()<<"\n";
      RayCollector<typ> collector(linealitySpace);
      int count = 0;
      for(auto &c:remainingCones) {
	if (count++ % 1000 == 0) 
	  std::cerr << "ray collector: " << remainingCones.size()-count << " cones remaining\n"; 
	if(!c.isEmpty())
	  {
	    auto rays=c.closure().getRays();
	    //              std::cerr<<"rays:"<<rays.toString();
	    for(int i=0;i<rays.getHeight();i++)
	      collector.lookup(normalize(rays[i].toVector()));
	  }
      }
      std::cerr << "N_RAYS\n"<<collector.getRays().getHeight() <<"\n";
      //std::cerr<<extractFaceComplex(remainingCones,collector,linealitySpace,0); // (not a complex) 
      //auto fc = extractFaceComplex(remainingCones,collector,linealitySpace,0,true/*closure*/,true/*only max cones*/);
      //PostComplex<typ> postcomplex(fc,true); 
      PostComplex<typ> postcomplex(remainingCones,collector);
      //std::cerr << postcomplex; 
      
      std::fstream f;
      f.open("STASHED_CONES.out",std::fstream::out);
      postcomplex.serialize(f);
      f.close();
      //postcomplex.printComponents();
    }
    return 0;
  } // main
}; // class AntonApplication

static AntonApplication theApplication;
