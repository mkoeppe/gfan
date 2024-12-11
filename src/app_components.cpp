/*
 * app_components.cpp
 *
 *  Created on: 12 Mar 2020
 *      Authors: anders, anton
 */

#include <thread>
#include <execution>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "app_anton.h"
//#include <algorithm>

using namespace gfan;

template<class mvtyp>
bool isSliceGood(const PostComplex<mvtyp>& pc, Rational p, 
		 const vector<Rational>& minConePotential,
		 const vector<Rational>& maxConePotential,
		 const bool saveNoncomets
		 )
{
  std::cerr << "--STARTED potential = " << p << std::endl;   
  auto nd = p.numeratorDenominator();
  mvtyp n(nd.first.toInt()), d(nd.second.toInt());
  PostComplex<mvtyp> slicedComplex(pc.dimension); // the new complex stores the slice of the old one
  RayCollector<mvtyp> rayCollector(gfan::Matrix<mvtyp>(0,pc.dimension));
  for(int count = pc.cones.size()-1; count>=0; count--) {
    if (count % 1000 == 0)
      std::cerr << "[ " << p << ", threadID=" << std::this_thread::get_id() << "] -- it remains " << count << " cones to slice\n"; 
    auto c = pc.cones[count];
    if(not (p<minConePotential[count]) and not (maxConePotential[count]<p)) {
      // this is an ad hoc implementation of "double-description method"
      PostCone intersectionLabels; // labels in slicedComplex
      std::vector<RayLabel> negativeLabels; // labels in pc
      std::vector<RayLabel> positiveLabels; // labels in pc	  
      for(auto r : c) {
	auto ray = pc.rays[r];
	auto diff = d * ray[pc.dimension-1] - n * ray[0];
	if(diff.isPositive())
	  positiveLabels.push_back(r);
	else if(diff.isNegative())
	  negativeLabels.push_back(r); 
	else {
	  auto newRay = ray.toVector();
	  intersectionLabels.insert(rayCollector.lookup(newRay));
	}
      }
      for(auto a : positiveLabels)
	for(auto b : negativeLabels) {
	  auto A = pc.rays[a];
	  auto B = pc.rays[b];
	  auto Ad = A[0], An = A[pc.dimension-1];
	  auto Bd = B[0], Bn = B[pc.dimension-1];
	  auto tA = - d * Bn + n * Bd;
	  auto tB = d * An - n * Ad;
	  assert(tA.isPositive());
	  assert(tB.isPositive());
	  auto newRay = normalize(tA*A.toVector()+tB*B.toVector());
	  intersectionLabels.insert(rayCollector.lookup(newRay));
	}
      // sanity check: at least one ray is finite!
      if (not intersectionLabels.empty()) {
	std::cerr << "[" << intersectionLabels.size() << "]";
	if (std::all_of(intersectionLabels.begin(),intersectionLabels.end(),
			[&](RayLabel r){return rayCollector.rays[r][0].isZero();}
			))
	  {
	    std::cerr << "error: no finite ray in this intersection!" << endl;
	    for (auto r : intersectionLabels)
	      std::cerr << rayCollector.rays[r].toVector().toString() << std::endl;
	  }
	else {
	  slicedComplex.cones.push_back(intersectionLabels);
	  if(not (p<minConePotential[count]) and not (maxConePotential[count]<p)) { // always true
	  } else { 
	    std::cerr << p << ": " << not (p<minConePotential[count]) << " and " << not (maxConePotential[count]<p) << " = " << (not (p<minConePotential[count]) and not (maxConePotential[count]<p)) << std::endl; 
	    std::cerr << "-- min = " << minConePotential[count] << ", max = " << maxConePotential[count] << " for cone:\n";
	    for (auto r : c)
	      std::cerr << pc.rays[r].toVector().toString() << std::endl;
	    std::cerr << "  intersection cone:\n";
	    for (auto r : intersectionLabels)
	      std::cerr << rayCollector.rays[r].toVector().toString() << std::endl; 
	  }
	}
      }
    }//if(...)
  }//for(count...)
  slicedComplex.rays = rayCollector.rays;
  slicedComplex.printComponents(slicedComplex.connectedComponents(),saveNoncomets);
  std::cerr << "--FINISHED potential = " << p << std::endl;   
  return true;//!!!
}

template<class mvtyp>
void computePotentials(
				  PostComplex<mvtyp>& pc, 
				  std::set<Rational>& potentials,
				  Rational& minPotential, 
				  Rational& maxPotential,
				  std::vector<Rational>& minConePotential, 
				  std::vector<Rational>& maxConePotential
				  )
{
  for(int r=0; r < pc.nRays(); r++)  
    if(not pc.rays[r][0].isZero()) {
      Rational lastCoord = pc.rays[r][pc.dimension-1].toInt64();
      Rational potential = lastCoord / pc.rays[r][0].toInt64();
      if (r==0) minPotential = maxPotential = potential;
      if (potential < minPotential) 
	minPotential = potential;
      else if (maxPotential < potential)
	maxPotential = potential;
      potentials.insert(potential);
    }
  for(const auto& c : pc.cones) {
    Rational min = maxPotential+1, max = minPotential-1;
    for(RayLabel r : c) {
      auto firstCoord = pc.rays[r][0];
      Rational lastCoord = pc.rays[r][pc.dimension-1].toInt64();
      if(firstCoord.isZero()) {
	if(lastCoord.sign()<0) // here the assumption that firstCoord < 0 is crucial!!!
	  max = maxPotential;
	else if(lastCoord.sign()>0)
	  min = minPotential;
      } else {
	Rational potential = lastCoord / firstCoord.toInt64();
	if (potential < min) 
	  min = potential; 
	if (max < potential)
	  max = potential;
      }
    }
    minConePotential.push_back(min);
    maxConePotential.push_back(max);
  }
}

template<class mvtyp>
void connectedComponentsOfSlices(PostComplex<mvtyp>& pc, bool saveNoncomets, bool potentialsOutputAndStop)
{
  std::set<Rational> potentials;
  Rational minPotential, maxPotential;
  std::vector<Rational> minConePotential, maxConePotential;
  computePotentials(pc,potentials,minPotential,maxPotential,minConePotential,maxConePotential);
  std::cerr << "potentials: " << potentials << std::endl;

  /*std::vector<std::thread> workers;
  for (auto& p : potentials) 
  workers.push_back(std::thread(
    [&]() //begin lambda
    {
      isSliceGood(pc,p,minConePotential,maxConePotential,saveNoncomets);
    }
    ));  
  std::for_each(workers.begin(), workers.end(), [](std::thread &t) 
		{
		  t.join();
		});
  */
  if (potentialsOutputAndStop) {
    std::for_each(
		  std::execution::seq, 
		  potentials.begin(),
		  potentials.end(),
		  [](Rational p)
		  {
		    auto nd = p.numeratorDenominator();
		    int n(nd.first.toInt()), d(nd.second.toInt());
		    std::filesystem::create_directory("POTENTIALS");
		    std::stringstream filename;		    
		    filename << "POTENTIALS/" << n << "over" << d;
		    std::ofstream of(filename.str());
		    if (of.fail())
		      std::cerr << "can't open " << filename.str() << std::endl;
		    of << n << std::endl << d << std::endl;
		    of.close();
		  });
  } else {
    std::for_each(
		  std::execution::par,//fails to trigger parallelism
		  potentials.begin(),
		  potentials.end(),
		  [pc,minConePotential,maxConePotential,saveNoncomets](Rational p)
		  {
		    isSliceGood(pc,p,minConePotential,maxConePotential,saveNoncomets); 
		  });
  }
}

class ComponentsApplication : public GFanApplication
{
  SimpleOption saveMaximalCones;
  SimpleOption potentialSliceOption;
  IntegerOption potentialNumerator;
  IntegerOption potentialDenominator;
  SimpleOption potentialsOutputAndStop;
  SimpleOption saveNoncomets;
  IntegerOption nThreadsOption;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This application reads in a polyhedral complex stored as a collection of its maximal cones. Each cone is generated by a subset of rays stored as vectors of integers. The first coordinate is treated specially: it is the parameter t of the field and has valuation 1. It computes the set of connected components of the complex obtained by intersecting the input complex with the hyperplane t=1.\n";
  }
  ComponentsApplication():
    saveMaximalCones("--save-maximal-cones","Save the (post)complex containing only maximal cones ."),
    potentialSliceOption("--potential-slice","Invoke slicing of the complex via fixed-potential hyperplanes."),
    potentialNumerator("--potential-numerator","Fixed potential numerator (used with --potential-slice)."),
    potentialDenominator("--potential-denominator","Fixed potential denominator (used with --potential-slice)."),
    potentialsOutputAndStop("--potential-output-and-stop","If --potential-slice, then output the potentials of slices and stop."),
    saveNoncomets("--save-noncomets","Create *.out files for all components with are not comets."),
    nThreadsOption("-j","Set the number of cpu threads to be used.",8)
  {
    registerOptions();
  }

  const char *name()
  {
    return "_components";
  }


  int main() 
  {
    using INTEGER = CircuitTableInteger;
    //using INTEGER = CircuitTableInt64;
    PostComplex<INTEGER> postcomplex(std::cin,saveMaximalCones.getValue());   
    if(saveMaximalCones.getValue()) {
      std::fstream f;
      f.open("MAXIMAL_CONES.out",std::fstream::out);
      postcomplex.serialize(f);
      f.close();
    }
    if (potentialSliceOption.getValue()) {
      int d = potentialDenominator.getValue();
      int n = potentialNumerator.getValue();
      if (d==0)  
	connectedComponentsOfSlices(postcomplex,saveNoncomets.getValue(),potentialsOutputAndStop.getValue());
      else { 
	std::set<Rational> potentials;
	Rational minPotential, maxPotential;
	std::vector<Rational> minConePotential, maxConePotential;
	computePotentials(postcomplex,potentials,minPotential,maxPotential,minConePotential,maxConePotential); // the above can be stored or loaded once and reused 
	isSliceGood(postcomplex,Rational(n)/Rational(d),minConePotential,maxConePotential,saveNoncomets.getValue()); 
      }
    }
    else {
      auto components = postcomplex.connectedComponents();
      postcomplex.printComponents(components,saveNoncomets.getValue()); 
    }
    return 0;
  } // main
}; // class ComponentsApplication

static ComponentsApplication theApplication;
