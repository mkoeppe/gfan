#ifndef GFANLIB_app_anton_H_
#define GFANLIB_app_anton_H_

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include "parser.h"
#include "gfanapplication.h"
#include "linalg.h"
#include "gfanlib_tableau.h"
#include "gfanlib_hypersurfaceintersection.h"
#include "gfanlib_circuittableinteger.h"
#include "gfanlib_q.h"
#include "wallideal.h"
//#include "rational.h"

template <class Os, class K>
  Os& operator<<(Os& os, const std::set<K>& v) {
  os << '[' << v.size() << "] {";
  bool o{};
  for (const auto& e : v)
    os << (o ? ", " : (o = 1, " ")) << e;
  return os << " }\n";
}


namespace gfan{
  
using RayLabel = int;
using ConeLabel = int;
struct Component : set<RayLabel> {
  int number;
};
using PostCone = std::set<RayLabel>;

template<class mvtyp> 
class PostComplex{
  Matrix<mvtyp> selectRays(const PostCone& c) const {
    Matrix<mvtyp> M(0,rays.getWidth());
    for (auto r : c)
      M.appendRow(rays[r]);
    return M;
  }

public:  
  int dimension;
  Matrix<mvtyp> rays; // rays are rows in this matrix 
  std::vector<PostCone> cones; // maximal cones

  // CONTRUCTORS		  
 PostComplex(vector<HalfOpenCone<mvtyp>> remainingCones, 
	     RayCollector<mvtyp> collector)  
   : dimension(collector.getRays().getWidth()),
     rays(collector.getRays())
  {
    for(auto& c : remainingCones) {
      //if(complex.isMaximal(c)) 
      auto r = c.closure().getRays();
      std::set<RayLabel> cone;
      for(int i=0;i<r.getHeight();i++) {
	RayLabel index=collector.lookup(normalize(r[i].toVector()));
	cone.insert(index);
      }
      cones.push_back(cone);
    }
  }

  PostComplex(const int d)  
   : dimension(d),
     rays(Matrix<mvtyp>(0,d))
  {
  }
  
  PostComplex(const Complex<mvtyp>& complex)  
   : dimension(complex.n),
     rays(complex.vertices)
  {
    for(auto& c : complex.cones) 
      cones.push_back(std::set<RayLabel>(c.indices.begin(),c.indices.end()));
  }
  
  bool includes(const PostCone& a, const PostCone& b)
  {
    return std::includes(a.begin(),a.end(),b.begin(),b.end());
  }

  PostComplex(istream& is, bool filterMaximal = true)
  {
    std::string line;
    std::string lead; 
  
    is >> lead; //"DIMENSION"
    is >> dimension;
    std::cerr << "dimension = " << dimension << std::endl;   

    is >> lead; //"RAYS"
    rays = gfan::Matrix<mvtyp>::readMatrix(is,dimension);
    std::cerr << "#rays = " << nRays() << std::endl;   

    std::vector<bool> isConeMaximal;  
    std::vector<std::set<ConeLabel>> maxConesContainingRay(nRays()); 

    while(std::getline(is,line)) { // this assumes that the file ends with after "cones"
      std::istringstream iss(line);
      iss >> lead; //"CONE"
      if (lead != "CONE") {
	std::cerr << "--skipping line\n" << /*line << */ std::endl;
	continue;
      }
      PostCone c;
      RayLabel r;
      while(iss >> r) 
	c.insert(r);
      cones.push_back(c);
      if (cones.size()%1000==0)
	std::cerr << cones.size() << " cones processed" << std::endl;
      bool isMaximal = true;
      if (filterMaximal) {
	set<ConeLabel> smallerCandidates;
	for(auto rL : c)
	  smallerCandidates.merge(maxConesContainingRay[rL]);
	for(auto cL : smallerCandidates)
	  if (includes(c,cones[cL])) {
	    isConeMaximal[cL] = false;
	    for(auto rL : cones[cL]) 
	      maxConesContainingRay[rL].erase(cL);
	  }
	set<ConeLabel> largerCandidates(maxConesContainingRay[*c.begin()]);
	for (auto it = ++c.begin(); it != c.end(); ++it) {
	  set<ConeLabel> intersection;
	  set_intersection(largerCandidates.begin(), largerCandidates.end(), 
			   maxConesContainingRay[*it].begin(), maxConesContainingRay[*it].end(),
			   std::inserter(intersection, intersection.begin()));
	  largerCandidates = intersection;
	}
	isMaximal = largerCandidates.empty(); 
	/*auto it = cones.begin(); 
	while(isMaximal and it != cones.end()) {
	  if (includes(c,*it)) {
	    it = cones.erase(it);
	    std::cerr << "erased non-maximal" << std::endl;
	  }
	  else if (includes(*it,c)) {
	    isMaximal = false;
	    std::cerr << "detected non-maximal" << std::endl;
	  }
	  else ++it;  
	  }*/
      }
      isConeMaximal.push_back(isMaximal);
      if (isMaximal) 
	for(auto rL : c) 
	  maxConesContainingRay[rL].insert(cones.size()-1);
    }
    std::vector<PostCone> maxCones;
    for(int i=0; i<cones.size(); i++)
      if(isConeMaximal[i]) 
	maxCones.push_back(cones[i]); 
    cones = maxCones;
  }
  
  // SERVICE METHODS
  int nRays() { return rays.getHeight(); }
  bool isFinite(RayLabel r) { return not rays[r][0].isZero(); }

  void serialize(ostream& os) const
  {
    os << "DIMENSION " << dimension << std::endl;
    os << "RAYS" << std::endl;
    os << rays;
    os << std::endl;
    for(auto& c : cones) {
      os << "CONE";
      for(auto& r : c) 
	os << " " << r;
      os << std::endl;
    }
  }
  
  PostComplex subcomplex(const Component& component) {
    std::vector<RayLabel> oldLabels(component.begin(),component.end());
    std::vector<RayLabel> newLabels(nRays(),-1); // -1 = has no new label
    for(int i=0; i<oldLabels.size(); i++)
      newLabels[oldLabels[i]] = i;
    PostComplex ret(dimension);
    for(auto r : oldLabels)
      ret.rays.appendRow(rays[r]);
    for(auto& c: cones){
      if(std::all_of(c.begin(),c.end(),[&](RayLabel r){return newLabels[r]>=0;})) {
	PostCone newCone;
	for(auto r : c) 
	  newCone.insert(newLabels[r]);
	ret.cones.push_back(newCone);
      }
    }
    return ret;     
  }

  std::set<shared_ptr<Component>> connectedComponents()
  {
    using ComponentLabel = RayLabel;
    std::set<shared_ptr<Component>> components;
    std::unordered_map<RayLabel,std::shared_ptr<Component>> m;
    for(RayLabel r=0; r < nRays(); r++) 
      if (isFinite(r)) {
	Component newComponent;
	newComponent.insert(r);
	newComponent.number = r;
	auto newComponentPtr = make_shared<Component>(newComponent);
	components.insert(newComponentPtr);
	m[r] = newComponentPtr;
      }
    for(const auto& c : cones) {
      auto firstFinite = std::find_if(std::begin(c), std::end(c), [&](RayLabel r){return isFinite(r);});
      if (firstFinite == std::end(c)) 
	std::cerr << "no finite ray!!!" << std::endl;
      else {
	auto componentPtr = m[*firstFinite];
	for(auto r : c) {
	  if(isFinite(r)) {
	    auto rComponentPtr = m[r];
	    if(rComponentPtr != componentPtr) { 
	      // std::cerr<< "merging " << rComponentPtr->number << " and " << componentPtr->number << std::endl;
	      for(auto s : *rComponentPtr)
		if(isFinite(s)) 
		  m[s] = componentPtr;
	      componentPtr->merge(*rComponentPtr);
	      /*
	      for(auto comp: components) 
		std::cerr << comp->number << " ";
	      std::cerr << std::endl;
	      */
	      components.erase(rComponentPtr);
	      //std::cerr<< "#components = " << components.size() << std::endl;
	      /*
		for(auto comp: components) 
		std::cerr << comp->number << " ";
		std::cerr << std::endl;
	      */
	    } 
	  } else componentPtr->insert(r);
	}
      }
    }
    return components;
  }

  void printComponents(std::set<shared_ptr<Component>> components, bool saveNoncomets)
  {
    std::cerr << "-- Found " << components.size() << " component(s)\n";
    bool greatSuccess = true; 
    int failedCount = 0;
    for (auto component : components) {
      //std::cerr  << *component;
      std::cout << "#rays = " << std::setfill('0') << std::setw(4) << component->size() << std::endl;
      gfan::Matrix<mvtyp> unboundedDirections(0,dimension);
      for(RayLabel r : *component)
        if(rays[r][0].isZero())
          unboundedDirections.appendRow(rays[r]);
      std::cout << "#unbounded = " << std::setfill('0') << std::setw(4) << unboundedDirections.getHeight() << std::endl;      
      GeneratedCone<mvtyp> recessionCone(unboundedDirections.transposed()); // global recession cone
      int dimLineality = recessionCone.getDimensionOfLinealitySpace();
      std::cout << "Lineality dim: " << dimLineality << std::endl
		<< "Cone dim: " << std::setfill('0') << std::setw(2) << recessionCone.getDimension() << std::endl;	
      if (dimLineality == 0)
	std::cerr << "GREAT SUCCESS!!! :)\n";
      else {
	std::cerr << "FAILED... :(\n";
	if(saveNoncomets) {
	  std::fstream f;
	  f.open(std::to_string(component->number) + ".component.out",std::fstream::out);
	  subcomplex(*component).serialize(f);
	  f.close();
	}
      }
    }
  }
};//end class Postcomplex

template<class mvtyp> 
  ostream& operator<<(ostream& os, const PostComplex<mvtyp>& pc) 
  { 
    //if (pc.C!=nullptr) os << pc.C->toString(FPF_cones|FPF_maximalCones); 
    pc.serialize(os);
    return os;
  }

} // namespace gfan

#endif /* GFANLIB_app_anton_H_ */
