/*
 * divisionobject.h
 *
 *  Created on: May 24, 2023
 *      Author: anders
 */

#ifndef SRC_DIVISIONOBJECT_H_
#define SRC_DIVISIONOBJECT_H_

#include <queue>
#include <set>
#include <memory>
//#include <ranges>
//#include <range/v3/all.hpp>
#include "log.h"
#include "polynomial.h"
#include "packedmonomial.h"
#include "wallideal.h"
#include "gfanlib_frequencytable.h"
#include "gfanlibglue.h"
#include "gfanlib_tableau.h"

//abstract class
class DivisionObject{
public:
	//In general term orderings are not needed for finding normal forms.
	//However, things are speeded up if one is known - and the remainder can be computed degree by degree.
	//A good specialisations of this abstract class will however us LP to compute its own weight vector
	//to be tie-broken lexicographically thereby providing
	//a compatible term order. This is because comparisons can then be made using 64-bit machine
	//integer compoarisons on the packed exponent vector.
	//It may be a good idea to make the weight vector close to a grading in which the divisors
	//are homogeneous to give lower upper bounds on the exponent vectors needed through the division.
	//Therefore it is expected to add such argument to the constructor in the future.
	DivisionObject(PolynomialSet &g_, TermOrder const &termOrder){};
	virtual void reduce(Polynomial &f)=0;
	// It is allowed to return the zero vector if no grading is used. This should cause no trouble as this grading is only used as a hint
	virtual IntegerVector getUsedGrading()const //makes only sense for some DivisionObjects. It is probably better to make a subclass containing this method.
	{
		assert(0);
		return IntegerVector();
	}

//Would be very useful:
	//void disableNumber();
	//void replaceNumber();
};

class DivisionObjectNaive: public DivisionObject{
	PolynomialSet &g;
	TermOrder const &termOrder;
public:
	DivisionObjectNaive(PolynomialSet &g_, TermOrder const &termOrder_):
		DivisionObject(g_,termOrder_),
		g(g_),
		termOrder(termOrder_)
	{
	}
	void reduce(Polynomial &f);
	virtual IntegerVector getUsedGrading()const
	{
		return IntegerVector(g.getRing().getNumberOfVariables());
	}
};


template<int NWORDS>
class Dividend{
	std::priority_queue<
	    std::pair<PackedMonomial<NWORDS>,FieldElement>
	> terms;

};

template<int NWORDS>
class DivisionPolynomial{
public:
	vector<pair<PackedMonomial<NWORDS>,FieldElement>> terms;
//	vector<PackedMonomial<2>> monomials;
//	vector<FieldElement> coefficients;
public:
	DivisionPolynomial(PacMan const &pacMan, Polynomial const &f, IntegerVector const &w)//should w be stored in pacMan instead?
	{
		int nTerms=f.numberOfTerms();
		terms.reserve(nTerms);
		for(auto i=f.terms.begin();i!=f.terms.end();i++)
			terms.push_back(pair<PackedMonomial<NWORDS>,FieldElement>(PackedMonomial<NWORDS>(dotLong(w,i->first.exponent),i->first.exponent,pacMan),i->second));
//		monomials.reserve(nTerms);
//		coefficients.reserve(nTerms);
		//sort(monomials.begin(),monomials.end());
		sort(terms.begin(),terms.end(),[](const pair<PackedMonomial<NWORDS>,FieldElement> & a, const pair<PackedMonomial<NWORDS>,FieldElement> & b) -> bool
				{return b.first < a.first;});
//		ranges::v3::sort(ranges::view::zip(monomials, coefficients));
		//ranges::sort(ranges::view::zip(monomials, coefficients),
		 //                std::less<>{},
		  //               [](const auto& t) -> decltype(auto) { return std::get<0>(t); });
	}
};

//static gfan::FrequencyTable compF("operator<");
// Object for priority queue
template<int NWORDS>
class MultiTerm{
public:
	MultiTerm *next;
private:
	DivisionPolynomial<NWORDS> const *p;
	PackedMonomial<NWORDS> m;
//	FieldElement coefficient;
	int i;
	shared_ptr<FieldElement> coefficientp; // should be changed to unique_ptr, but this seems to require implementing a better priority queue
	//scalar multiplier
	//monomial multiplier
public:
//	MultiTerm():
//		coefficientp(nullptr)
//	{
//	}
	MultiTerm(DivisionPolynomial<NWORDS> &p_, PackedMonomial<NWORDS> const &m_, FieldElement coefficient_):
		p(&p_),
		m(m_),
//		coefficient(coefficient_),
		coefficientp(make_shared<FieldElement>(coefficient_)),
		i(0),
		next(0)
	{
	}
	bool isZero()const
	{
		return i==p->terms.size();
	}
	int	numberOfTermsLeft()const
	{
		return p->terms.size()-i;
	}
	bool operator<(MultiTerm const &b)const
	{
//		compF.record();
//		assert(!isZero());
//		assert(!b.isZero());
		return p->terms[i].first.sumCompare(m, b.p->terms[b.i].first, b.m);
	}
	PackedMonomial<NWORDS> getLeadingMonomial()const
	{
//		assert(!isZero());
		return m.sumNoOverflow(p->terms[i].first);
	}
	bool isLeadingMonomialEqualTo(PackedMonomial<NWORDS> const &b)const
	{
//		assert(!isZero());
		return b.isEqualToSum(p->terms[i].first,m);
	}
	FieldElement getLeadingCoefficient()const
	{
//		assert(!isZero());
		return (*coefficientp)*p->terms[i].second;
	}
	void addLeadingCoefficientTo(FieldElement &v)const
	{
		v.madd(*coefficientp,p->terms[i].second);
	}
	void popLeadingTerm()
	{
		assert(!isZero());
		assert(i<p->terms.size());
		i++;
	}
	void print(PacMan const &pacMan)const
	{
		cerr<<"Number of terms:"<<p->terms.size()-i;
		debug<<m.extractExponent(pacMan)<<"\nTimes";
		for(int j=i;j<p->terms.size();j++)
		{
			debug<<" "<<p->terms[j].first.extractExponent(pacMan)<<"\nTimes";
		}
	}
	//	bool isLeadingMonomialEqualTo();
};

template<int NWORDS>
class PriorityQueue
{
public:
	using Priority=PackedMonomial<NWORDS>;
	using Data=unique_ptr<MultiTerm<NWORDS>>;
private:
	class Element{
	public:
		Priority priority;
		Data data;
		Element(Priority priority_, Data&& data_):
			priority(priority_),
			data(std::move(data_))
		{
//			assert(data_.get()==nullptr);
		}
		Element()
		{
		}
	};
	vector<Element> queue;

	void insertAtEmptySpot(int emptySpotIndex, Priority p, Data && d)
	{
		int parentIndex=emptySpotIndex>>1;
		if(parentIndex && queue[parentIndex].priority<p)
			{
				queue[emptySpotIndex]=std::move(queue[parentIndex]);
				insertAtEmptySpot(parentIndex,p,std::move(d));
			}
		else
			queue[emptySpotIndex]=Element(p,std::move(d));
	}
	void __attribute__ ((noinline)) adjustEntry(int i, Element && value)//ith index is empty, but we pretend it contains value//by moving entry down in the tree
	{
		int child1=i+i;
		int child2=child1+1;
		if(child1>=queue.size())//no child 1
		{
			queue[i]=std::move(value);
//			assert(value.data.get()==nullptr);
			return;
		}
		if(child2>=queue.size())//no child 2
		{
			if(value.priority<queue[child1].priority)
			{
				queue[i]=std::move(queue[child1]);
				queue[child1]=std::move(value);
//				assert(value.data.get()==nullptr);
				return; // There is no need to recurse, since child1 cannot have children without having a sibling.
			}
			queue[i]=std::move(value);
//			assert(value.data.get()==nullptr);
			return;
		}
		else//child 1 and child 2
		{
			int child=child1+(queue[child1].priority<queue[child2].priority);
			if(value.priority<queue[child].priority)
			{
				std::move(queue[i])=std::move(queue[child]);//?
				adjustEntry(child,std::move(value));
			}
			else
				queue[i]=std::move(value);
//		assert(value.data.get()==nullptr);
		}
	}
public:
	void push(Priority p, Data &&d)
	{
		queue.resize(queue.size()+1);
		insertAtEmptySpot(queue.size()-1,p,std::move(d));
	}
	bool isEmpty()const
	{
		return queue.size()==1;
	}
	PriorityQueue():
		queue(1)//first entry is unused
	{
	}
	MultiTerm<NWORDS>* topDataGet()
	{
		return queue[1].data.get();
	}
	Priority topPriority()const
	{
		return queue[1].priority;
	}
	void changeTop(Priority p)
	{
		adjustEntry(1,std::move(Element(p,std::move(queue[1].data))));
	}
	Data pop()
	{
		Data ret=std::move(queue[1].data);
		auto temp=std::move(queue.back());
		queue.pop_back();
		adjustEntry(1,std::move(temp));
		return ret;
	}
	void consistencyCheck()
	{
		for(int i=1;i<queue.size();i++)
			assert(queue[i].data.get());
	}
};

template<int NWORDS>
class PriorityQueue2 // We try to implement this more efficiently without unique_ptr
{
	const int hashShift=8;
	const int numberOfHashEntries=0x10000>>hashShift;
public:
	using Priority=PackedMonomial<NWORDS>;
	vector<vector<pair<Priority,int> > > hashVector;
	using Data=MultiTerm<NWORDS>*; //Is it possible to not make this a pointer?
//	using Data=vector<MultiTerm<NWORDS>*>;
private:
	void releaseHash(Priority p, int i)
	{
		int hash=p.hash()>>hashShift;
		//std::cerr<<"Releasing Hash"<<hash<<" "<<i<<"\n";
		int index=-1;
		auto &v=hashVector[hash&(numberOfHashEntries-1)];
		for(int i=0;i<v.size();i++)
			if(p==v[i].first)
				{
					index=i;
					break;
				}
		assert(index>=0);
		//for(auto &a:v)std::cerr<<a.second<<"\n";
		//std::cerr<<"Erasing"<<index<<"\n";
		v.erase(v.begin()+index);
		//for(auto &a:v)std::cerr<<a.second<<"\n";
		//std::cerr<<"Done releasing\n";
	}
	void recordHash(Priority p, int i)
	{
		int hash=p.hash()>>hashShift;
//		std::cerr<<"Recording Hash"<<hash<<" "<<i<<"\n";
		hashVector[hash&(numberOfHashEntries-1)].push_back(pair<Priority,int>(p,i));
	}
	class Element{
	public:
		Priority priority;
		Data data;
		Element(Priority priority_, Data data_):
			priority(priority_),
			data(data_)
		{
//			assert(data_.get()==nullptr);
		}
		Element()
		{
		}
	};
	vector<Element> queue;

	void insertAtEmptySpot(int emptySpotIndex, Priority p, Data d)
	{
//		std::cerr<<"Inserting at"<<emptySpotIndex<<"\n";
		int parentIndex=emptySpotIndex>>1;
		if(parentIndex && queue[parentIndex].priority<p)
			{
				releaseHash(queue[parentIndex].priority,parentIndex);
				queue[emptySpotIndex]=queue[parentIndex];
				recordHash(queue[emptySpotIndex].priority,emptySpotIndex);
				insertAtEmptySpot(parentIndex,p,d);
			}
		else
		{
			queue[emptySpotIndex]=Element(p,d);
			recordHash(queue[emptySpotIndex].priority,emptySpotIndex);
		}
//		print();
	}
	void __attribute__ ((noinline)) adjustEntry(int i, Element const& value)//ith index is empty, but we pretend it contains value//by moving entry down in the tree
	{
		int child1=i+i;
		int child2=child1+1;

//		if(child2<queue.size()) // both children
		if(__builtin_expect(child2<queue.size(),1)) // both children
		{
			int child=child1+(queue[child1].priority<queue[child2].priority);
//			if(value.priority<queue[child].priority)
			if(__builtin_expect(value.priority<queue[child].priority,1))
			{

				releaseHash(queue[child].priority,child);
				queue[i]=queue[child];
				recordHash(queue[i].priority,i);

				adjustEntry(child,value);
			}
			else
			{
				queue[i]=value;
				recordHash(queue[i].priority,i);
			}
//		assert(value.data.get()==nullptr);
		}
		else
		{
			if(child1>=queue.size())//no child 1
			{
				queue[i]=value;
				recordHash(queue[i].priority,i);
//				assert(value.data.get()==nullptr);
				return;
			}
			else //no child 2
			{
				if(value.priority<queue[child1].priority)
				{
					releaseHash(queue[child1].priority,child1);
					queue[i]=queue[child1];
					recordHash(queue[i].priority,i);
					queue[child1]=value;
					recordHash(queue[child1].priority,child1);
	//				assert(value.data.get()==nullptr);
					return; // There is no need to recurse, since child1 cannot have children without having a sibling.
				}
				queue[i]=value;
				recordHash(queue[i].priority,i);
	//			assert(value.data.get()==nullptr);
				return;
			}
		}
	}
public:
	void dumpInfo()const
	{
			std::cerr<<"Queue size:"<<queue.size()-1<<"  ";
			set<PackedMonomial<NWORDS>> a;
			uint64_t nrep=0;
			for(int i=1;i<queue.size();i++)
			{
				a.insert(queue[i].priority);
				nrep+=queue[i].data->numberOfTermsLeft();
			}
			std::cerr<<"Number of represented terms:"<<nrep<<"  ";
			std::cerr<<"Unique initial terms:"<<a.size()-1<<"\n";
	}
	void maybeDumpInfo()const
	{
		return;
		thread_local int64_t t;
		t++;
		if(t>10000000)
		{
			dumpInfo();
			t=0;
		}
	}
	int hashIndex(Priority p)const // use hash to find index in priority queue with given priority
	{
		auto &v=hashVector[(p.hash()>>hashShift)&(numberOfHashEntries-1)];
		for(int i=0;i<v.size();i++)
			if(p==v[i].first)return v[i].second;
		return 0;
	}
	void consistencyCheck()const
	{
		for(int i=0;i<numberOfHashEntries;i++)
			if(hashVector[i].size())
		{
				std::cerr<<"Checking "<<i<<":";
			for(int j=0;j<hashVector[i].size();j++)
			{
				std::cerr<<" "<<j<<" ";
				int h=queue[hashVector[i][j].second].priority.hash()>>hashShift;
				if(h!=i)
				{
					std::cerr<<"Inconsistency in hash vector "<<i<<" entry "<<j<<"\n";
					assert(0);
				}
			}
			std::cerr<<"\n";
		}

		for(int i=1;i<queue.size();i++)
		{
			for(int j=i+1;j<queue.size();j++)
				if(queue[i].priority==queue[j].priority)
				{
					std::cerr<<"Entry "<<i<<" and "<<j<<" have the same priority!\n";
					assert(0);
				}
		}
	}
	void print()const
	{
		std::cerr<<"-------------------------------------------\n";
		std::cerr<<"Printing PriorityQueue:\n";
		std::cerr<<"Number of priorities:"<<queue.size()-1<<"\n";
		for(int i=1;i<queue.size();i++)
		{
			std::cerr<<i<<":\n";
			std::cerr<<"Priority:\n";
			queue[i].priority.print();
			std::cerr<<"\nHash:"<<(queue[i].priority.hash()>>hashShift)<<"\n";
			std::cerr<<"Data:\n";
			Data p=queue[i].data;
			for(int i=0;i<3;i++)
			{
				std::cerr<<std::hex<<(int64_t)p<<"\n";
				if(!p)break;
				p=p->next;
			}
		}
		std::cerr<<"hashVector:\n";
		int totalSize=0;
		for(int i=0;i<numberOfHashEntries;i++)
			if(hashVector[i].size())
		{
			std::cerr<<i<<":"<<hashVector[i].size()<<" ";
			totalSize+=hashVector[i].size();
			for(int j=0;j<hashVector[i].size();j++)
				std::cerr<<hashVector[i][j].second<<" ";
			std::cerr<<"\n";
		}
		std::cerr<<"total"<<totalSize<<"\n";
		assert(totalSize==queue.size()-1);
		std::cerr<<"Done printing PriorityQueue\n";
		std::cerr<<"-------------------------------------------\n";
		consistencyCheck();
	}
	int findIndex(Priority p)const //return position in queue of given priority
	{
		return hashIndex(p); //A naive implementation follows below

//		print();
		int i=1;
		int ret=0;
		while(i<queue.size())
		{
			if(queue[i].priority==p)
				{
					ret=i;
					break;
				}
			i++;
		}
//		std::cerr<<"ret:"<<ret<<"\n";
//		std::cerr<<"hashIndex(p):"<<hashIndex(p)<<"\n";

		assert(ret==hashIndex(p));
		return ret;
	}
	void push(Priority p, Data d) //Data will be owned by PriorityQueue2
	{
//		std::cerr<<"Push\n";
		int i=findIndex(p);
		if(i)
		{
//			std::cerr<<"Found Index\n";
			d->next=queue[i].data;
			queue[i].data=d;
		}
		else
		{
//			std::cerr<<"Resizing\n";
			queue.resize(queue.size()+1);
			assert(d->next==0);
			insertAtEmptySpot(queue.size()-1,p,d);
		}
	}
	bool isEmpty()const
	{
		maybeDumpInfo();
		return queue.size()==1;
	}
	PriorityQueue2():
		queue(1),//first entry is unused
		hashVector(numberOfHashEntries)
	{
	}
	~PriorityQueue2()
	{
		assert(isEmpty());//if not, delete should be called on elements.
	}
	Data topDataGet()
	{
		return queue[1].data;
	}
	Priority topPriority()const
	{
		return queue[1].priority;
	}
/*	void changeTop(Priority p)?????????????
	{
		adjustEntry(1,Element(p,queue[1].data));
	}*/
	void popNoDelete()
	{
//		print();
//		std::cerr<<"PopNoDelete\n";
		Data next=queue[1].data->next;
		if(next)
		{
//			std::cerr<<"Not null\n";//If more data
			queue[1].data=next;
		}
		else
		{
			releaseHash(queue[1].priority,1);
			if(0)
			{
				std::cerr<<"hashVector:\n";
				int totalSize=0;
				for(int i=0;i<numberOfHashEntries;i++)
					if(hashVector[i].size())
				{
					std::cerr<<i<<":"<<hashVector[i].size()<<" ";
					totalSize+=hashVector[i].size();
					for(int j=0;j<hashVector[i].size();j++)
						std::cerr<<hashVector[i][j].second<<" ";
					std::cerr<<"\n";
				}

			}
			if(queue.size()>2)
			{
//				std::cerr<<"More than one element\n";
				auto temp=queue.back();
				releaseHash(queue.back().priority,queue.size()-1);
				queue.pop_back();
				adjustEntry(1,temp);
			}
			else
			{
//				std::cerr<<"One element only\n";
				assert(queue.size()==2);
				queue.pop_back();
			}
		}
	}
	void pop()
	{
//		std::cerr<<"Pop\n";
		Data next=queue[1].data->next;
		if(next)
		{
			delete queue[1].data;
			queue[1].data=next;
		}
		else
		{
			releaseHash(queue[1].priority,1);
			delete queue[1].data;
			if(queue.size()>2)
			{
				auto temp=queue.back();
				releaseHash(queue.back().priority,queue.size()-1);
				queue.pop_back();
				adjustEntry(1,temp);
			}
			else
			{
				assert(queue.size()==2);
				queue.pop_back();
			}
		}
	}
	void consistencyCheck()
	{
//		for(int i=1;i<queue.size();i++)
//			assert(queue[i].data.get());
	}
};


static IntegerVector weightVector(PolynomialSet const &g, IntegerVector const *wHint=0)
{
	if(wHint && g.checkMarkings(*wHint))
		{
//		debug<<"Hint used.\n";
		return *wHint;
		}
	int n=g.getRing().getNumberOfVariables();
	auto inequalities=wallInequalities(g);
	for(int i=0;i<n;i++)
		inequalities.push_back(IntegerVector::standardVector(n,i));

	log1 debug<<"Number of inequalities"<<(int)inequalities.size()<<"\n";
#if 1
	gfan::Cone<gfan::CircuitTableInt32> C2(
			gfan::convertMatrix<gfan::CircuitTableInt32>(gfan::rowsToIntegerMatrix2(inequalities/*,n*/)).transposed()
			);
	return toIntegerVector(C2.getRelativeInteriorPoint());
#else
	PolyhedralCone C(inequalities,IntegerVectorList());
	return C.getRelativeInteriorPoint();
#endif
}


//static gfan::FrequencyTable reduceFT("num poly");



template<int NWORDS, class FieldElementImplementation=void>
class DivisionObjectFast : public DivisionObject{
	PacMan pacMan;
	vector<DivisionPolynomial<NWORDS>> g;
	vector<PackedMonomial<NWORDS>> gInitial;
	IntegerVector w;
	int64 wmax;
	PolynomialSet original;
	static int64 degreeBound(PolynomialSet &g_, IntegerVector const &w)
	{
		int64 ret=0;
		for(auto &f:g_)
			ret=max(ret,f.degree(w));
		return ret;
	}
	static vector<int> degreeBoundToBitBounds(int64 bound, IntegerVector const &w)
	{
		IntegerVector ret=w;
		for(int i=0;i<w.size();i++)
			ret[i]=(bound/w[i]);
		return PacMan::bitsNeeded(ret);
	}
	/** N0tt neeeded! Computes what w-degree bound guarantees that any monomial with at most
	 this w-degree will fit in the packed monomial structure specified by bitBounds.
	 */
	int64 bitBoundsToDegreeBound(vector<int> const &bitBounds, IntegerVector const &w)
	{
		int64 ret=1000000000;
		for(int i=0;i<bitBounds.size();i++)
			ret=min(ret,w[i]*((((int64)1)<<bitBounds[i])-1));
		return ret;
	}
/*
	vector<int> bounds(PolynomialSet &g_)
	{
		PolynomialRing r_=g_.getRing();
		IntegerVector maxexp(r_.getNumberOfVariables());
		for(auto &i:g_)
			maxexp=max(maxexp,i.degreeVector());
		vector<int> maxexp2;for(int i=0;i<maxexp.size();i++)maxexp2.push_back(maxexp[i]);
		return maxexp2;
	}*/
	void computeDegreeBoundsAndBuildPacked(int64 neededDegree=0)
	{
		int64 b=degreeBound(original,w);
	//	cerr<<"degreeBound"<<b<<"\n";
		if(neededDegree>b)b=neededDegree;
	//	cerr<<"degreeBound"<<b<<"\n";
		auto bb=degreeBoundToBitBounds(b,w);
	//	for(auto&a:bb)cerr<<" "<<a;
	//	cerr<<"\n";
wmax=b;
//		wmax=bitBoundsToDegreeBound(bb,w);// <--problemet må være her!
//		cerr<<"wmax"<<wmax<<"\n";
pacMan=PacMan(original.getRing(),bb,31-16,true/*,2*/);
//pacMan=PacMan(original.getRing(),bb,31-20,true,2);
		g.clear();
		g.reserve(original.size());
		for(auto &f:original)g.emplace_back(pacMan,f,w);
		gInitial.clear();
		gInitial.reserve(g.size());
		for(int i=0;i<g.size();i++)
			{
				assert(g[i].terms.size());
				gInitial.push_back(g[i].terms[0].first);
			}

//		static gfan::FrequencyTable nTermsFT("nTermsDivisors");
//		for(auto &f:original)nTermsFT.record(f.numberOfTerms());
//	debug<<"Grading:"<<w<<"\n";
//		cerr<<"Wmax:"<<wmax<<"\n";
	}
public:

	void print(PriorityQueue2<NWORDS> /*const*/ &a, PacMan const &pacman)
	{
		cerr<<"Printing queue (first element)\n";
		a.topDataGet()->print(pacMan);
//		for(auto &b:a.c)
		{
//			b.print(pacMan);
		}
		cerr<<"Done printing queue\n";
	}

	DivisionObjectFast(PolynomialSet &g_, IntegerVector const *wHint, TermOrder const &termOrder):
//		pacMan(g_.getRing(),bounds(g_),31),
		w(weightVector(g_,wHint)),
		original(g_),
		DivisionObject(g_,termOrder)
	{
/*		debug<<original;
		original.mark(termOrder);
		debug<<original;
		autoReduce(&original, termOrder);
		debug<<original;*/
//		original.sort( []( const Polynomial &a, const Polynomial &b ) { return a.numberOfTerms() < b.numberOfTerms(); } );
		original.sort(PolynomialCompareMarkedTerms(TotalDegreeTieBrokenTermOrder(termOrder)));
//		original.sort( [=]( const Polynomial &a, const Polynomial &b ) { return a.degree(w) < b.degree(w); } );
//		original.sort( [=]( const Polynomial &a, const Polynomial &b ) { return a.totalDegree()(w) > b.totalDegree(w); } );





		assert(g_.checkMarkings(w));
		computeDegreeBoundsAndBuildPacked();
//		debug<<g_;
//		pacMan.print();
//		debug.printVector(w);pout.printNewLine();
//		g.reserve(g_.size());
//		for(auto &f:g_)g.emplace_back(pacMan,f,w);
	}
	IntegerVector getUsedGrading()const
	{
		return w;
	}
	void reduce(Polynomial &f)
	{
//		cerr<<"reduced on "<<typeid(FieldElementImplementation).name()<<"\n";
		assert(pacMan.wordData.size()>=2);
		assert(g.size()==gInitial.size());
		if(f.isZero())return;
		if(f.degree(w)>wmax)
		{
			computeDegreeBoundsAndBuildPacked(max((wmax+1)*2,f.degree(w)));
	//		cerr<<f.degree(w)<<">"<<wmax<<"\n";
	//		debug<<w<<"\n";
	//		debug<<f<<"\n";
	//		pacMan.print();
		}
		assert(f.degree(w)<=wmax);
//		debug<<"Reduce on"<<f.toString(false)<<"\n";
log3	debug<<"Reduce on poly with "<<(int)f.terms.size()<<" terms\n";
		DivisionPolynomial<NWORDS> F(pacMan,f,w);

		PriorityQueue2<NWORDS> thePriorityQueue;
//		thePriorityQueue.print();
		f=f.getRing().zero();//should this code be able to throw?
		auto temp=make_unique<MultiTerm<NWORDS>>(F,PackedMonomial<NWORDS>(),f.getRing().getField().zHomomorphism(1));
		PackedMonomial<NWORDS> tempPriority=temp->getLeadingMonomial();
//		a.push(tempPriority,std::move(temp));
		thePriorityQueue.push(tempPriority,temp.get());temp.release();
		while(!thePriorityQueue.isEmpty())
		{
			//	reduceFT.record(a.size());
			//FieldElementImplementation
			auto c=f.getRing().getField().zHomomorphism(0);

//			FieldElementImplementation C;

			PackedMonomial<NWORDS> m=thePriorityQueue.topPriority();
			log3{
				static int i;i++;
				if(i>10000)
					{
					debug<<m.extractExponent(pacMan)<<(int)m.extractWeight(pacMan)<<"\n";
					i=0;
					}
			}
//			a.top().print(pacMan);
			while(!thePriorityQueue.isEmpty()&&thePriorityQueue.topPriority()==m)
			{
	//			a.consistencyCheck();
				MultiTerm<NWORDS> *temp(thePriorityQueue.topDataGet());
//				c+=temp->getLeadingCoefficient();
				temp->addLeadingCoefficientTo(c);
				temp->popLeadingTerm();
				if(!temp->isZero())
				{
					thePriorityQueue.popNoDelete();
					temp->next=0;
					thePriorityQueue.push(temp->getLeadingMonomial(),temp);
				}
				else
				{
					thePriorityQueue.pop();
				}
			}
			if(c.isZero())continue;
			int i=g.size();
			static gfan::FrequencyTable T("DIVISON STEPS");
			T.record(1);
			for(int I=0;I<g.size();I++)
			{
//				if(g[i].terms[0].first.dividesNoAdditionalBit(m,pacMan))break;
//				if(g[i].initialMonomial.dividesNoAdditionalBit(m,pacMan))break;
				if(gInitial[I].dividesNoAdditionalBit(m,pacMan))
#if 0
					{
						i=I;
						break;
					}
#else
				{
					if(i==g.size())
						i=I;
					else
/*					{
						if(g[I].terms.size()<3 || g[i].terms.size()<3)
						{
							if(g[I].terms.size()<g[i].terms.size())
								i=I;
						}
						else
						{
							if(g[I].terms.size()<2)
								i=I;
							else if(g[I].terms.size()>=2)
								{
									if(g[i].terms.size()>1)
										if(g[i].terms[0].first.sumCompare(
												g[I].terms[1].first,
												g[I].terms[0].first,
												g[i].terms[1].first)
												)
										{i=I;}
								}

						}
					}
					if(0)*/
					{
						if(g[I].terms.size()<2)
							i=I;
						else if(g[I].terms.size()>=2)
//						if(g[i].terms.size()>g[I].terms.size())
//						i=I;
//						else
//							if(g[i].terms.size()==g[I].terms.size())
							{
								if(g[i].terms.size()>1)
									if(g[i].terms[0].first.sumCompare(
											g[I].terms[1].first,
											g[I].terms[0].first,
											g[i].terms[1].first)
											)
									{i=I;}
							}
					}
				}
#endif
			}
			if(i!=g.size())
			{
				c=c.getField()->zHomomorphism(0)-c;
				m.divideBy(g[i].terms[0].first,pacMan);
				unique_ptr<MultiTerm<NWORDS>> temp=make_unique<MultiTerm<NWORDS>>(g[i],m,c*g[i].terms[0].second.inverse());
				temp->popLeadingTerm();
				if(!temp->isZero())
				{
					auto tempPriority=temp->getLeadingMonomial();
//					a.push(tempPriority,std::move(temp));
					thePriorityQueue.push(tempPriority,temp.get());temp.release();
				}
			}
			else
			{
//				debug<<"Nothing divides:"<<Term(c,Monomial(f.getRing(),m.extractExponent(pacMan)))<<"\n";
				f+=Term(c,Monomial(f.getRing(),m.extractExponent(pacMan)));
			}
		}
//	debug<<"New remainder:"<<f<<"\n";
	}
};

unique_ptr<class DivisionObject> DivisionObjectFactory(Field const &k, PolynomialSet &g_, IntegerVector const *wHint, TermOrder const &termOrder);

#endif /* SRC_DIVISIONOBJECT_H_ */
