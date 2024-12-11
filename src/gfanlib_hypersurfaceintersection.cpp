/*
 * gfanlib_hypersurfaceintersection.cpp
 *
 *  Created on: 12 Mar 2020
 *      Author: anders
 */

#include <mutex>
#include "gfanlib_hypersurfaceintersection.h"
#include "gfanlib_paralleltraverser.h"
#include "gfanlib_memoryresource.h"

namespace gfan{
	template<class typ>
	void codimension1ConesB(Matrix<typ> nonstrict, Matrix<typ> strict, pmrvector<HalfOpenCone<typ> > &collector, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
	{
		assert(nonstrict.getHeight()==strict.getHeight());
		HalfOpenCone<typ> C(nonstrict,Matrix<typ>(strict.getHeight(),0,mr2),strict,mr2);
//		std::cerr<<"DUAL"<<C.lifted.dualCone.toString()<<"\n";
//		std::cerr<<"NonStrict:\n"<<nonstrict.toString()<<"Strict:\n"<<strict.toString()<<"\n";
//		for(int i=0;i<nonstrict.getWidth();i++)
		for(int i=nonstrict.getWidth()-1;i>=0;i--)
		{
//			std::cerr<<"Iterated:"<<C.toString()<<"\n";
			HalfOpenCone<typ> K(mr);K=C;
//			std::cerr<<"DUAL"<<K.lifted.dualCone.toString()<<"\n";
			K.makeEquation(i);
//			std::cerr<<"DUAL"<<K.lifted.dualCone.toString()<<"\n";
			if(!K.isEmpty())
			{
//				std::cerr<<"ADDING:\n:";/*<<*/
//				K.toString();//<<"\n";						//This line has side effects
				collector.emplace_back(std::move(K));
			}
//			if(i+1<nonstrict.getWidth())
			if(i>0)
			{
//				std::cerr<<"C BEFORE"<<C.toString()<<"\n";
				C.makeStrictInequality(i);
//				std::cerr<<"C AFTER"<<C.toString()<<"\n";
			}
		}
	}
/*	void codimension1ConesC(Matrix<CircuitTableInt32> choosable, Matrix<CircuitTableInt32> nonstrict, Matrix<CircuitTableInt32> strict, Matrix<CircuitTableInt32> equations, std::vector<HalfOpenCone<CircuitTableInt32> > &collector, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
	{
		HalfOpenCone<CircuitTableInt32> C(combineLeftRight(choosable,nonstrict),strict,equations);
		for(int i=0;i<choosable.getWidth();i++)
		{
//			std::cerr<<"Iterated:"<<C.toString()<<"\n";
			HalfOpenCone<CircuitTableInt32> K=C;
			K.makeEquation(i);
			if(!K.isEmpty())
				collector.emplace_back(std::move(K));
			if(i+1<chooseable.getWidth())
				C.makeStrictInequality(i);
		}
	}*/
	void codimension1Cones(Matrix<CircuitTableInt32> nonstrict, Matrix<CircuitTableInt32> strict, std::vector<HalfOpenCone<CircuitTableInt32> > &collector, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
//			ResourceInvariant a(mr2);
//			std::cerr<<"codim1"<<nonstrict.toString()<<strict.toString()<<"\n";
			if(nonstrict.getWidth()>0)
			{
				//only do the following if equation defines a facet?

										HalfOpenCone<CircuitTableInt32> temp(//HERE
												nonstrict.submatrixColumns(0,nonstrict.getWidth()-1,mr),
												nonstrict.submatrixColumns(nonstrict.getWidth()-1,nonstrict.getWidth(),mr),
												strict,mr2);//TO HERE
					//					std::cerr<<temp.toString()<<"contains origin"<<temp.containsOrigin()<<"\n";
										collector.emplace_back(temp,mr
										);

										//				std::vector<HalfOpenCone<CircuitTableInt32> > ret=
						codimension1Cones(
						nonstrict.submatrixColumns(0,nonstrict.getWidth()-1,mr2),
						combineLeftRight(strict,nonstrict.submatrixColumns(nonstrict.getWidth()-1,nonstrict.getWidth(),mr2),mr2),
						collector,mr,mr2
						);

//						std::cerr<<"NONSTRICT:"<<									combineOnTop(nonstrict,
//								-Matrix<CircuitTableInt32>::rowVectorMatrix(nonstrict[nonstrict.getHeight()-1]))<<"ADDING"<<collector.back().toString()<<"\n\n";
			}
//			std::cerr<<"codim1returning\n";
//			else
//			return std::vector<HalfOpenCone<CircuitTableInt32> >();
		}
	IntegerVectorList cast(Matrix<CircuitTableInt32> const A)
	{
		IntegerVectorList ret;
		for(int i=0;i<A.getHeight();i++)
		{
			IntegerVector a(A.getWidth());
			for(int j=0;j<a.size();j++)
				a[j]=A[i][j].v;
			ret.push_back(a);
		}
		return ret;
	}
	template<class typ>
	bool compare(Vector<typ> const &order, Vector<typ> const &v)
	{
		auto p=dot(order,v);
		if(p.isNegative())return 1;
		if(p.isPositive())return 0;
		for(int i=0;i<v.size();i++)
		{
			if(v[i].isNegative())return 1;
			if(v[i].isPositive())return 0;
		}
		return 0;
	}
	template<class typ>
	bool compare(Matrix<typ> const &order, Vector<typ> const &a, MR *mr=get_default_resource())
	{
		for(int j=0;j<order.getHeight();j++)
		{
			int s=compare(order[j].toVector(mr),a);
			if(s<0)return true;
			if(s>0)return false;
		}
		return false;
	}
	template<class typ>
	bool vertexCompare(Vector<typ> const &a, Vector<typ> const &b, Matrix<typ> const &order)
	{
		for(int j=0;j<order.getHeight();j++)
		{
			int s=dot(a-b,order[j]).sign();
			if(s<0)return true;
			if(s>0)return false;
		}
		return false;
	}
	template<class typ>
	bool isPerpendicular(Vector<typ> const &a, Matrix<typ> const &generators)
	{
		for(int j=0;j<generators.getHeight();j++)
			if(dot(a,generators[j]).sign()!=0)return false;
		return true;
	}
	template<typename typ>
	int vertexIndex(Matrix<typ> const &exponents, Matrix<typ> const &order, int &numberOfTimesAttained, MR *mr=get_default_resource())
	{
		assert(exponents.getHeight());
		int best=0;
//		std::cerr<<"newBest:"<<0;
		numberOfTimesAttained=1;
		for(int i=1;i<exponents.getHeight();i++)
		{
			if(vertexCompare(exponents[best].toVector(mr),exponents[i].toVector(mr),order))
			{
//				std::cerr<<"newBest:"<<i;
				best=i;
				numberOfTimesAttained=1;
			}
			else if(!vertexCompare(exponents[i].toVector(mr),exponents[best].toVector(mr),order))
			{
				numberOfTimesAttained++;
//				std::cerr<<"equalvalue"<<i;
			}
		}
		return best;
	}
//	HalfOpenCone<CircuitTableInt32> && T(HalfOpenCone<CircuitTableInt32> && a){stackResource.dump();return std::move(a);}

// ---------------------
// VertexPairToEdgeTable
// ---------------------

	VertexPairToEdgeTable::VertexPairToEdgeTable(int numberOfVertices):
			table(numberOfVertices,numberOfVertices),
			nextFree(0)
	{
		for(int i=0;i<numberOfVertices;i++)
			for(int j=0;j<numberOfVertices;j++)
				table[i][j]=-1;
	}

	int VertexPairToEdgeTable::allocate(int i, int j)
	{
		table[i][j]=nextFree;
		table[j][i]=nextFree;
		return nextFree++;
	}

	int VertexPairToEdgeTable::lookup(int i, int j)const
	{
		return table[i][j];
	}

	string VertexPairToEdgeTable::toString()const
	{
		return table.toString();
	}

	// ---------------------
	// RelationTableLayout
	// ---------------------
	RelationTableLayout::RelationTableLayout():
			nextFree(0)
	{
	}
	void RelationTableLayout::allocate(int size)
	{
		startOffsets.push_back(nextFree);
		nextFree+=size;
		endOffsets.push_back(nextFree);
	}
	int RelationTableLayout::numberOfUsedBits()const
	{
		return nextFree;
	}

	// -------------
	// RelationTable
	// -------------
	RelationTable::RelationTable(RelationTableLayout const &layout_):
			layout(layout_),
			data(layout_.numberOfUsedBits())
	{
			for(int i=0;i<data.size();i++)data[i]=true;
	}
	RelationTable RelationTable::intersection(RelationTable const &b)const
	{
		assert(b.layout.numberOfUsedBits()==layout.numberOfUsedBits());
		RelationTable ret(layout);
		for(int i=0;i<layout.numberOfUsedBits();i++)
			ret.data[i]=data[i]&&b.data[i];
		return ret;
	}
	void RelationTable::mark(int polytopeIndex, int edgeIndex, bool b)
	{
		data[layout.startOffsets[polytopeIndex]+edgeIndex]=b;
	}
	bool RelationTable::lookUp(int polytopeIndex, int edgeIndex)const
	{
		return data[layout.startOffsets[polytopeIndex]+edgeIndex];
	}
	int RelationTable::count(int polytopeIndex)const
	{
		return std::count(data.begin()+layout.startOffsets[polytopeIndex],data.begin()+layout.endOffsets[polytopeIndex],true);
	}
	int RelationTable::bestPolytopeIndex(vector<int> const &used)const
	{
		int best;int bestIndex=-1;
		assert(used.size()==layout.startOffsets.size());
		for(int i=0;i<used.size();i++)
			if(!used[i])
				if(bestIndex==-1 || count(i)<best)
				{
					bestIndex=i;
					best=count(i);
				}
		assert(bestIndex!=-1);
		return bestIndex;
	}


	// ------------
	template<typename typ>
	Matrix<typ> normalConeInequalities(Matrix<typ> const &vertices, int i)
		{
			Matrix<typ> ret(0,vertices.getWidth());
			ret.reserveRows(vertices.getHeight()-1);
			for(int j=0;j<vertices.getHeight();j++)
				if(j!=i)
					ret.appendRow(vertices[i].toVector()-vertices[j].toVector());
//			std::cerr<<"Returning:"<<ret<<"\n";
			return ret.transposed();
		}
	template<typename typ>
	HalfOpenCone<typ> relativeInterior(Cone<typ> C)
		{
			return HalfOpenCone<typ>(
					Matrix<typ>(C.getAmbientDimension(),0),
					C.getOrthogonalComplement(),
					C.getFacetNormals());
		}

	// ------------
	// PolytopeData
	// ------------
	template<typename typ>
	PolytopeData<typ>::PolytopeData(Matrix<typ> const &generators_):
			edgeTable(0),
			generators(generators_)
	{
			// Find vertices
			for(int i=0;i<generators_.getHeight();i++)
				{
					Cone<typ> C(normalConeInequalities(generators,i));
//					std::cerr<<"Ambient dimension"<<C.getAmbientDimension()<<C.getDimension()<<"\n";
					if(C.getDimension()==C.getAmbientDimension())
					{
						vertices.push_back(i);
//						vertexNormalCones.push_back(C);
					}
				}
//			std::cerr<<"NumberOfVertices"<<vertices.size()<<"\n";
			edgeTable=VertexPairToEdgeTable(vertices.size());
			for(int i=0;i<vertices.size();i++)
				for(int j=0;j<i;j++)
				{
					Cone<typ> C(combineLeftRight(normalConeInequalities(generators,i),normalConeInequalities(generators,j)));
					{
					if(C.getDimension()==C.getAmbientDimension()-1)
						{
							edgeNormalCones.push_back(relativeInterior(C));
							edgeTable.allocate(i,j);
							edges.push_back(pair<int,int>(i,j));
						}
					}
				}
//			std::cerr<<edgeTable.toString();
	}
	template<typename typ>
	void PolytopeData<typ>::initializeRelationTables(RelationTableLayout const &layout)
	{
		for(auto &e:edges)relationTables.emplace_back(layout);
	}
	template<typename typ>
	string PolytopeData<typ>::toString(RelationTableLayout const &layout)const
	{
		stringstream s;
		for(auto &r:relationTables)
		{
			int I=0;
			for(int i=0;i<r.data.size();i++)
			{
				if(I<layout.startOffsets.size()&&layout.startOffsets[I]==i){s<<"|";I++;}
				s<<r.data[i];
			}
			s<<"\n";
		}
		s<<"\n";
		return s.str();
	}


	// ------------------------
	// PolytopeIntersectionData
	// ------------------------
	template<typename typ>
	std::string PolytopeIntersectionData<typ>::toString()
	{
		stringstream s;
		for(auto &p:polytopes)
			s<<p.toString(layout)<<"\n";
		return s.str();
	}
	template<typename typ>
	PolytopeIntersectionData<typ>::PolytopeIntersectionData(vector<Matrix<typ>> const &polytopes_)
	{
		polytopes.reserve(polytopes_.size());
		for(auto &p:polytopes_)
			{
			polytopes.emplace_back(PolytopeData<typ>(p));
			layout.allocate(polytopes.back().edges.size());
			}
		for(auto &p:polytopes)
			p.initializeRelationTables(layout);
		for(int P=0;P<polytopes.size();P++)
			for(int Q=0;Q<P;Q++)
			{
				for(int E=0;E<polytopes[P].edges.size();E++)
					for(int F=0;F<polytopes[Q].edges.size();F++)
					{
						HalfOpenCone<typ> temp=intersection2(polytopes[P].edgeNormalCones[E],polytopes[Q].edgeNormalCones[F],&stackResource,&stackResource2);
						if(temp.isEmpty(&stackResource))
						{
							polytopes[P].relationTables[E].mark(Q,F,false);
							polytopes[Q].relationTables[F].mark(P,E,false);
						}
					}
			}
//			std::cerr<<toString();
	}
	template<typename typ>
	RelationTable PolytopeIntersectionData<typ>::getTableOfEdge(int polytopeIndex, int u, int v)const
	{
		return polytopes.at(polytopeIndex).relationTables[polytopes.at(polytopeIndex).edgeTable.lookup(u,v)];
	}
	template<typename typ>
	int PolytopeIntersectionData<typ>::bestPolytopeIndex(vector<int> const &used, HalfOpenCone<typ> C)const//using current cone
	{
			assert(0);
		int bestIndex=-1;
		int lowestNumber=0;
		for(int i=0;i<used.size();i++)
			if(!used[i])
			{
				int number=0;

	/*			if(0)
					for(auto & D:polytopes[i].vertexNormalCones)
				{
					auto E=D;
					HalfOpenCone<CircuitTableInt32> T(E.getInequalities(),E.getEquations(),Matrix<CircuitTableInt32>(E.getAmbientDimension(),0));
					auto I=intersection2(T,C);
					if(!I.isEmpty())number++;
				}
				else*/ if(1)
				for(auto & D:polytopes[i].edgeNormalCones)
				{
//						auto DD=D;
					auto E=D;//D.closure();
					if(0){//takes span of closure which might be easier to intersect with C. However, this gives a very bad heuristic.
						auto F=E.closure();
						E=HalfOpenCone<typ>(F.getOrthogonalComplement(),Matrix<typ>(F.getAmbientDimension(),0),Matrix<typ>(F.getAmbientDimension(),0));
					}
//						HalfOpenCone<CircuitTableInt32> T(E.getInequalities(),E.getEquations(),Matrix<CircuitTableInt32>(E.getAmbientDimension(),0));
					auto I=intersection2(E,C);
					if(!I.isEmpty())number++;
				}
				else
					for(auto & D:polytopes[i].edgeNormalCones)
					{
						auto DD=D;
						auto E=DD.closure();
						HalfOpenCone<typ> T(E.getInequalities(),E.getEquations(),Matrix<typ>(E.getAmbientDimension(),0));
						auto I=intersection2(T,C);
						if(!I.isEmpty())number+=I.getDimension();//*I.getDimension()*I.getDimension()*I.getDimension();
					}

				if(bestIndex==-1 || number<lowestNumber)
				{
					bestIndex=i;
					lowestNumber=number;
				}
			}
		assert(bestIndex!=-1);
		return bestIndex;
	}

	template<typename typ>
	int PolytopeIntersectionData<typ>::bestPolytopeIndex2(vector<int> const &used, RelationTable const &RT)const //Experiment 3
	{
			assert(0);
		int best;int bestIndex=-1;
		assert(used.size()==layout.startOffsets.size());
		for(int i=0;i<used.size();i++) ///polytope
			if(!used[i])
			{
				int c=0;
				for(int j=0;j<polytopes[i].relationTables.size();j++)//edge
					if(RT.lookUp(i,j))
				{
					//c+=1000; // It seems better to use the immediate RT and tie-break with the intersection
					auto &RT2=polytopes[i].relationTables[j];
					auto RT3=RT.intersection(RT2);
					bool hasZeroRow=false;
					for(int k=0;k<RT3.layout.startOffsets.size();k++)
						if((k!=i)&&(!used[k]))
						if(RT3.count(k)==0)
							hasZeroRow=true;
					if(!hasZeroRow)c++;//On cyclic examples some of the other 3 variants of this line perform better
				}

				if(bestIndex==-1 || c<best)
				{
					bestIndex=i;
					best=c;
				}
			}
		assert(bestIndex!=-1);
		return bestIndex;
	}

	template<typename typ>
	void print(vector<HalfOpenCone<typ> > &f, bool enableAssert=true)
	{
		int numberContainingOrigin=0;
		std::cout<<"************************************** start\n";
		int I=0;
		for(vector<HalfOpenCone<CircuitTableInt32> >::iterator i=f.begin();i!=f.end();i++,I++)
		{
			if(i->containsOrigin())numberContainingOrigin++;
			bool b=i->isEmpty();
			std::cout<<"-=-=-=-=-==-=-=-=-=-=-=-=- cone number "<<I<<":"<<f.size()<<"\n"
					<<i->toString()<<((b)?"EMPTY":"NOT EMPTY")
					<<i->lifted.dualCone.toString()
					<<"-=-=-=-=-==-=-=-=-=-=-=-=-\n";
		}
		std::cout<<"************************************** end\n";
		if(enableAssert)assert(numberContainingOrigin<=1);
	}

	template<typename typ>
	vector<HalfOpenCone<typ>> restrictedTropicalHypersurface(HalfOpenCone<typ> &cone, Matrix<typ> const &exponents, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
//		ResourceInvariant a(mr);
		MR *mr3=get_default_resource();
		// Problems to fix:
		// We dont just need to compute a spanning tree, but an entire orientation of the graph
		// For every vertex we should be able to find the outgoing edges and corresponding inequalities
		int cdim=cone.getDimension(mr);
//		std::cerr<<exponents;
//		std::cerr<<"Contains origin"<<cone.containsOrigin()<<"\n";
		vector<HalfOpenCone<typ> > ret;
		vector<HalfOpenCone<typ> > refinement;refinement.reserve(exponents.getHeight());
		pmrvector<bool> isFullDim(mr2);isFullDim.reserve(exponents.getHeight());
		pmrvector<int> comesFromVertexNumber(mr2);comesFromVertexNumber.reserve(exponents.getHeight());
		int numberOfFullDim=0;
		assert(!cone.isEmpty(mr));
//		std::cerr<<"TAKING CLOSURE\n";
		auto C=cone.closure(mr2,mr);
//		std::cerr<<"DONE TAKING CLOSURE\n";
//		std::cerr<<C.dualCone.resourceString()<<"\n";
		auto generatorsOfSpanOfC=C.getSpan(mr2,mr);

		C.getNumberOfFacets();//marks the facets - thereby simplifying computations later. Figure out if the facets can be figure out quicker - for example simplicial cones
		Matrix<int> adjacencyMatrix2(exponents.getHeight(),exponents.getHeight(),mr2);
		{
			pmrvector<int> toExplore(mr2);toExplore.reserve(exponents.getHeight());
			pmrvector<int> markedNumber(exponents.getHeight(),mr2);

			auto firstRelativeInteriorPoint=C.getRelativeInteriorPoint(mr2);
			auto firstOrder=combineOnTop(Matrix<typ>::rowVectorMatrix(firstRelativeInteriorPoint,mr),generatorsOfSpanOfC,mr2);
			{
				int number;
				int optimalVertexIndex=vertexIndex(exponents,firstOrder,number,mr);

				toExplore.push_back(optimalVertexIndex);
				markedNumber[optimalVertexIndex]=number;
			}
			while(!toExplore.empty())
			{
				int current=toExplore.back();
				toExplore.pop_back();
				Matrix<typ> empty(exponents.getWidth(),0,mr2);//HERE
				Matrix<typ> eq=C.getEquations(mr2,mr);
				Matrix<typ> ineq=C.getInequalities(mr2,mr);
#if 0
				Matrix<CircuitTableInt32> ineqAdditional(exponents.getHeight()-1,exponents.getWidth(),mr);
				for(int j=0;j<exponents.getHeight();j++)
					if(j!=current)
						ineqAdditional[j-(j>current)]=operatorMinus(exponents[current].toVector(mr),exponents[j].toVector(mr),mr2);

				{
					Matrix<CircuitTableInt32> ineqAdditional2(0,exponents.getWidth(),mr);ineqAdditional2.reserveRows(exponents.getHeight()-1);
					for(int i=0;i<ineqAdditional.getHeight();i++)
						if(!cone.knownToBeContainedInOpenHalfSpace(ineqAdditional[i].toVector(mr)))
							ineqAdditional2.appendRow(ineqAdditional[i].toVector(mr));
					ineqAdditional=ineqAdditional2;
				}
#else
				Matrix<typ> ineqAdditional(0,exponents.getWidth(),mr);ineqAdditional.reserveRows(exponents.getHeight()-1);
				for(int j=0;j<exponents.getHeight();j++)
					if(j!=current)
					{
						auto temp=operatorMinus(exponents[current].toVector(mr),exponents[j].toVector(mr),mr2);
						if(!cone.knownToBeContainedInOpenHalfSpace(temp))
//						if(!cone.lifted.isImplied(concatenation(temp,Vector<typ>(1,mr),mr)))
							ineqAdditional.appendRow(temp);
					}
#endif
				ineq=combineLeftRight(ineq,ineqAdditional.transposed(mr),mr2);
				Cone<typ> C2(ineq,eq,mr2,mr);
				auto te=HalfOpenCone<typ>(ineq/*ineq*/,empty,empty,mr2,mr);
				// A possible optimisation is to replace the line above with the line below.
//				auto te=HalfOpenCone<typ>(ineqAdditional.transposed(mr/*2*/),empty,empty,mr2,mr);
				auto C3=intersection2(te,cone,mr2,mr);
				auto facetNormals=C2.getFacetNormals(mr2,mr);
				for(int j=0;j<facetNormals.getWidth();j++)
					if(!compare(firstOrder,facetNormals.column(j,mr2),mr2))
				{
					Matrix<typ> theEquation=facetNormals.submatrixColumns(j,j+1,mr2);
					Cone<typ> facet(ineq,combineLeftRight(eq,theEquation,mr),mr2);//TO HERE
					auto relintF=facet.getRelativeInteriorPoint(mr2);
					if(C.containsInRelativeInterior(relintF))
					{
//						if(rigtigretning)
						{
							auto relativeInteriorPoint=C2.getRelativeInteriorPoint(mr2);//sometimes this can be reused
							auto generatorsOfSpanOfFacet=facet.getSpan(mr2,mr);
							auto order=combineOnTop(
									combineOnTop(Matrix<typ>::rowVectorMatrix(relintF,mr),
											Matrix<typ>::rowVectorMatrix(-relativeInteriorPoint,mr),mr),generatorsOfSpanOfFacet,mr2);
							int number;
							int optimalVertexIndex=vertexIndex(exponents,order,number,mr);
							adjacencyMatrix2[current][optimalVertexIndex]=1;
							adjacencyMatrix2[optimalVertexIndex][current]=1;
							if(!markedNumber[optimalVertexIndex])
							{
								toExplore.push_back(optimalVertexIndex);
								markedNumber[optimalVertexIndex]=number;
							}
						}
					}
				}
				refinement.push_back(C3);
				comesFromVertexNumber.push_back(current);
				isFullDim.push_back(markedNumber[current]>=2);
				numberOfFullDim++;
			}
		}


		assert(refinement.size());
		int ncones=refinement.size();


		Matrix<int> adjacencyMatrix(ncones,ncones,mr2);
		for(int i=0;i<ncones;i++)
			for(int j=0;j<ncones;j++)
				adjacencyMatrix[i][j]=adjacencyMatrix2[comesFromVertexNumber[i]][comesFromVertexNumber[j]];

		//			std::cerr<<"ADJ"<<adjacencyMatrix.toString();
			int root=0;
			int bestScore=0;
			for(int i=0;i<ncones;i++)
			{
				int sum=0;
				for(int j=0;j<ncones;j++)if(adjacencyMatrix[i][j]&&!isFullDim[j])sum++;
				int score=100*isFullDim[i]+sum;
//				std::cerr<<score<<"\n";
				if(score>bestScore)
				{
					root=i;
					bestScore=score;
				}
			}

			auto w=refinement[root].getRelativeInteriorPoint(mr2);


			auto orderW=combineOnTop(Matrix<typ>::rowVectorMatrix(w,mr),generatorsOfSpanOfC,mr2);

			// The following seems to be dead data/codes
#if 0
			pmrvector<int> order(ncones,mr2);
			for(int i=0;i<ncones;i++)order[i]=-1;
			order[root]=0;
			int numberOfOrderedVertices=1;
			while(numberOfOrderedVertices<ncones)
			{
				bool found=false;
				bool foundgood=false;
				int I,J;
				for(int i=0;i<ncones;i++)
					if(order[i]!=-1)
						for(int j=0;j<ncones;j++)
							if(order[j]==-1)
//								if(adjacencyMatrix[i][j])
									if(adjacencyMatrix[i][j]&&(!found||!foundgood))
								{
									found=true;
									I=i;
									J=j;
									i=ncones;j=ncones;//break
									//foundgood=isFullDim[j]&&isFullDim[i];
								}
				assert(found);
				if(found)
				{
					order[J]=numberOfOrderedVertices;
					numberOfOrderedVertices++;
				}
			}
#endif

//			a.assertIsInvariant();

			//			for(int i=0;i<order.size();i++)std::cerr<<order[i];
//			std::cerr<<"\n";
//			std::cerr<<"REFINEMENT\n";
			//print(refinement,false);
			auto zero=Vector<typ>(exponents.getWidth(),mr2);
			for(int i=0;i<refinement.size();i++)
			{
//				ResourceInvariant a(mr);
				Matrix<typ> strict(0,exponents.getWidth(),mr2);strict.reserveRows(refinement.size()-1);
				Matrix<typ> nonstrict(0,exponents.getWidth(),mr2);nonstrict.reserveRows(exponents.getHeight()-1+refinement.size());
//				{//				auto mrTemp=ResourceWrapper(mr);
				for(int j=0;j<exponents.getHeight();j++)
					if(j!=comesFromVertexNumber[i])//this does not look like the general way of making COpen the interior
						nonstrict.appendRow(operatorMinus(exponents[comesFromVertexNumber[i]].toVector(mr),exponents[j].toVector(mr),mr));
//				}
				//{				auto mrTemp=ResourceWrapper(mr);
				for(int j=0;j<refinement.size();j++)
					if(adjacencyMatrix[i][j])
					{
						assert(i!=j);
						auto W=operatorMinus(exponents[comesFromVertexNumber[i]].toVector(mr),exponents[comesFromVertexNumber[j]].toVector(mr),mr);
						if(isPerpendicular(W,orderW))continue;
//						if(compare(zero,W))
						if(!vertexCompare(W,zero/*W-W*/,orderW))
							nonstrict.appendRow(operatorMinus(exponents[comesFromVertexNumber[i]].toVector(mr),exponents[comesFromVertexNumber[j]].toVector(mr),mr));
						else
							strict.appendRow(operatorMinus(exponents[comesFromVertexNumber[i]].toVector(mr),exponents[comesFromVertexNumber[j]].toVector(mr),mr));

					}
				//}
				if(1){
					auto A=HalfOpenCone<typ>(nonstrict.transposed(mr),Matrix<typ>(0,exponents.getWidth(),mr).transposed(mr),strict.transposed(mr),mr2,mr);//HERE TO HERE
					A=intersection2(A,cone,mr2,mr);
					Matrix<typ> HiPri(0,exponents.getWidth(),mr2);HiPri.reserveRows(nonstrict.getHeight());
					Matrix<typ> LoPri(0,exponents.getWidth(),mr2);LoPri.reserveRows(nonstrict.getHeight());
					for(int j=0;j<nonstrict.getHeight();j++)
					{
						auto v=nonstrict[j].toVector(mr);
						if(A.validInequalityDefinesExistingFacet(v,mr))
							HiPri.appendRow(v);
						else
							LoPri.appendRow(v);
					}
					nonstrict=combineOnTop(LoPri,HiPri,mr2);
				}
if(1)				{
//					std::cerr<<"Old:"<<nonstrict<<"\n";
					HalfOpenCone<typ> temp(nonstrict.transposed(mr),Matrix<typ>(0,exponents.getWidth(),mr).transposed(mr),strict.transposed(mr),mr2);
					auto temp2=temp.facetDefining();//HERE
//					int nfacets=0;
//					for(int i=0;i<temp.size();i++)if(temp[i])nfacets++;
					Matrix<typ> nonstrict2(0,exponents.getWidth(),mr2);
					for(int i=0;i</*temp2.size()*/nonstrict.getHeight();i++)
						if(temp2[i])
							nonstrict2.appendRow(nonstrict[i].toVector(mr));
					nonstrict=nonstrict2;
//					std::cerr<<"New:"<<nonstrict<<"\n";
				}

				if(isFullDim[i])
				{
					auto A=HalfOpenCone<typ>(nonstrict.transposed(mr),Matrix<typ>(0,nonstrict.getWidth(),mr).transposed(mr),strict.transposed(mr),mr2,mr);

					ret.reserve(ret.size()+1);
					ret.emplace_back(std::move(intersection2(A,cone,mr,mr2)),mr);//FAILS IF ENABLED
					if(ret.back().isEmpty(mr))ret.pop_back();
				}
				else
				{
					pmrvector<HalfOpenCone<typ> > temp(mr2);temp.reserve(nonstrict.getHeight());
					codimension1ConesB(nonstrict.transposed(mr),strict.transposed(mr),temp,mr2,mr);//This causes a a lot of memory ussage

					ret.reserve(ret.size()+temp.size());
					for(int j=0;j!=temp.size();j++)
					{
//						ret.reserve(ret.size()+1);
						ret.emplace_back(std::move(intersection2(temp[j],cone,mr2,mr)),mr);//FAILS IF ENABLED
						if(ret.back().isEmpty(mr))ret.pop_back();
					}
				}
			}

			return ret;
		}
// ---------------
// ProgressCounter
// ---------------
	void ProgressCounter::newIterator(int c){bounds.push_back(c);counter.push_back(0);}
	void ProgressCounter::iterate()
	{
		counter.back()++;
		if(counter.back()==bounds.back())
		{
			counter.pop_back();
			bounds.pop_back();
		}
	}
	string ProgressCounter::toString()const
	{
		stringstream s;
		s<<"Bounds:";
		for(auto &b:bounds)s<<std::setw(3)<<b;
		s<<"\n";
		s<<"Counts:";
		for(auto &b:counter)s<<std::setw(3)<<b;
		s<<"\n";
		return s.str();
	}


	template<typename typ>
	class CommonRefinementTraverser: public Traverser
	{
		StackResource stackResource;
		StackResource stackResource2;
		vector<HalfOpenCone<typ>> coneStack;
		vector<vector<int>> usedStack;
		vector<RelationTable> relationTableStack;
		vector<vector<HalfOpenCone<typ>>> coneOptionsStack;
		vector<int> chosenStack;

		vector<Matrix<typ> > exponents;
		std::function<void (HalfOpenCone<typ>&)> collector;
		std::function<bool (HalfOpenCone<typ>&)> filter;
		PolytopeIntersectionData<typ> &data;
//		vector<int> const &edgeCounts;
		bool filterOK;
		int depth;
		CommonStatistics *statistics;
		public:
		void setupOptionStack()
		{
//			std::cerr<<"Setup\n";
			if(depth!=data.polytopes.size() && filterOK)
			{
				if(statistics)statistics->numberOfIntermediateVertices++;
				int chosen=relationTableStack.back().bestPolytopeIndex(usedStack.back());
				chosenStack.push_back(chosen);
				usedStack.push_back(usedStack.back());usedStack.back()[chosen]=1;

				if(1)
				{
				auto temp=restrictedTropicalHypersurface((coneStack.back()),exponents[chosen],&stackResource,&stackResource2/*,mr,mr2 WHAT DO WE DO ABOUT ALLOCATORS?*/);
				vector<HalfOpenCone<typ>> temp2(temp.size());
				temp2.resize(distance(temp2.begin(),copy_if(temp.begin(),temp.end(),temp2.begin(),[](HalfOpenCone<typ> &c)->bool{return !c.isEmpty();})));
				coneOptionsStack.push_back(temp2);//can one of these be empty? Should we filter?
				}
				else
				coneOptionsStack.push_back(
						restrictedTropicalHypersurface(
								coneStack.back(),
								exponents[chosen],
								&stackResource,
								&stackResource2
							)
				);
			}
			else
			{
				coneOptionsStack.push_back(vector<HalfOpenCone<typ>>());
				usedStack.push_back(usedStack.back());// we need to push something, so that we have something to pop
				chosenStack.push_back(0);
			}
//			std::cerr<<"Done setup\n";
		}
		CommonRefinementTraverser(
				HalfOpenCone<typ> cone,
				vector<int> &used,
				vector<Matrix<typ> > const &exponents_,
				std::function<void (HalfOpenCone<typ>&)> collector_,
//				vector<int> const &edgeCounts_,
				PolytopeIntersectionData<typ> /*const*/ &data_,
				RelationTable const &relationTable,
				std::function<bool (HalfOpenCone<typ>&)> filter_,
				CommonStatistics *statistics_):
					coneStack{cone},
					usedStack{used},
					relationTableStack{relationTable},
					collector(collector_),
					exponents(exponents_),
					data(data_),
//					edgeCounts(edgeCounts_),
					filter(filter_),
					depth(0),
					stackResource(10160000),
					stackResource2(10160000),
					statistics(statistics_)
					{
						filterOK=filter(coneStack.back());
						setupOptionStack();
					}

		virtual int  getEdgeCountNext( void )
		{
			return coneOptionsStack.back().size();
		}

	  virtual int  moveToNext( int   index,
	                           bool  collect_info )
	  {
//		  std::cerr<<"Next\n";
		  coneStack.push_back(coneOptionsStack.back()[index]);


		  int chosen=chosenStack.back();

			 Vector<typ> v=coneStack.back().getRelativeInteriorPoint(&stackResource,&stackResource2);
			 set<int> indices;
			 for(int i=0;i<data.polytopes[chosen].edgeNormalCones.size();i++)
			 {
				 auto &T=data.polytopes[chosen].edgeNormalCones[i];
				 if(T.closureOfNonEmptyConeContains(v,&stackResource))indices.insert(i);
			 }
//				 assert(indices.size()>0); THIS CAN FAIL - what does that mean?

			 if(indices.size()==1)
				 relationTableStack.push_back(relationTableStack.back().intersection(
						 data.polytopes[chosen].relationTables[*indices.begin()]
						 ));
			 else
				 relationTableStack.push_back(relationTableStack.back());




		  depth++;
		  filterOK=filter(coneStack.back());
		  setupOptionStack();
		  return 0;
	  }

	  virtual void  moveToPrev( int  index )
	  {
//		  std::cerr<<"Prev\n";
		  coneStack.pop_back();
		  usedStack.pop_back();
		  chosenStack.pop_back();
		  relationTableStack.pop_back();
		  coneOptionsStack.pop_back();
		  depth--;
		  filterOK=true;
	  }

	  // This function will be called once for every state in a traversal.
//	  int collectCounter = 0;
	  virtual void  collectInfo( void ){
//	    if (++collectCounter % 1000 == 0)
//	      std::cout << "." << std::flush;
	    if(depth==usedStack.back().size())
	      if(depth!=0 || !coneStack.back().isEmpty())
	    	  if(filterOK)
	    	  collector(coneStack.back());
	  }

	  // Function for printing the state to cout for debug purposes.
	  virtual void  printState( void ){}
	};


	template<typename typ>
	void commonRefinement(
			HalfOpenCone<typ> cone,
			vector<int> &used,
			vector<Matrix<typ> > const &exponents,
			std::vector<HalfOpenCone<typ> > &collector,
			int &nonEmptyIntersections,
//			vector<int> const *edgeCounts,
			PolytopeIntersectionData<typ> /*const*/ &data,
			RelationTable const &relationTable,
			ProgressCounter &progress,
			CommonStatistics *commonStatistics,
			int j,
			std::function<bool (HalfOpenCone<typ>&)> filter)
		{
		// std::cerr<<"NumberOfPivotingStep at start:"<<numberOfPivotingSteps<<"\n";
		mutex m;
		//		std::function<void (HalfOpenCone<typ>&)> collector_
		vector<CommonRefinementTraverser<typ>> V;
		V.reserve(j);
		for(int i=0;i<j;i++)V.emplace_back(cone,used,exponents,
				//collector_
				[&collector,&m](HalfOpenCone<typ> const& c)mutable->void{lock_guard<std::mutex> lck{m};collector.push_back(c);}
				/*,*edgeCounts*/,data,relationTable,filter,commonStatistics);
/*		CommonRefinementTraverser<typ> T(cone,used,exponents,
				//collector_
				[&collector,&m](HalfOpenCone<typ> const& c)mutable->void{lock_guard<std::mutex> lck{m};collector.push_back(c);}
				,*edgeCounts,data,relationTable,filter);*/
//		traverse_simple(&T);
		if(j==1)
			traverse_simple(&V[0]);
		else
		{
			assert(j>0);
			vector<Traverser*> V2;
			for(int i=0;i<j;i++)V2.push_back(&V[i]);
			traverse_threaded(&V2[0],j,10,true);
		}

#if 0
		//does not return values on stack, so it makes sense to initialise mr and mr2 at this level
		if(!filter(cone))return;
				 MR *mr=&stackResource;
				 MR *mr2=&stackResource2;

		 if(!cone.isEmpty(mr))
		 {
			 nonEmptyIntersections++;
			 int sum=0;for(auto v:used)sum+=v;
			 if(sum==exponents.size())
			 {
				 collector.push_back(/*std::move*/(cone));
			 }
		 else
		 {
#if 0
			 auto tem=cone.closure(mr,mr2);
//				vector<int> dontChoose=used;
				auto w=tem.getRelativeInteriorPoint(mr);
#endif
				int chosen=0;
#if 0
				if(0)
				 for(chosen=0;chosen<exponents.size();chosen++){if(!used[chosen])break;}
			 else
			 {
				 int bestValue=10000000;
				 int best=0;
				 for(int i=0;i<exponents.size();i++)
					 if(!used[i])
				 {
					 auto A=tem.getInequalities(mr2,mr).transposed(mr);//HERE
					 auto B=tem.getEquations(mr2,mr).transposed(mr);
					 for(int j=1;j<exponents[i].getHeight();j++)
						 A.appendRow(exponents[i][0].toVector(mr2)-exponents[i][j].toVector(mr2));
					 Matrix<typ> empty(A.getWidth(),0);
					 Matrix<typ> e=combineOnTop(A,B,mr).transposed(mr2);//TO HERE
					 Cone<typ> s(empty,e,mr,mr2);
					 int D=s.getDimension();
					 int d=D;

if(1)					 {
						 int max=-10000000;
						 int min=10000000;
						 int max2=max;//second largest
						 int count=0;
						 for(int j=0;j<exponents[i].getHeight();j++)
						 {
							 int prod=dot(exponents[i][j],w).toInt64();
							 if(prod>max){max2=max;max=prod;count=1;}
							 else if(prod==max){count++;max2=max;}
							 else if(prod>max2){max2=prod;}
							 if(prod<min){min=prod;}
						 }
						 d=100*count;//+D;
						 if(max>min)
							 if(count<2)
						 {
							 d=100*count+100*(max2-min)/(max-min);
						 }

					 }
					 if(d<bestValue)
					 {
						 bestValue=d;
						 best=i;
						 dontChoose=used;for(int j=0;j<i;j++)dontChoose[j]=true;
						 dontChoose[i]=false;
					 }
					 else if(d==bestValue)
					 {
						 if((*edgeCounts)[i]>(*edgeCounts)[best])
						 {
							 bestValue=d;
							 best=i;
							 dontChoose=used;for(int j=0;j<i;j++)dontChoose[j]=true;
							 dontChoose[i]=false;
						 }
						 if((*edgeCounts)[i]==(*edgeCounts)[best])
						 {
							 dontChoose[i]=false;
						 }
						 else
							 dontChoose[i]=true;
					 }
					 else
						 dontChoose[i]=true;
				 }
				 chosen=best;
			 }
#endif
			 chosen=relationTable.bestPolytopeIndex(used);       //Experiment 1 in email to Jeff?
	//		 chosen=data.bestPolytopeIndex(used,cone); 			 //Experiment 2 in email to Jeff?
	//		 chosen=data.bestPolytopeIndex2(used,relationTable); //Experiment 3 in email to Jeff?

//			 chosen=relationTable.bestPolytopeIndex(dontChoose);


			 vector<HalfOpenCone<typ>>f=restrictedTropicalHypersurface((cone),exponents[chosen],mr,mr2);

#define NUMBER_OF_ITERATIONS_TO_EXECUTE_BEFORE_PRINT 1000
			 {static int i; i++;if(i==NUMBER_OF_ITERATIONS_TO_EXECUTE_BEFORE_PRINT){std::cerr<<progress.toString()<<"\n";i=0;}}
			 progress.newIterator(f.size());
			 for(int i=0;i<f.size();i++)
			 {
				 if(!f[i].isEmpty(mr))//How is it possible for one of these cones to be empty???
			 {
				 used[chosen]=1;

				 Vector<typ> v=f[i].getRelativeInteriorPoint(mr,mr2);
				 set<int> indices;
				 for(int i=0;i<data.polytopes[chosen].edgeNormalCones.size();i++)
				 {
					 auto &T=data.polytopes[chosen].edgeNormalCones[i];
					 if(T.closureOfNonEmptyConeContains(v,mr))indices.insert(i);
				 }
//				 assert(indices.size()>0); THIS CAN FAIL - what does that mean?

				 if(indices.size()==1)
				 {
					 commonRefinement(f[i],used,exponents,collector,nonEmptyIntersections,edgeCounts,data,
							 relationTable.intersection(
									 data.polytopes[chosen].relationTables[*indices.begin()]
									 ),progress,filter
					 );
				 }
				 else
					 commonRefinement(f[i],used,exponents,collector,nonEmptyIntersections,edgeCounts,data,relationTable,progress,filter);

				 used[chosen]=0;
			 }
				 progress.iterate();
			 }
		 }
			 }
#endif
		// std::cerr<<"NumberOfPivotingSteps at end:"<<numberOfPivotingSteps<<"\n";
		}
/*	template<typename typ>
	 ZMatrix toZMatrix(Matrix<typ> const &m)
	 {
		 ZMatrix ret(m.getHeight(),m.getWidth());
		 for(int i=0;i<m.getHeight();i++)
			 for(int j=0;j<m.getWidth();j++)
				 ret[i][j]=m[i][j].v;
		 return ret;
	 }*/
	template<typename typ>
	 string toString(std::vector<HalfOpenCone<typ> > &collector, int ambientDimension)
	 {
		 ZFan F(ambientDimension);
		 for(auto c:collector)
		 {
			 auto C=c.closure();
			 F.insert(ZCone(toZMatrix(C.getFacetNormals().transposed()),toZMatrix(C.getOrthogonalComplement().transposed())));
		 }
		 return F.toString();
	 }
};

#include "gfanlib_circuittableint.h"
#include "gfanlib_circuittableinteger.h"
template class gfan::PolytopeIntersectionData<gfan::CircuitTableInt32>;
template class gfan::PolytopeIntersectionData<gfan::CircuitTableInt64>;
template class gfan::PolytopeIntersectionData<gfan::CircuitTableInteger>;
template void gfan::commonRefinement<gfan::CircuitTableInt32>(HalfOpenCone<CircuitTableInt32> cone, vector<int> &used, vector<Matrix<CircuitTableInt32> > const &exponents, std::vector<HalfOpenCone<CircuitTableInt32> > &collector,int &nonEmptyIntersections/*,vector<int> const *edgeCounts*/, PolytopeIntersectionData<CircuitTableInt32> /*const*/ &data, RelationTable const &relationTable, ProgressCounter &progress, CommonStatistics *statistics, int j, std::function<bool (HalfOpenCone<CircuitTableInt32>&)> filter);
template void gfan::commonRefinement<gfan::CircuitTableInt64>(HalfOpenCone<CircuitTableInt64> cone, vector<int> &used, vector<Matrix<CircuitTableInt64> > const &exponents, std::vector<HalfOpenCone<CircuitTableInt64> > &collector,int &nonEmptyIntersections/*,vector<int> const *edgeCounts*/, PolytopeIntersectionData<CircuitTableInt64> /*const*/ &data, RelationTable const &relationTable, ProgressCounter &progress, CommonStatistics *statistics, int j, std::function<bool (HalfOpenCone<CircuitTableInt64>&)> filter);
template void gfan::commonRefinement<gfan::CircuitTableInteger>(HalfOpenCone<CircuitTableInteger> cone, vector<int> &used, vector<Matrix<CircuitTableInteger> > const &exponents, std::vector<HalfOpenCone<CircuitTableInteger> > &collector,int &nonEmptyIntersections/*,vector<int> const *edgeCounts*/, PolytopeIntersectionData<CircuitTableInteger> /*const*/ &data, RelationTable const &relationTable, ProgressCounter &progress, CommonStatistics *statistics, int j, std::function<bool (HalfOpenCone<CircuitTableInteger>&)> filter);

namespace gfan{
	int64 numberOfPivotingSteps;
};
//--------------------------------------------------------------- remove after more lambda expression experience
#include <functional>

void testabc(std::function<bool (int)> func)
{
	for(int i=0;i<10;i++)std::cerr<<func(i)?"True":"False";
}
void testabce()
{
	int a[10]={1,2,4,5,7,0,2,3,5,3};
	testabc([a](int z)->bool{return a[z]%3;});
}
