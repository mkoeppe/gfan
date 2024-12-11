/*
 * gfanlib_hypersurfaceintersection.h
 *
 *  Created on: 12 Mar 2020
 *      Author: anders
 */

#ifndef SRC_GFANLIB_HYPERSURFACEINTERSECTION_H_
#define SRC_GFANLIB_HYPERSURFACEINTERSECTION_H_


#include "gfanlib_tableau.h"
#include <functional>

namespace gfan{


class VertexPairToEdgeTable{
	Matrix<int> table;
	int nextFree;
public:
	VertexPairToEdgeTable(int numberOfVertices);
	int allocate(int i, int j);
	int lookup(int i, int j)const;
	string toString()const;
};

class RelationTableLayout{
public:
	vector<int> startOffsets;
	vector<int> endOffsets;
	int nextFree;
public:
	RelationTableLayout();
	void allocate(int size);
	int numberOfUsedBits()const;
};

class RelationTable{
public:
	vector<bool> data;
	RelationTableLayout const &layout;
public:
	RelationTable(RelationTableLayout const &layout_);
	RelationTable intersection(RelationTable const &b)const;
	void mark(int polytopeIndex, int edgeIndex, bool b);
	bool lookUp(int polytopeIndex, int edgeIndex)const;
	int count(int polytopeIndex)const;
	int bestPolytopeIndex(vector<int> const &used)const;
};

template<typename typ>
class PolytopeData{
public:
	Matrix<typ> generators;
	vector<HalfOpenCone<typ>> edgeNormalCones;
	VertexPairToEdgeTable edgeTable;
	vector<pair<int,int>> edges;
	vector<int> vertices;
	vector<RelationTable> relationTables;//one table for each edge
public:
	PolytopeData(Matrix<typ> const &generators_);
	void initializeRelationTables(RelationTableLayout const &layout);
	string toString(RelationTableLayout const &layout)const;
};

template<typename typ>
class PolytopeIntersectionData{
public:
	vector<PolytopeData<typ>> polytopes;
	RelationTableLayout layout;
	std::string toString();
	PolytopeIntersectionData(vector<Matrix<typ>> const &polytopes_);
	RelationTable getTableOfEdge(int polytopeIndex, int u, int v)const;
	int bestPolytopeIndex(vector<int> const &used, HalfOpenCone<typ> C)const;//using current cone
	int bestPolytopeIndex2(vector<int> const &used, RelationTable const &RT)const; //Experiment 3
};

class CommonStatistics{
public:
	std::atomic<int> numberOfIntermediateVertices;
	CommonStatistics():
		numberOfIntermediateVertices{0}
	{
	}
};

class ProgressCounter{
	vector<int> bounds;
	vector<int> counter;
public:
	void newIterator(int c);
	void iterate();
	string toString()const;
};

template<typename typ>
void commonRefinement(
		HalfOpenCone<typ> cone,
		vector<int> &used,
		vector<Matrix<typ> > const &exponents,
		std::vector<HalfOpenCone<typ> > &collector,
		int &nonEmptyIntersections,
//		vector<int> const *edgeCounts,
		PolytopeIntersectionData<typ> /*const*/ &data,
		RelationTable const &relationTable,
		ProgressCounter &progress,
		CommonStatistics *statistics=0,
		int j=8,
		std::function<bool (HalfOpenCone<typ>&)> filter=[](HalfOpenCone<typ>& c) -> bool { return true; });

template<typename typ>
Vector<int> fvector(std::vector<HalfOpenCone<typ> > &f)
{
	assert(f.size());
	Vector<int> s(f.begin()->lifted.getAmbientDimension());
	for(auto &c:f)
		s=s+c.fVector();
	return s;
}
};
#endif /* SRC_GFANLIB_HYPERSURFACEINTERSECTION_H_ */
