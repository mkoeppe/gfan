/*
 * gfanlib_tableau.h
 *
 *  Created on: May 5, 2017
 *      Author: anders
 */

#ifndef GFANLIB_TABLEAU_H_
#define GFANLIB_TABLEAU_H_

#include <algorithm>
#include <bitset>
#include "gfanlib_matrix.h"
#include "gfanlib_circuittableint.h"
#include "gfanlib_frequencytable.h"
#include "gfanlib_zfan.h"
#include "gfanlib_polymakefile.h"
#include "vektor.h"


#define NOEXCEPT_
/*
 * Simplex algorithm without objective function.
 * The original matrix does not have to be stored.
 * At any time a B denotes a submatrix of columns
 * indexed by basisIndices. This matrix also does
 * not have to be stored. The adjoint matrix adj(B)
 * is also unknown, but its product with the original
 * matrix is kept updated as the columns of the basis
 * B are exchanged. Schrijver version of Bland's rule is used.
 *
 *
 * My understanding of the code on July the 14th 2017:
 * 1) The interesting class is TableauSolver and its coneInfo() method, as that can basically be extended into the functionality of gfans PolyhedralCone objects.
 * 2) The tableau represents a cone given by generators.
 * 3) The generators are passed as a matrix to the constructor which will extend it by appending a row of ones and prepending an identity matrix.
 * 4) The purpose of adding the identity matrix is to get started with the pivoting. We may immediately use 0..d-1 as a basis.
 * 5) The purpose of adding the ones vector becomes appearant when studying the code for determining the lineality space. For other parts of the code, this row seems to be ignored.
 * 6) Each row is marked by a flag ignoreRows. If the flag is set for a row, then that row contains as its first d coordinates a combination of the original rows giving 0. When the row has been marked to be ignored, it will nolonger change. Combined the ignored rows form a basis of the orthogonal complement of the generated cone.
 * 7) Each column is marked by a flag inLineality indicating whether both directions of the vector are in the generated cone. After these flags have been determined, we know a generating set for the lineality space of the cone.
 * 8) Each column is marked by a flag nonExtreme, indicating whether the column is extreme. If colums generate the same ray modulo the lineality space, the intend is to mark only one of the generators as extreme.
 * 9) As the inLineality flags are discovered, columns are made permanent members of the basis. These permanent members will correspond to a basis of the lineality space.
 *
 * States to support
 * 0: none?
 * 1: basis of implied equations of polar cone <-> basis for the lineality space of the cone
 * 2: redundant inequalites of polar has been removed <-> extreme generators as described in 8) above
 * 3: canonical form
 *
 * We would like to support the following functions for the polar cone:
 * getDimension in state 1 = dimension of lineality space of primal
 * dimensionOfLinealitySpace in state 0.5 = dimension of orthogonal complement
 * getHalfSpaces = extreme rays of primal
 * getEquations = makes not so much sense
 * getImpliedEquations = marked by inLineality (or found as a permanent member of the basis)
 *
 * Maybe we should allow some generators to point in both directions. That would need some way of finding a basis among those rays and some how take care of this information in the initial steps.
 *
 * It seems that the info function should be split into several parts
 * Part A: Initialise without any computation
 * Part B: Figure out the orthogonal complement and the dimension
 * Part C: Find the dimension of the lineality space and mark a basis for it.
 * Part D: Find the extreme rays - only one for each direction
 *
 *
 * The confusing column d-1 and row d-1.
 * -------------------------------------
 * Most computations work in dimension d-1 -ignoring the last row of combinedMatrix.
 * When the row is ignored column d-1 is marked as the last basis element to make the basis have size d.
 * In this mode the fact that d-1 is listed in the basis should be ignored.
 * The last row is used only in the case where positive dependencies are sought to get genererators of the lineality space.
 *
 * Column d-1 is used for two things:
 * 1) checking if a vector is in the cone spanned by the columns given in the constructor.
 * 2) finding the lineality space
 *
 * The problem with having these two modes of operation is that the common multiplied determinants can come out of sync in different
 * parts of the matrix. Therefore, when entering the mode for finding lineality space, the entire last row should be set. On the other,
 * it seems that it is important that the last row is not updated by accident when not intended in order to avoid the possibility of making
 * illegal divisions in the last row.
 */

/*
 * Thoughts about finding relative interior points.
 *
 * Problem 1. Given a set of vectors, figure out if there is a facet of the generated cone with a normal vector having negative last coordinate and compute that facet normal.
 * Problem 2. Given a set of outer normals, find a point in the cone with with positive last coordinate.
 * Problem 3. Given a GeneratedCone, find a generic implied inequality i.e. a relative interior point of the dual cone.
 * Problem 4. Given a Cone, find a relative interior point.
 *
 * Problem 1 can be solved by assigning (0,0,-1) to d-1st column and move towards it
 * Problem 1 and Problem 2 are dual.
 * Problem 3 and Problem 4 are dual.
 * Problem 4 can be solved after finding the implied equations by adding a new coordinate being 1 facet inequalities and 0 for other. Then asking for a problem of type 2 to be solved.
 */


/*
 * Thoughts on deleting rows/columns
 *
 *           LLL
 *  B B   B   B
 * [0A0AEG0GG?d?] I
 * [0AdAE0000000] I
 * [0A0AEGdGG?0?]
 * [dA0AEG0GG?0?]
 * [            ]
 *
 * When columns are discovered to be in the lineality space - if they are not part of the basis, they should be forgotten (unless they are part of the first square submatrix).
 * This does not happen at the moment. An efficient implementation will permute the columns to make sure that undeleted columns appear together.
 *
 * When a row is marked to be ignored, it happens for one of two reasons:
 *  * The first entries of the row are in the orthogonal complement of the cone (last part are all zeros)
 *    In this case the row stays in the basis forever. It is unnecessary to create new zeros in this row as the remaining basis moves around - essentially because it does not matter what values column vectors attain at this row.
 *  * The first entries of the row together with an additional entry encode a positive combination giving zero - thereby
 *    giving a vector in the lineality space.
 *
 *  CLAIM:
 *   After a row is set to be ignored, its ignore flag will never be cleared and we only need that row for one purpose,
 *   namely the first d entries will be a generator for the orthogonal complement of the cone (if the row was ignored because of
 *   orthogonal complement). In particular this row will violate the requirement that the submatrix of basis columns is a scaled identity matrix.
 */
namespace gfan{

 class AtReturn
 {
 	std::string s;
 public:
 	AtReturn(const char *s_):s(s_){}
 	~AtReturn(){std::cerr<<s;}
 };

extern int64 numberOfPivotingSteps;

//template<class typ>
//class Matrix<typ,typename std::enable_if<has_member(typ,POD2),CircuitTableInt32POD>::type >{};

//constexpr bool hasPod(CircuitTableInt32 *){return true;}

template<typename mvtyp>
mvtyp gcd(Vector<mvtyp> const &v)
{
	assert(!v.isZero());
	int i;
	for(i=0;i<v.size();i++)if(!v[i].isZero())break;
	mvtyp ret=abs(v[i]);
	i++;
	while(i<v.size())ret=gcd(ret,v[i++]);
	return ret;
}

template<typename mvtyp>
Vector<mvtyp> normalize(Vector<mvtyp> const &v)
{
	Vector<mvtyp> ret=v;
	mvtyp d=gcd(ret);
	assert(!d.isZero());
	for(int i=0;i<ret.size();i++)ret[i]=ret[i]/d;
	return ret;
}

template<typename T>
string toString2(pmrvector<T> const&v)
{
	stringstream s;
	s<<"(";
	for(int i=0;i<v.size();i++)s<<(i?",":"")<<v[i];
	s<<")";
	return s.str();
}

template<typename T>
string toString2(vector<T> const&v)
{
	stringstream s;
	s<<"(";
	for(int i=0;i<v.size();i++)s<<(i?",":"")<<v[i];
	s<<")";
	return s.str();
}

static ZVector toZVector(gfan::Vector<int> const &m)
{
	ZVector ret(m.size());
	for(int i=0;i<m.size();i++)
		ret[i]=Integer(m[i]);
	return ret;
}

template<typename mvtyp> ZMatrix toZMatrix(Matrix<mvtyp> const &m)
{
	ZMatrix ret(m.getHeight(),m.getWidth());
	for(int i=0;i<m.getHeight();i++)
		for(int j=0;j<m.getWidth();j++)
			ret[i][j]=Integer(m[i][j].toInt64());
	return ret;
}
template<typename mvtyp> Matrix<mvtyp> fromZMatrix(ZMatrix const &m)
{
	Matrix<mvtyp> ret(m.getHeight(),m.getWidth());
	for(int i=0;i<m.getHeight();i++)
		for(int j=0;j<m.getWidth();j++)
			ret[i][j]=mvtyp(m[i][j].toInt());
	return ret;
}

static_assert(has_member(CircuitTableInt32,POD2),"");

template <class mvtyp> class Tableau{
	private:
	public:
	pmrvector<bool> isColumnInvalid;
		/*
		 * We disallow the columns to be made invalid as they usually store the transformation or target vector.
		 * Ideally the columns should be shuffled to the end of the matrix to get more efficient code. That will
		 * happen in a later version. For now the entries of the columns are just set to zero.
		 * There are several purposes of setting entries to zero:
		 *    freeing memory
		 *    efficiency of dealing with zeros rather than arbitrary precision ints
		 *    one cannot just stop some of the updates of a column without setting the column to zero as that would eventually lead to non-exact division and overflows
		 */
		void makeColumnInvalid(int j)
		{
			if(!isColumnInvalid[j])
			if(j>=combinedMatrix.getHeight())
			{
				isColumnInvalid[j]=true;
				for(int i=0;i<combinedMatrix.getHeight();i++)
					combinedMatrix[i][j]=mvtyp(0);
			}
//			computeRowBounds();//!!!!! It is maybe not necessary to recompute the bounds - but useful for debugging
		}
		bool CHECKCOL(int j)const{
			return true;
			if(isColumnInvalid[j])
			{
				std::cerr<<(const_cast<Tableau*>(this))->toString();
				std::cerr<<"INVALID INDEX:"<<j<<"\n";
				std::cerr<<"invalid:"<<toString2(isColumnInvalid)<<"\n";
			}
			return !isColumnInvalid[j];
		}
		Matrix<mvtyp> combinedMatrix; // Adj(B)*M
		Vector<mvtyp> negativeRowBounds;
		/* Initially no rows are ignored.
		 * As the rows are marked as ignored, the cannot be changed. The only valid use of a row is to read off the first entries to obtain vectors perpendicular to the cone.
		 *
		 *  In the cone class there are two reasons that a row may be ignored:
		 *    dimension is reduced when finding orthogonal complement
		 *    and when finding lineality
		 */
		pmrvector<bool> ignoreRows;
	public:
		int getWidth()const{return combinedMatrix.getWidth();}
		int getHeight()const{return combinedMatrix.getHeight();}
		int signAt(int i, int j)const{return combinedMatrix.UNCHECKEDACCESS(i,j).sign();/*return combinedMatrix[i][j].sign();*/}
		void assignEntry(int i, int j, mvtyp const &v){combinedMatrix[i][j]=v;}
		pmrvector<int> basisIndices;
		pmrvector<bool> inBasis;
		mvtyp determinantOfBasis;
		bool ignoreLastRow;//indicates which mode we are in. The variable is set to true on construction and only changed when the lineality space is found.
	public:
		void computeRowBounds()
		{//CHECKCOLNEEDED
			assert(negativeRowBounds.size()==combinedMatrix.getHeight());
			for(int i=0;i<negativeRowBounds.size();i++)
				if(!ignoreRows[i])
					negativeRowBounds[i]=mvtyp::computeNegativeBound(&combinedMatrix[i][0],combinedMatrix.getWidth());
		}
		void computeLastRowBound() //as above but only for the last row
		{//CHECKCOLNEEDED
			assert(negativeRowBounds.size()==combinedMatrix.getHeight());
			int i=negativeRowBounds.size()-1;
			if(!ignoreRows[i])
				negativeRowBounds[i]=mvtyp::computeNegativeBound(&combinedMatrix[i][0],combinedMatrix.getWidth());
		}
//		Tableau& operator=(Tableau &&a)=default;
		Tableau& operator=(Tableau &&a)
		{
//			std::cerr<<"Tableau move assignment\n";
			combinedMatrix=std::move(a.combinedMatrix);
			negativeRowBounds=std::move(a.negativeRowBounds);
			isColumnInvalid=std::move(a.isColumnInvalid);
			ignoreRows=std::move(a.ignoreRows);
			basisIndices=std::move(a.basisIndices);
			inBasis=std::move(a.inBasis);
			determinantOfBasis=a.determinantOfBasis;
			ignoreLastRow=a.ignoreLastRow;
			return *this;
		}
		Tableau& operator=(Tableau const &a)=default;
		Tableau(MR *mr):
			combinedMatrix(mr),
			negativeRowBounds(0,mr),
			isColumnInvalid(0,0,mr),
			ignoreRows(0,0,mr),
			basisIndices(0,0,mr),
			inBasis(0,0,mr)
		{
//			std::cerr<<"Tableau constructor\n";
		}

		// Template for conversion
		template <typename otherTyp>
		explicit Tableau(Tableau<otherTyp> const &a, MR *mr=get_default_resource()):
			combinedMatrix(Matrix<mvtyp>(a.combinedMatrix,mr)),
			negativeRowBounds(Vector<mvtyp>(a.negativeRowBounds,mr)),
			isColumnInvalid(a.isColumnInvalid,mr),
			ignoreRows(a.ignoreRows,mr),
			basisIndices(a.basisIndices,mr),
			inBasis(a.inBasis,mr),
			determinantOfBasis(a.determinantOfBasis),
			ignoreLastRow(a.ignoreLastRow)
			{
			}

		void swapRows(int a, int b)//may screw up sign of determinant if called just once
		{
			combinedMatrix.swapRows(a,b);
			swap(negativeRowBounds[a],negativeRowBounds[b]);
			swap(ignoreRows[a],ignoreRows[b]);
			swap(basisIndices[a],basisIndices[b]);
		}
		Tableau(Matrix<mvtyp> const &M, bool appendIdentity, bool appendAllOnes=false, MR *mr=get_default_resource()):
			basisIndices(M.getHeight()+appendAllOnes,0,mr),
			ignoreLastRow(true),
			combinedMatrix(mr),
			negativeRowBounds(M.getHeight()+appendAllOnes,mr),
			isColumnInvalid(M.getWidth()+M.getHeight()+appendAllOnes,0,mr),
			ignoreRows(M.getHeight()+appendAllOnes,0,mr),
			inBasis(M.getWidth()+M.getHeight()+appendAllOnes,0,mr)
		{
			if(appendAllOnes)
			{
				combinedMatrix=Matrix<mvtyp>(M.getHeight()+1,M.getHeight()+1+M.getWidth(),mr);
				combinedMatrix.setSubMatrix(0,M.getHeight()+1,M.getHeight(),getWidth(),M);
				for(int i=0;i<M.getHeight()+1;i++)combinedMatrix[i][i]=1;
				for(int i=M.getHeight()+1;i<getWidth();i++)combinedMatrix[M.getHeight()][i]=1;
//				Matrix<mvtyp> M2=M;
//				M2.appendRow(Vector<mvtyp>::allOnes(M2.getWidth()));
//				combinedMatrix=combineLeftRight(M2.identity(M.getHeight()+1),M2,mr);
			}
			else
			{
				assert(0);// Not to say that this code is bad... but it is never called...
				combinedMatrix=combineLeftRight(M.identity(M.getHeight()),M,mr);											//MEM
			}
			assert(inBasis.size()==getWidth());
			for(int i=0;i<M.getHeight()+appendAllOnes;i++){basisIndices[i]=i;inBasis[i]=true;}
			determinantOfBasis=1;
			computeRowBounds();
		}
		Tableau(Tableau<mvtyp> const &a, MR *mr=get_default_resource()):
			Tableau(mr)
		{
			operator=(a);
		}
		Tableau(Tableau<mvtyp> &&a, MR *mr)NOEXCEPT_:
			Tableau(mr)
		{
			operator=(std::move(a));
		}
		Tableau(Tableau<mvtyp> &&other)NOEXCEPT_=default;  //CHECK IF THIS WORKS
	public:
		void regret(int i, int j)//Should call regret of superclass and do the last lines of adjustment here  //MAY OVERFLOW
		{
	//		std::cerr<<"Regret("<<i<<","<<j<<") called on:\n"<<combinedMatrix.toString()<<"\n";
	//		std::cerr<<"Regret In"<<matrixToString(combinedMatrix)<<"\n";
			assert(i!=j);
			int d=getHeight();
			if(!inBasis[i+d])
			{
	//			std::cerr<<"Case A\n";
				assert(CHECKCOL(i+d));assert(CHECKCOL(j+d));
				for(int k=0;k<getHeight()-1;k++)  // <------------------------------------
					combinedMatrix[k][i+d]=combinedMatrix[k][i+d]-combinedMatrix[k][j+d];
//				std::cerr<<"Result:\n"<<combinedMatrix.toString()<<"\n";
			}
			else
			{
				int I=-1;
				for(int k=0;k<basisIndices.size();k++)if(i+d==basisIndices[k]){I=k;}
				assert(I!=-1);
				if(signAt(I,j+d))
				{
//					std::cerr<<"Case B\n";
					exchange(I,j+d);
					assert(CHECKCOL(i+d));assert(CHECKCOL(j+d));
					for(int k=0;k<getHeight()-1;k++)  // <------------------------------------
						combinedMatrix[k][i+d]=combinedMatrix[k][i+d]-combinedMatrix[k][j+d];
				}
				else
				{
//					std::cerr<<"Case C\n";
					//Nu boer determinant gaa op i alle indgange i jte soejle
					for(int k=0;k<d-1;k++)// <---------------------------------------------
						if(k!=I)
						{
							assert(CHECKCOL(i+d));assert(CHECKCOL(j+d));
							combinedMatrix[k][i+d]=-combinedMatrix[k][j+d];
							assert(combinedMatrix[k][i+d].toInt64()%combinedMatrix[I][i+d].toInt64()==0);
							for(int a=0;a<getWidth();a++)
//								dotDivAssign(combinedMatrix[k][a],)
							{
								assert(CHECKCOL(a));assert(CHECKCOL(i+d));
								combinedMatrix[k][a]=combinedMatrix[k][a]-(combinedMatrix[k][i+d]/combinedMatrix[I][i+d])*combinedMatrix[I][a];
							}
//							combinedMatrix.madd(I,CircuitTableInt32(-combinedMatrix[k][i+d+1].v/combinedMatrix[I][i+d+1].v),k);
						}
				}
			}
			computeRowBounds();
		}
		std::string toString(Matrix<mvtyp> const *extra=0)
		{
//			std::cerr<</*combinedMatrix.toString()+*/toString2(basisIndices)+"\n";
			Matrix<mvtyp> firstRow(1,combinedMatrix.getWidth());
			for(int i=0;i<basisIndices.size();i++){firstRow[0][basisIndices[i]]=mvtyp(i);}
			if(extra)
			{
				firstRow=combineOnTop(*extra,firstRow);
			}

			Matrix<mvtyp> colors(combinedMatrix.getHeight(),combinedMatrix.getWidth());
			for(int i=0;i<colors.getHeight();i++)
				for(int j=0;j<colors.getWidth();j++)
				{
					int color=16*7;
					if(ignoreRows[i])
						color-=16;
					if(inBasis[j])
						color-=16*2;
//					if(inLineality[j])
//						color-=16*4;
					if(basisIndices[i]==j)
						color-=16*4-2;

					colors[i][j]=color;
				}


			return "Tableau:\n"+
			matrixToString/*<Matrix<mvtyp> >*//*(combineOnTop(firstRow,combinedMatrix))*/(combinedMatrix,&colors)+
			"Basis:\n"+toString2(basisIndices)+"\nDeterminant:\n"+determinantOfBasis.toString()+"\n";
		}
		void exchange(int i, int j) // make assigment basis[i]=j, updating matrix  //MAY OVERFLOW
		{
			assert(CHECKCOL(j));
//			std::cerr<<i<<" "<<j<<"BEFORE:"<<toString();
			mvtyp detMultiplier=combinedMatrix[i][j];
			typename mvtyp::Divisor divisorObject(determinantOfBasis);
//			std::cerr<<"Det:"<<determinantOfBasis.toString()<<"\n";
//			std::cerr<<divisorObject.toString()<<"\n";
			for(int k=0;k<getHeight()-ignoreLastRow;k++)
				if(!ignoreRows[k])
				if(k!=i)
				{
					assert(CHECKCOL(j));
					mvtyp temp=-combinedMatrix[k][j];
					mvtyp boundK=negativeRowBounds[k];
					mvtyp boundI=negativeRowBounds[i];
//					std::cerr<<"k"<<k<<" "<<temp.toString()<<" "<<boundK.toString()<<" "<<boundI.toString()<<"\n";
					//CHECKCOL NEEDED
//					computeRowBounds();//!!!!!!!!!!!!!!
					auto newBound=(temp.isZero())?
							mvtyp::scaleVector(&combinedMatrix[k][0],combinedMatrix[i][j],divisorObject,getWidth(),boundK)
							:
							mvtyp::dotDivVector(&combinedMatrix[k][0],&combinedMatrix[i][0],combinedMatrix[i][j],temp,divisorObject,getWidth(),boundK,boundI);
					negativeRowBounds[k]=newBound;
				}
			inBasis[basisIndices[i]]=false;
			basisIndices[i]=j;
			inBasis[j]=true;
			determinantOfBasis=detMultiplier;
//			std::cerr<<"AFTER:"<<toString();
//			std::cerr<<"Det:"<<determinantOfBasis.toString()<<"\n";
		}
	};

	template<class mvtyp> class GeneratedCone:public Tableau<mvtyp>{
	public:
		using Tableau<mvtyp>::getWidth;
		using Tableau<mvtyp>::getHeight;
		using Tableau<mvtyp>::inBasis;
		using Tableau<mvtyp>::ignoreLastRow;
		using Tableau<mvtyp>::basisIndices;
		using Tableau<mvtyp>::determinantOfBasis;//!
		using Tableau<mvtyp>::signAt;
		using Tableau<mvtyp>::combinedMatrix;//!
		using Tableau<mvtyp>::ignoreRows;
		using Tableau<mvtyp>::assignEntry;//!
		using Tableau<mvtyp>::swapRows;//
		//using Tableau<mvtyp>::exchange;//
		using Tableau<mvtyp>::computeRowBounds;
		using Tableau<mvtyp>::computeLastRowBound;
		using Tableau<mvtyp>::CHECKCOL;
		using Tableau<mvtyp>::makeColumnInvalid;
	private:
		int state;
		/*
		 * The states are:
		 * 0: initial state
		 * 1: the span of the cone is known (a basis can be requested, and so can a basis for the orthogonal complement)
		 * 2: the lineality space is known and also a generic supporting hyperplane
		 * 3: the extreme rays are known
		 * ..: the set of facets is known  - not clear if such information should be stored.
		 *     It is also not clear if they should be forced to be in the span of the cone,
		 *     they should be made primitive or their (unique) sum should be stored.
		 */
		int ambientDimension;//valid if state>=0
		int dimension;//valid if state>=1
		// Markings for rows:
		pmrvector<bool> inOrthogonalComplement;//valid if state>=1
		// Markings for columns:
		pmrvector<bool> nonExtreme;//knownToBeNonExtreme
		pmrvector<bool> inLineality;
		Vector<mvtyp> temporaryNormal;
		/*const*/ Matrix<mvtyp> originalMatrix; // we almost don't have to store this matrix. When for example rays for the cone are requested, it however is easiest to return a submatrix of this stored matrix rather than changing basis.
		static uint64_t bitsSet(pmrvector<bool> const &v)
		{
			uint64_t r=0;
			for(int i=0;i<v.size();i++)if(v[i])r++;
			return r;
		}
	public:
		int getState()const  // This function should only be used for choosing strategies.
		{
			return state;
		}
		void exchange(int i, int j)
		{
			if(0)
			{
				static FrequencyTable widthHist("Width");
				widthHist.record(combinedMatrix.getWidth());
				static FrequencyTable heightHist("Height");
				heightHist.record(combinedMatrix.getHeight());
				int nInvalidColumns=0;
				for(int i=0;i<combinedMatrix.getWidth();i++)
					if(Tableau<mvtyp>::isColumnInvalid[i])nInvalidColumns++;
				static FrequencyTable nInvalidHist("Number of invalid columns");
				nInvalidHist.record(nInvalidColumns);
				static FrequencyTable nInLinealitySpace("Number of linealityspace columns");
				nInLinealitySpace.record(bitsSet(inLineality));
				static FrequencyTable nonExtremeHist("Known to be non-extreme");
				nonExtremeHist.record(bitsSet(nonExtreme));
				static FrequencyTable ignoreRowsHist("ignoreRows");
				ignoreRowsHist.record(bitsSet(ignoreRows));
				static FrequencyTable inOrthogonalComplementHist("inOrthogonalComplement");
				inOrthogonalComplementHist.record(bitsSet(inOrthogonalComplement));

			}

			Tableau<mvtyp>::exchange(i,j);
		}
		/**
		 * This method checks if it happens that a column vector is trivially redundant because it is non-negative
		 * in the current basis. If so, it is marked as nonExtreme.
		 */
		void __attribute__ ((noinline)) redundancyCheckAtCurrentBasis()
		{
			return;
			int lowerIndex=Tableau<mvtyp>::getHeight();
			int upperIndex=Tableau<mvtyp>::getWidth();

			if(state>=1 /*&& ignoreLastRow*/)
			{
				for(int j=lowerIndex;j<upperIndex;j++)
					if(!nonExtreme[j] && !inBasis[j] && !inLineality[j])
					{
						bool isNonNegative=true;
						for(int i=0;i<getHeight()-ignoreLastRow;i++)
								if(!ignoreRows[i])
								if(signAt(i,j)==-determinantOfBasis.sign())
								{
									isNonNegative=false;
									break;
								}
						if(isNonNegative)
						{
//							std::cerr<<"Simplifying column "<<j<<" (column is non-extreme):\n";
							nonExtreme[j]=true;
							makeColumnInvalid(j);
						}
					}
			}
		}
		/* Equivalent of subtracting column j from column i in original data (under certain conditions)
		 * Notice that the indices correspond to column indices of combineLeftRight(M,L).
		 * That means that since the identity was preprended, the indices must be adjusted.
		 *
		 * For this method to be called several conditions must be satisfied:
		 * * The dimension of the cone is preserved
		 * * column i and j are extreme rays (could they be positive scalar multiple of other columns?)
		 * * The lineality space is preserved
		 * Consequences:
		 * * new vector is not in lineality space
		 * * new vector is still extreme
		 */
		void regret(int i, int j)
		{
			Tableau<mvtyp>::regret(i,j);
			for(int k=0;k<originalMatrix.getHeight();k++)
				originalMatrix[k][i]-=originalMatrix[k][j];
			//HERE THE STATE SHOULD DROP TO 2, SINCE SOME EXTREME RAYS COULD HAVE BECOME NONEXTREME
			if(state>2)state=2;
		}
		std::string toString()
		{
			std::stringstream s;
			s<<"state:"<<state<<"\n";
//			s<<"numberOfOriginalLineGenerators"<<numberOfOriginalLineGenerators<<"\n";
			s<<"nonExtreme:";for(int i=0;i<getWidth();i++)s<<int(nonExtreme[i]);
/*			{
				Matrix<mvtyp> temp(0,getHeight());
				for(int i=0;i<getWidth();i++)if(!nonExtreme[i])temp.appendRow(combinedMatrix.column(i));
				s<<matrixToString<Matrix<mvtyp> >(temp.transposed());
			}*/
			s<<"inLinealin:";for(int i=0;i<getWidth();i++)s<<inLineality[i];s<<"\n";
			s<<"inBasis:   ";for(int i=0;i<getWidth();i++)s<<int(inBasis[i]);s<<"\n";
			s<<"ignoreRows:      ";for(int i=0;i<getHeight();i++)s<<int(ignoreRows[i]);s<<"\n";
			s<<"inOrthComplement:";for(int i=0;i<getHeight();i++)s<<int(inOrthogonalComplement[i]);s<<"\n";

			Matrix<mvtyp> extra(4,getWidth());
			for(int i=0;i<getWidth();i++)
			{
				extra[0][i]=i;
				extra[1][i]=mvtyp(nonExtreme[i]);
				extra[2][i]=mvtyp(inLineality[i]);
				extra[3][i]=mvtyp(inBasis[i]);
			}



			return Tableau<mvtyp>::toString(&extra)+s.str()+"originalMatrix:"+matrixToString(originalMatrix);
//			return "Tableau:\n"+combinedMatrix.toString()+"Basis:\n"+vectorToString(basisIndices)+"\nDeterminant:\n"+determinantOfBasis.toString()+"\n";
		}
		void save(std::ostream &s)const
		{
			s<<"("<<ambientDimension<<"\n";
			s<<","<<getRayGenerators().transposed();
			s<<","<<getLineGenerators().transposed()<<")\n";
		}
		static GeneratedCone load(std::istream &s)
		{
			auto a=parseSequence(s);
			assert(a.size()==3);
			stringstream ss(a[0]);
			int ambientDimension;ss>>ambientDimension;
			stringstream ss1(a[1]);
			stringstream ss2(a[2]);
			auto rays=Matrix<mvtyp>::readMatrix(ss1,ambientDimension).transposed();
			auto lines=Matrix<mvtyp>::readMatrix(ss2,ambientDimension).transposed();
			return GeneratedCone(rays,lines);
		}
		bool isColumnZero(int i)//index to combinedMatrix, ignoreRows are ignored  REPLACE BY CODE BELOW
		{
			for(int a=0;a<ambientDimension;a++)
				if(!ignoreRows[a])
					if(signAt(a,i))return false;
			return true;
		}
		int rowIndexOfNonZeroEntryInColumn(int i)//index to combinedMatrix, ignoreRows are ignored
		{
			for(int a=0;a<ambientDimension;a++)
				if(!ignoreRows[a])
					if(signAt(a,i))return a;
			return -1;
		}
/*
 * The constructor only accepts a matrix with columns being generators of the cone.
 * A future version should allow also to specify lines that are known to be in the cone.
 */
		explicit GeneratedCone(Matrix<mvtyp> const &M, Matrix<mvtyp> const &L, MR *mr=get_default_resource(), MR *mr2=get_default_resource()):
			Tableau<mvtyp>(combineLeftRight(M,L,mr2), true,true, mr),//MEM
			originalMatrix(combineLeftRight(M,L,mr),mr),//MEM
			nonExtreme(M.getHeight()+M.getWidth()+L.getWidth()+1,0,mr),
			inLineality(M.getHeight()+M.getWidth()+L.getWidth()+1,0,mr),
			temporaryNormal(M.getHeight(),mr),
			state(0),
			inOrthogonalComplement(M.getHeight()+1,0,mr),
			ambientDimension(M.getHeight())//,
		{
			// We make sure that known part of the lineality space is part of basis
			for(int i=getWidth()-L.getWidth();i<getWidth();i++)
				{
				inLineality[i]=true;
				int a=rowIndexOfNonZeroEntryInColumn(i);
				if(a>=0)
				{
					exchange(a,i);
					ignoreRows[a]=true;
				}
				makeColumnInvalid(i);
				}
//			std::cerr<<"CONSTRUCTING from\n"<<M.transposed().toString();
//			std::cerr<<toString()<<"\n";
		}
		GeneratedCone(MR *mr)NOEXCEPT_://is this constructor needed - and how should the initialisation be?
			Tableau<mvtyp>(mr),
			originalMatrix(mr),
			nonExtreme(mr),
			inLineality(mr),
			temporaryNormal(0,mr),
			state(0),
			inOrthogonalComplement(mr),
			ambientDimension(0)//,
//			numberOfOriginalLineGenerators(0)
		{
//			std::cerr<<"GeneratedCone constructor\n";
		}
		GeneratedCone()NOEXCEPT_//is this constructor needed - and how should the initialisation be?
		{
		}
		explicit GeneratedCone(Matrix<mvtyp> const &M, MR *mr=get_default_resource()):
 			Tableau<mvtyp>(M, true,true,mr),
			originalMatrix(M,mr),
			nonExtreme(M.getHeight()+M.getWidth()+1,0,mr),
			inLineality(M.getHeight()+M.getWidth()+1,0,mr),
			temporaryNormal(M.getHeight(),mr),
			state(0),
			inOrthogonalComplement(M.getHeight()+1,0,mr),
			ambientDimension(M.getHeight())//,
//			numberOfOriginalLineGenerators(0)
		{
		}
//		GeneratedCone &operator=(GeneratedCone&& other)NOEXCEPT_=default; //not generated by default???

		GeneratedCone &operator=(GeneratedCone&& other)NOEXCEPT_
		{
//			std::cerr<<"GeneratedCone move assignment\n";
			Tableau<mvtyp>::operator=(std::move(other));
			originalMatrix=std::move(other.originalMatrix);
			nonExtreme=std::move(other.nonExtreme);
			inLineality=std::move(other.inLineality);
			temporaryNormal=std::move(other.temporaryNormal);
			state=other.state;
			inOrthogonalComplement=std::move(other.inOrthogonalComplement);
			ambientDimension=other.ambientDimension;
			dimension=other.dimension;
			return *this;
		}
		GeneratedCone &operator=(const GeneratedCone &other)NOEXCEPT_
		{
//			std::cerr<<"move init4?"<<&other<<"  "<<this<<"\n";
			if(&other==this)return *this;
//			GeneratedCone temp(other);
//			std::swap(temp,*this);
//			return *this;

			Tableau<mvtyp>::operator=(other);
			originalMatrix=other.originalMatrix;
			nonExtreme=other.nonExtreme;
			inLineality=other.inLineality;
			temporaryNormal=other.temporaryNormal;
			state=other.state;
			inOrthogonalComplement=other.inOrthogonalComplement;
			ambientDimension=other.ambientDimension;
			dimension=other.dimension;
			return *this;
		}
		explicit GeneratedCone(GeneratedCone const &a, MR *mr=get_default_resource())NOEXCEPT_:
			GeneratedCone(mr)
		{
			operator=(a);
		}
		GeneratedCone(GeneratedCone &&other)NOEXCEPT_=default;  //CHECK IF THIS WORKS
/*		GeneratedCone(GeneratedCone &&other)NOEXCEPT_
			:GeneratedCone(other.get_memory_resource())
			 {
//				std::cerr<<"move init5?"<<&other<<"  "<<this<<"\n";
				operator=(std::move(other));
			 }*/
		GeneratedCone(GeneratedCone &&other, MR *mr)NOEXCEPT_
			:GeneratedCone(mr)
			 {
//				std::cerr<<"move init5?"<<&other<<"  "<<this<<"\n";
				operator=(std::move(other));
			 }
		~GeneratedCone()NOEXCEPT_{
/*				std::cerr<<"GeneratedCone destructor:\n";
				if(temporaryNormal.size()>0)
					std::cerr<<(int*)&(temporaryNormal[0])<<"\n";*/

		}
		// Template for conversion
		template<typename> friend class GeneratedCone;
		template <typename otherTyp>
		explicit GeneratedCone(GeneratedCone<otherTyp> const &a, MR *mr=get_default_resource()):
		Tableau<mvtyp>(a,mr),
		state(a.state),
		ambientDimension(a.ambientDimension),
		dimension(a.dimension),
		inOrthogonalComplement(a.inOrthogonalComplement,mr),
		nonExtreme(a.nonExtreme,mr),
		inLineality(a.inLineality,mr),
		temporaryNormal(a.temporaryNormal,mr),
		originalMatrix(a.originalMatrix,mr)
			{
			}
	private:
		bool loose(int basisIndex, int lowerIndex, int upperIndex, pmrvector<bool> const &disallowed)
		{
//			std::cerr<<"loose\n";
			int basisIndex2=-1;
			for(int i=0;i<basisIndices.size();i++)if(basisIndex==basisIndices[i]){basisIndex2=i;}
			if(ignoreLastRow){assert(basisIndex2!=getHeight()-1);}

//			std::cerr<<"loose called"<<basisIndex<<"lower"<<lowerIndex<<"upper"<<upperIndex<<"s\n";
//			std::cerr<<toString();
			for(int i=lowerIndex;i<upperIndex;i++)
				if(!disallowed[i]&&!inLineality[i])
				if((!inBasis[i])&&signAt(basisIndex2,i))
				{
//					std::cerr<<"EXCHANGE!\n";

					exchange(basisIndex2,i);
					return true;
				}
//			std::cerr<<"returning false\n";
			return false;
		}
		bool doesIthsLiftsHigherThanSimplicialCone(int i)const
		{
			assert(!inBasis[i]);
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			for(int j=upperIndex-1;j>=lowerIndex;j--)
			{
				if(inBasis[j])
				{
					int J=0;
					for(;J<getHeight()-ignoreLastRow;J++)if(basisIndices[J]==j)break;
					assert(J!=getHeight());
					if(!ignoreRows[J])
					{
						assert(CHECKCOL(i));
						int t=(determinantOfBasis.sign()*combinedMatrix[J][i].sign());
						if(t<0)return true;
						if(t>0)return false;
					}
				}
				if(j==i) //nothing to check since ith is lifted positively
				{
					return true;//break; //If this is reached we are ok
				}
			}
			assert(0);
		}
		void checkInLexTriangulation()
		{
			for(int i=combinedMatrix.getWidth();i<combinedMatrix.getHeight();i++)
				if(!inBasis[i])
					if(!doesIthsLiftsHigherThanSimplicialCone(i))
					{
						std::cerr<<i<<"th does not lift higher\n";
						std::cerr<<toString();
						std::cerr<<"\n";
						assert(0);
					}
//			std::cerr<<"ok\n";
		}

		bool compareLex(int i, int j, int violatedIndex, int detSign, int ignoredTargetIndex) // Some signs may be wrong in this code:
		{
			assert(signAt(violatedIndex,i)==signAt(violatedIndex,j));
			assert(!inBasis[i]);
			assert(!inBasis[j]);
			assert(signAt(violatedIndex,ignoredTargetIndex));

#if 1
			int min=(i<j)?i:j;
			for(int k=0;k<min;k++)
				if(inBasis[k])
				{
					int I=0;
					while(basisIndices[I]!=k)I++;
					if(!ignoreRows[I] && (!ignoreLastRow||I!=combinedMatrix.getHeight()-1))
					{
						assert(CHECKCOL(i));assert(CHECKCOL(j));
						auto ival=combinedMatrix[violatedIndex][j]*combinedMatrix[I][i];
						auto jval=combinedMatrix[violatedIndex][i]*combinedMatrix[I][j];
						if(ival<jval)return signAt(violatedIndex,ignoredTargetIndex)==1;
						if(jval<ival)return signAt(violatedIndex,ignoredTargetIndex)==-1;
					}
				}
			return (i<j)^(signAt(violatedIndex,min)!=detSign);//^(signAt(violatedIndex,ignoredTargetIndex)==-1);
#else
			int max=(i<j)?j:i;
			for(int k=inBasis.size()-1;k>max;k--)
				if(inBasis[k])
				{
					int I=0;
					while(basisIndices[I]!=k)I++;
					if(!ignoreRows[I] && (!ignoreLastRow||I!=combinedMatrix.getHeight()-1))
					{
						auto ival=combinedMatrix[violatedIndex][j]*combinedMatrix[I][i];
						auto jval=combinedMatrix[violatedIndex][i]*combinedMatrix[I][j];
						if(ival<jval)return signAt(violatedIndex,ignoredTargetIndex)==-1;
						if(jval<ival)return signAt(violatedIndex,ignoredTargetIndex)==1;
					}
				}
			return (i<j)^(signAt(violatedIndex,max)==detSign)^(signAt(violatedIndex,ignoredTargetIndex)==-1);
#endif
		}
		/**
		 * Single pivot step. The replacing rule here should come from lexicograph perturbation and should guarantee anti-cycling rule.
		 * Will replace violatedIndex in the basis with something else.
		 */
		bool replaceLex(int violatedIndex, int lowerIndex, int upperIndex, pmrvector<bool> &ignoredColumns, int ignoredTargetIndex)
		{
//			std::cerr<<"call\n";
			int candidate=-1;
			auto detSign=determinantOfBasis.sign();

			for(int i=lowerIndex;i<upperIndex;i++)
				if(!ignoredColumns[i] && !inLineality[i] && i!=ignoredTargetIndex && !inBasis[i])
					if(detSign==-signAt(violatedIndex,i))
					{
						if(candidate==-1)
							candidate=i;
						else
						{
							if(compareLex(candidate,i,violatedIndex,detSign,ignoredTargetIndex))
								candidate=i;
						}
					}
	//		std::cerr<<"ret\ncandi:"<<candidate<<"\n"<<toString();
			if(candidate==-1)return false;
	//		std::cerr<<"Exchanging ViolatedIndex"<<violatedIndex<<"basisIndices[violatedIndex]"<<basisIndices[violatedIndex]<<"candidate"<<candidate<<"\n";
			exchange(violatedIndex,candidate);
			redundancyCheckAtCurrentBasis();//possibility of marking column as non-extreme
			return true;
		}
		/**
		 * Single pivot step. The replacing rule here seems to be Bland's rule and should guarantee anti-cycling rule.
		 * Will replace violatedIndex in the basis with something else.
		 */
		bool replaceBland(int violatedIndex, int lowerIndex, int upperIndex, pmrvector<bool> &ignoredColumns, int ignoredTargetIndex)
		{
			auto detSign=determinantOfBasis.sign();
			for(int i=lowerIndex;i<upperIndex;i++)
				if(!ignoredColumns[i] && !inLineality[i] && !nonExtreme[i] && i!=ignoredTargetIndex && !inBasis[i])
					if(detSign==-signAt(violatedIndex,i))
					{
						exchange(violatedIndex,i);
						redundancyCheckAtCurrentBasis();//possibility of marking column as non-extreme
						return true;
					}
			return false;
		}
		pmrvector<int> shadowVertexPerturbation;
		void setShadowVertexStart()
		{
			shadowVertexPerturbation=basisIndices;//alloc
		}
		/*
		 * Combined with replaceBland() above, this function is an implementation of Schrijver version of Bland's rule (Theorem 7.1 in Schrijver's book).
		 * It attempts to change the current simplicial cone into a simplicial cone that will contain the column vector
		 * indexed by targetIndex. This may fail, in which case a non-negative number is returned indexing a separating
		 * hyperplane. On success -1 is returned.
		 *
		 * Notice that the algorithm works even in the case where the configuration is not pointed.
		 */
		int towards(int targetIndex, int lowerIndex, int upperIndex, pmrvector<bool> &ignoredColumns)
		{
			int numberOfIterations=0;
			bool shadowVertex=false;true;
			bool bland=true;false;
			assert(!bland||!shadowVertex);
//			std::cerr<<"BLAND:"<<bland<<"\n";
			if(shadowVertex)setShadowVertexStart();
//			computeRowBounds();//!!!!!!!!
			while(1)
			{
//				checkInLexTriangulation();

				if(0)
				{
					if(!(numberOfIterations&511))std::cerr<<"Number of iterations:"<<numberOfIterations<<"\n";
				numberOfIterations++;
if(0)				if(numberOfIterations>=500000)
				{
//					printf("%c[%d", 0x1B, 32);
					std::cerr<<setColor(1)<<"Looping forever on:"<<this->toString()<<"\n";
					if(numberOfIterations>=500010)
						assert(0);
				}
				}
				//check whether targetIndex column is in simplex of basis
				int violatedIndex2=1000000000;
				int violatedIndex=-1;

				if(shadowVertex)
				{
					for(int i=0;i<getHeight()-ignoreLastRow;i++)
						if(!ignoreRows[i])
							if(determinantOfBasis.sign()*signAt(i,targetIndex)==-1)
							{
								if(violatedIndex==-1)
									violatedIndex=i;
								else
								{
									bool isBetter=false;
									for(int perturbationIndex=0;perturbationIndex<shadowVertexPerturbation.size();perturbationIndex++)
									{
										int d=combinedMatrix.getHeight();
										CHECKCOL(d-1);
										CHECKCOL(shadowVertexPerturbation[perturbationIndex]);
										int s=determinantOfBasis.sign()*determinantSign3(combinedMatrix[violatedIndex][d-1],combinedMatrix[i][shadowVertexPerturbation[perturbationIndex]],combinedMatrix[i][d-1],combinedMatrix[violatedIndex][shadowVertexPerturbation[perturbationIndex]]);
										if(s<0){isBetter=true;break;}
										if(s>0){isBetter=false;break;}
									}
									if(isBetter)violatedIndex=i;
								}
							}
				}
				else
					for(int i=0;i<getHeight()-ignoreLastRow;i++)
						if(!ignoreRows[i])
						{
							if(determinantOfBasis.sign()*signAt(i,targetIndex)==-1)
								if(basisIndices[i]<violatedIndex2)
									{
										violatedIndex=i;
										violatedIndex2=basisIndices[i];
									}
						}
//				std::cerr<<"vio:"<<violatedIndex<<"\n";
				if(violatedIndex==-1)
				{
//					std::cerr<<"Inside cone!\n";
					return -1;
				}
				else
				{

					if(bland?!replaceBland(violatedIndex,lowerIndex,upperIndex,ignoredColumns,targetIndex):
							 !replaceLex  (violatedIndex,lowerIndex,upperIndex,ignoredColumns,targetIndex))
					{
//						std::cerr<<"Separating hyperplane found:\n";
//						std::cerr<<combinedMatrix.submatrix(violatedIndex,0,violatedIndex+1,getHeight()).toString()<<"\n";
						return violatedIndex;
					}
					numberOfPivotingSteps++;
				}
			}
			return -2;
		}

		/**
		 * Below is the private code that takes the cone data between states. The reason for
		 * making these private is that the user should not deal with the states. Rather
		 * when requesting a property of the cone, it will automatically change state so that
		 * the property becomes known.
		 */
		void findOrthogonalComplementAndDimension()
		{
//			std::cerr<<"TEST\n";
//			checkInLexTriangulation();
//			checkInLexTriangulation();
//			std::cerr<<"findOrthogonalComplementAndDimension()\n";
			assert(ignoreLastRow);
			assert(state>=0);
			/*
			 * If the cone is lower dimensional, to make it full dimensional, some rows must be marked as ignored.
			 * When this has been done, generators for the orthogonal complement can be read off.
			 */
			int dim=getHeight()-1;
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			int d=getHeight();
			for(int i=0;i<d-1;i++)
				if(!ignoreRows[i])//could already be set if part of known lineality space
			{
//					checkInLexTriangulation();
				int j;
				for(j=lowerIndex;j<upperIndex;j++)
					if(signAt(i,j))
						break;
				if(j==upperIndex)
				{
//					std::cerr<<"FLAT MUST IGNORE ROW"<<i<<"\n";
//					std::cerr<<"CAUSE:"<<combinedMatrix.submatrix(i,0,i+1,d).toString()<<"\n";
					ignoreRows[i]=true;
					inOrthogonalComplement[i]=true;
					dim--;
				}
				else
				{
//					checkInLexTriangulation();
//					std::cerr<<"Ex"<<i<<j<<"\n";
					this->exchange(i,j);
					redundancyCheckAtCurrentBasis();
//					checkInLexTriangulation();
				}
			}
			dimension=dim;
//			checkInLexTriangulation();
//			redundancyCheckAtCurrentBasis();
		}
		void subtractOneFromIthEntryOfLastRowOfOriginalCombinedMatrix(int i)  //MAY OVERFLOW
		{
			int d=getHeight();
			int upperIndex=getWidth();
			assert(CHECKCOL(i));assert(CHECKCOL(d-1));
			for(int a=0;a<d;a++)combinedMatrix[a][i]-=combinedMatrix[a][d-1];
			if(inBasis[i])
			{
				int index=0;
				for(index=0;index<d;index++)if(basisIndices[index]==i)break;
				assert(index!=d);

				for(int j=0;j<upperIndex;j++)
				{
					assert(CHECKCOL(j));
					combinedMatrix[d-1][j]+=combinedMatrix[index][j];
				}

				{
					//check that the column is correct:
					for(int i=0;i<d;i++)
						if(i==index)
						{
							assert(CHECKCOL(basisIndices[index]));
							assert((combinedMatrix[index][basisIndices[index]]-determinantOfBasis).isZero());
						}
						else
						{
							assert(CHECKCOL(basisIndices[index]));
							assert((combinedMatrix[i][basisIndices[index]]).isZero());
						}
				}
			}
			computeRowBounds();
		}
		void correctLastRowForBasisElementI(int i)  //MAY OVERFLOW
		{
			int d=getHeight();

			assert(CHECKCOL(basisIndices[i]));
			if(!combinedMatrix[d-1][basisIndices[i]].isZero())
			{
				for(int j=0;j<getWidth();j++)
					if(!nonExtreme[j] && !inLineality[j])// We should not correct columns that we may ignore
				{
					assert(CHECKCOL(j));
					combinedMatrix[d-1][j]=combinedMatrix[d-1][j]-combinedMatrix[i][j];
				}
			}
		}
		/*
		 * The ith column is assumed to be in the lineality space. If this is nonzero outside current
		 * linealityspace, this function will force the ith column into the basis basis.
		 * Last row is ignored.
		 */
		void forceLinealityToBasis(int i)
		{
			int d=getHeight();
			int nonZeroIndex=0;
			for(nonZeroIndex=0;nonZeroIndex<d-1;nonZeroIndex++)
				if(!ignoreRows[nonZeroIndex])
					if(signAt(nonZeroIndex,i))break;
			if(nonZeroIndex!=d-1)
			{
				exchange(nonZeroIndex,i);
				ignoreRows[nonZeroIndex]=true;
			}
		}
		/*
		 * Last row is ignored.
		 */
		void markZeroColumnsAsInLineality(int lowerIndex, int upperIndex)
		{
			for(int i=lowerIndex;i<upperIndex;i++)
				if(!inLineality[i])
				{
					bool isZeroColumn=true;
					for(int j=0;j<getHeight()-1;j++)
						if(!ignoreRows[j])if(signAt(j,i))isZeroColumn=false;
					if(isZeroColumn)
					{
						inLineality[i]=true;
						makeColumnInvalid(i);
					}
				}
		}
		void findBasisForLinealitySpaceAndItsDimension()//  MAY OVERFLOW (IN LAST LOOP)
		{// Should be allowed to assign entry of matrix and to change sign of column
			assert(state>=1);
			ignoreLastRow=false;
			// during this method d-1 can leave the basis, but it is put back in at the end.
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			int d=getHeight();


			// use information already stored in inLineality
			for(int i=d;i<getWidth();i++)
				if(inLineality[i])
					forceLinealityToBasis(i);


			// We need to setup the last row
//			assert(basisIndices[d-1]==d-1);
//			for(int j=d-1;j<upperIndex;j++)combinedMatrix[d-1][j]=determinantOfBasis;

//			computeRowBounds();//!!!!!!!
			int violatedIndex=-2;
			do
			{
				violatedIndex=-2;
//				std::cerr<<"BL:"<<toString();
//				std::cerr<<"BLA\n";
				//we need to check if some is in span of added vector
				markZeroColumnsAsInLineality(lowerIndex,upperIndex);
				bool somethingNotKnownToBeInLineality=false;
//				for(int i=lowerIndex;i<upperIndex;i++)
//					if(!inLineality[i]){somethingNotKnownToBeInLineality=true;break;}
//				if(!somethingNotKnownToBeInLineality)break;
				for(int j=0;j<upperIndex;j++)assignEntry(d-1,j,(j==d-1 || (j>=d && !inLineality[j]))?determinantOfBasis:mvtyp(0));
				for(int i=0;i<d-1;i++)assignEntry(i,d-1,0);
//				computeRowBounds();// not needed
				for(int i=0;i<d-1;i++)if(!inLineality[basisIndices[i]])correctLastRowForBasisElementI(i);
				computeLastRowBound();
				//				std::cerr<<facetStatus();
				if(loose(d-1,lowerIndex,upperIndex,nonExtreme))
				{
	//				std::cerr<<"WE LOST ADDED\n";
	//				std::cerr<<":"<<toString();

					violatedIndex=towards(d-1,lowerIndex,upperIndex,nonExtreme);
//					std::cerr<<"FOUND SIMPLEX"<<foundSimplex<<"\n";
//					std::cerr<<"violated index"<<violatedIndex<<"\n";
//					std::cerr<<toString();

					if(violatedIndex!=-1)
						{
//						assert(0);
						break;
						}
					// Now we have found a positive circuit involving the basis elements with non-zero coefficient in special column and the special column
					// We must mark all that are in span
//					vector<bool> goesToLineality(upperIndex);
					for(int i=lowerIndex;i<upperIndex;i++)
					{
//						goesToLineality[i]=0;
						bool isInSpan=true;
						for(int j=0;j<d;j++)
							if(!ignoreRows[j])
								if(signAt(j,d-1)==0)
									if(signAt(j,i))isInSpan=false;
						if(isInSpan)
							{
//								if(inLineality[i]==false)goesToLineality[i]=1;
								inLineality[i]=true;
								makeColumnInvalid(i);
							}
					}
					// We must swap in the special column

//					std::cerr<<"MUSTSWAPINSPECIAL\n"<<toString();
					// find special rows
					int a;
					for(a=0;a<d;a++)if(!ignoreRows[a])if(signAt(a,d-1))break;
					int b;
					for(b=a+1;b<d;b++)if(!ignoreRows[b])if(signAt(b,d-1))break;

//					std::cerr<<"a"<<a<<"b"<<b<<"d"<<d<<"\n";

					if(a==d || b==d)
					{
						std::cerr<<toString2(basisIndices)<<"\n";
//						std::cerr<<toString();
					}
					assert(a!=d);
					assert(b!=d);

					if(b!=d-1)
					{
						//we now cycle a,b,d-1
						swapRows(a,b);
						swapRows(b,d-1);
						b=d-1;
					}
					// now a and b are non-zero rows of special column
					// we now know that b=d-1
					// we let a be the row to ignore and let b=d-1 be special


					int makePermanentIndex=a;

//					assert(newSpecialPos<d);
					//					assert(combinedMatrix[d-1][d-1].isNonZero());
//					int makePermanentIndex;
//					for(makePermanentIndex=0;makePermanentIndex<d;makePermanentIndex++)if(combinedMatrix[makePermanentIndex][d-1].isNonZero())break;
//					std::cerr<<toString();
//					std::cerr<<"EXCHANGE!\n";
					exchange(d-1,d-1);
//					exchange(newSpecialPos,d-1);
//					std::cerr<<toString();

					// We must make one of the other circuit members a permanent member of the basis by setting ignoreRow
//					assert(0);
					assert(makePermanentIndex!=d-1);
					assert(makePermanentIndex!=d);
					ignoreRows[makePermanentIndex]=true;
//					std::cerr<<"IGNORING\n";
					for(int i=0;i<getHeight()-1;i++)
						if(inLineality[basisIndices[i]])
							{
								ignoreRows[i]=true;// lineality space could be higher dimensional!
							}

					//WE CHANGE THE ORIGINAL MATRIX M by making zeros in the lower appended row at the places where the discovered
					//lineality space vectors were to make sure that in next iteration we will consider new vectors.
					//The update has two steps. Change M by subtracting the d-1st column from the lineality columns.
					// Because this may destroy the property that the selected columns are a scaled standard basis we thereafter
					// need to make sure that the matrix again becomes a scaled identity matrix in the basis columns

/*
					for(int i=lowerIndex;i<upperIndex;i++)
						if(goesToLineality[i])
							subtractOneFromIthEntryOfLastRowOfOriginalCombinedMatrix(i);
*/

//					std::cerr<<toString();
				}
				else
				{
					violatedIndex=d-1;
//					std::cerr<<"WE DID NOT LOOSE ADDED\n";
//					std::cerr<<toString()<<"\n";
//					std::cerr<<"WHAT DO WE DO NOW?\n"<<"It seems that at thic point the lineality space is zero-dimensional, and we should break out of the loop.\n";

//					std::cerr<<"We are in this case because the d-1st column is not even in the span of the last columns and therefore cannot be in the non-negative span of the columns.\n";
//					assert(0);

					break;
				}
//				std::cerr<<"L---------------------\n";
			}while(1);
//			std::cerr<<"E---------------------\n";

			// There are two ways to break out of the loop above. In both cases violatedIndex indexes a row containing a certificate
			// that the lineality space cannot be larger. We swap that row to the last row and keep the last row forever.

			if(violatedIndex!=d-1)swapRows(violatedIndex,d-1);

			assert(signAt(d-1,d-1));
			// We make sure that d-1 is in basis at position d-1.
			exchange(d-1,d-1);

			assert(CHECKCOL(d-1));
			if(combinedMatrix[d-1][d-1].isPositive())for(int i=0;i<d;i++)
				{
				assert(CHECKCOL(i));
				combinedMatrix[d-1][i]=-combinedMatrix[d-1][i];
				}
//			assert(isLastRowGenericNormal());

			ignoreLastRow=true;
			computeRowBounds();
		}
	public:
		/*
		 * Equivalent to letting ith column of the original matrix generate a line of the original matrix.
		 * Assumptions:
		 * * the column is not already in the lineality space
		 */
				void regretMakeLine(int i, bool tryToReachState2=false)//it turns out to be better to not increase the state.
				{
					int d=getHeight();
					inLineality[i+d]=true;
					makeColumnInvalid(i+d);
					auto oldState=state;
					state=0;
					if(oldState>0)
					{ //the orthogonal complement will not change and therefore this call is unnecessary
						state=1;
					}
					if(tryToReachState2)
						if(oldState>1)
						{
#if 0
							findBasisForLinealitySpaceAndItsDimension();//We know that the lineality is only extended by ith vector
#else
							forceLinealityToBasis(i+d);
							markZeroColumnsAsInLineality(getHeight(),getWidth());
#endif
							state=2;
						}
				}
				/*
				 * Some thoughts about the above.
				 * Ideally, only the ith column enters the lineality space
				 *          and nonExtreme columns remain nonExtreme
				 * It then seems that making ith column part of the basis and marking it as inLineality would suffice.
				 *
				 *
				 */
	private:
		void findExtremeRays()
		{
//			computeRowBounds();//!!!!!!!!
			assert(state>=2);

			redundancyCheckAtCurrentBasis();///???

			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			int d=getHeight();

			int neededExtremeRays=getDimension()-getDimensionOfLinealitySpace();

			for(int i=0;i<d;i++)nonExtreme[i]=true;
			for(int i=lowerIndex;i<upperIndex;i++)if(inLineality[i])nonExtreme[i]=true;

			int remaining=0;
			for(int i=lowerIndex;i<upperIndex;i++)if(nonExtreme[i]==false)remaining++;

			static int i;
	if(0)		if(!(i++&65535))			{
				int lowerBoundOnExtreme=this->getDimension()-this->getDimensionOfLinealitySpace();
				int candidates=0;
				for(int i=lowerIndex;i<upperIndex;i++)if((!inLineality[i])&&(!nonExtreme[i]))candidates++;
				std::cerr<<"Ambient "<<getAmbientDimension()<<
					    "\tDimension"<<getDimension()<<
					    "\tLineality"<<getDimensionOfLinealitySpace()<<
					    "\tCandidates"<<candidates<<
						"\tLowerBound"<<lowerBoundOnExtreme<<
						"\tW"<<this->originalMatrix.getWidth()<<
						toString()<<"\n"<<std::endl;
			}

	if(0)
	{
			static FrequencyTable widthHist("Width");
			widthHist.record(combinedMatrix.getWidth());
			static FrequencyTable heightHist("Height");
			heightHist.record(combinedMatrix.getHeight());
			static FrequencyTable remainingHist("Remaining");
			remainingHist.record(remaining);
			static FrequencyTable neededHist("Needed");
			neededHist.record(neededExtremeRays);
			int nInvalidColumns=0;
			for(int i=0;i<combinedMatrix.getWidth();i++)
				if(Tableau<mvtyp>::isColumnInvalid[i])nInvalidColumns++;

			int numberOfNonNegativeRows=0;
			auto norestriction=inLineality;
			static FrequencyTable numberOfNonNegativeRowsHist("nonnegativerows");
			for(int i=0;i<combinedMatrix.getHeight()-ignoreLastRow;i++)
			{
				bool isOk=true;
				for(int j=0;j<combinedMatrix.getWidth();j++)
					if(!norestriction[j])
						if(determinantOfBasis.sign()*combinedMatrix[i][j].sign()<0)
							isOk=false;
				if(isOk)
					{
						numberOfNonNegativeRows++;
						//
	//					for()
	//					if(.isPositive())norestriction[]=true;
	//					i=-1;continue;
					}
			}
			numberOfNonNegativeRowsHist.record(numberOfNonNegativeRows);

		}

			for(int i=lowerIndex;i<upperIndex;i++)
				if(!nonExtreme[i])
			{
				//we wish to check if ith column is extreme
//				for(int j=0;j<inBasis.size();j++)std::cerr<<inBasis[j];
	//				computeRowBounds();//!!!!!!!!
				if((neededExtremeRays>=remaining) || inBasis[i] && !loose(i,lowerIndex,upperIndex,nonExtreme))//we first must remove it from the basis to check
				{
					nonExtreme[i]=false;//not needed
					neededExtremeRays--;
					remaining--;
				}
				else
				{
					bool foundSimplex=towards(i,lowerIndex,upperIndex,nonExtreme)==-1;//<<"\n";
					if(foundSimplex)
						{
							nonExtreme[i]=true;
							makeColumnInvalid(i);
							remaining--;
						}
					else
					{
						remaining--;
						neededExtremeRays--;
					}
				}
			}
if(0)			assert(neededExtremeRays<=0);
if(0)			{
				int nextreme=0;
				for(int i=lowerIndex;i<upperIndex;i++)if(nonExtreme[i]==false)nextreme++;
				static FrequencyTable neededHist("nExtreme");
				neededHist.record(nextreme);
			}
		}
	public:
		/*
		 * Mark a set of extreme rays among all generators. (The reason for returning a vector with entries for all generators (rathar than just the ray generators is to make the interface cleaner).
		 */
		pmrvector<bool> markingOfExtremeRaysAmongAllGenerators()
		{
			ensureStateAsMinimum(3);
//			assert(this->getWidth()==ambientDimension+1+numberOfOriginalRayGenerators+numberOfOriginalLineGenerators);
			pmrvector<bool> ret(originalMatrix.getWidth());
			for(int i=0;i<ret.size();i++)
				ret[i]=!nonExtreme[i+ambientDimension+1];
			return ret;
		}
		void setSpecialColumnTransformed(Vector<mvtyp> const &v)  //MAY OVERFLOW  //This function also sets rows marked to be ignored because isSpecialColumnOutsideSpan uses these values
		{
			int d=getHeight();

			for(int i=0;i<d-1;i++)
			{
				mvtyp s=0;
				for(int j=0;j<d-1;j++)
				{
					assert(CHECKCOL(j));
					s+=v[j]*combinedMatrix[i][j];
				}
				assert(CHECKCOL(d-1));
				combinedMatrix[i][d-1]=s;
			}
			computeRowBounds();
		}
		bool isSpecialColumnOutsideSpan(int *violatingRowIndex=0)
		{
			int d=getHeight();
			for(int i=0;i<d-1;i++)
				if(inOrthogonalComplement[i])
					if(signAt(i,d-1))
					{
						if(violatingRowIndex)*violatingRowIndex=i;
						return true;
					}
			return false;
		}
/*		bool isSpecialColumnInLinealitySpace(int *violatingRowIndex=0) //disabled, but may work
		{
			assert(state>=2);
			int d=getHeight();
			for(int i=0;i<d-1;i++)
				if(inLineality[i])
					if(signAt(i,d-1))
					{
						if(violatingRowIndex)*violatingRowIndex=i;
						return false;
					}
			return true;
		}*/
		bool areColumnsDependent(int i, int j)//indices to columns of combinedMatrix  //MAY OVERFLOW
		{//Should be able to test if 2x2 determinant is non-zero
		/* Unfortunately the following beautiful piece of code needs more precision:
		      typ pp=dot(p,p);
		      typ qq=dot(q,q);
		      typ pq=dot(p,q);
		      return pq*pq==pp*qq;
		*/
			if(inLineality[i]||inLineality[j])return true;//it should be illegal to call this on columns that are made invalid ALSO columns that are known to be non-extreme
			assert(state>=2);
//            std::cerr<<"i"<<i<<"j"<<j<<"n"<<ambientDimension<<"\n";
			int n=ambientDimension;
			int a;
			for(a=0;a<n;a++)
				if(!ignoreRows[a])
					if(signAt(a,i))break;

			if(a==n)return true;
			assert(CHECKCOL(i));
			assert(CHECKCOL(j));
			if(combinedMatrix[a][j].isZero())return isColumnZero(j);
			mvtyp A=combinedMatrix[a][i];
			mvtyp B=combinedMatrix[a][j];
			for(int b=0;b<n;b++)
				if(!ignoreRows[b])
					if(determinantSign1(A,combinedMatrix[b][j],B,combinedMatrix[b][i])){return false;}
			return true;
		}
		/*
		 * While findExtremeRays marks exactly one column for each geometric ray, this routine must work differently.
		 *
		 * If we are in state 3, we could compare the given vector to v.
		 *
		 * If we are in state 2, we first "rewrite" v modulo the lineality space - to be more precise we transform v using
		 * the first square matrix. Then we should really project away the coordinates of the lineality space, but rather
		 * we just ignore them - they are marked in vector difference ignoreRows-inOrthogonalComplement.
		 * Then we check if v is in the span of the of the cone. This is done by checking the coordinates inOrthogonalComplement.
		 * Finally, we need first try to mark v and columns proportional to v to be ignored. Then if one of these columns is in the
		 * basis, we loose it. If it is impossible to loose, then we have an extreme ray. If we manage to loose it, then we need to
		 * check if v is in the cone of the remaining columns. This is done similarly to contains. If not v must be extreme, otherwise
		 * it is not.
		 *
		 * If v is not contained in the cone then false is returned.
		 * If v is contained in the lineality space then false is returned.
		 *
		 * The name of the function suggests that the given vector must be contained. The current implementation does not require this - but maybe it should.
		 *
		 */
		bool containedVectorIsExtreme(Vector<mvtyp> const &v)
		{
			ensureStateAsMinimum(2);

			int d=getHeight();

			assert(basisIndices[d-1]==d-1);
			setSpecialColumnTransformed(v);
//			std::cerr<<"isExtreme on"<<v<<"with this equal"<<toString()<<"\n";
			if(isSpecialColumnOutsideSpan())
				return false;
			if(isColumnZero(d-1))
				return false;

			auto disallowed=nonExtreme;
			// Unfortunately we need to find all dependent columns. Instead of doing it one coulumn at a time like this we could do a pivoting step.
			for(int i=d;i<getWidth();i++)
				if(!nonExtreme[i]&&areColumnsDependent(d-1,i))disallowed[i]=true;
			for(int i=0;i<d-1;i++)
				if(!ignoreRows[i])
				{
					assert(CHECKCOL(basisIndices[i]));
					while(areColumnsDependent(d-1,basisIndices[i]))
						{
							disallowed[basisIndices[i]]=true;
							if(!loose(basisIndices[i],d,getWidth(),disallowed))return true;
						}
				}
			return -1!=towards(d-1,d,getWidth(),disallowed);
		}
		bool isGenericSupportingHyperplane(Vector<mvtyp> const &v)const
		{
			assert(state>=1);
			int d=getHeight();
			auto prod=(Matrix<mvtyp>::rowVectorMatrix(v)*originalMatrix);
			for(int i=0;i<originalMatrix.getWidth();i++)
			{
				if(prod[0][i].isZero()!=inLineality[i+d])return false;
				if(prod[0][i].isNegative())return false;
			}
			return true;
		}
		bool isImplied(Vector<mvtyp> inequality)
		{
			for(int i=0;i<originalMatrix.getWidth();i++)
			{
				if(this->inLineality[i+getHeight()])
				{
					if(originalMatrix.columnIDot(i,inequality).isNonZero())return false;
				}
				else
				{
					if(originalMatrix.columnIDot(i,inequality).isNegative())return false;
				}
			}
			return true;
		}
		/*
		 * The following two methods together return lines resp. rays generating the cone.
		 * No other things are guaranteed for the returned values. Moreover, the property
		 * is only guaranteed for consecutive calls not interleaved with other calls.
		 *
		 * The generators are columns of the matrices.
		 */
		Matrix<mvtyp> getLineGenerators(MR *mr=get_default_resource())const
		{
			int d=getHeight();
			// If the lineality space has been found (state>=2) then it is guaranteed that a subset of the basis generates the lineality space
			// If not, we return the entire set of known lineality generators.
			if(state>=2)
				return originalMatrix.submatrixColumns(
						[this,d](int i)->bool{return inLineality[i+d]&&inBasis[i+d];},
						mr);
			return originalMatrix.submatrixColumns(
					[this,d](int i)->bool{return inLineality[i+d];},
					mr);

		}
		Matrix<mvtyp> getRayGenerators(MR *mr=get_default_resource())const
		{
			int d=getHeight();
			return originalMatrix.submatrixColumns(
					[this,d](int i)->bool{return (!inLineality[i+d])&&(!nonExtreme[i+d]);},
					mr);
		}

		friend GeneratedCone sum(GeneratedCone &a, GeneratedCone &b, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			Matrix<mvtyp> temp=combineLeftRight(a.getLineGenerators(mr),b.getLineGenerators(mr),mr2);
//			{static int i;i++;if(i%100000==0)std::cerr<<temp.toString()<<"\nR:"<<a.getRayGenerators(mr2)<<b.getRayGenerators(mr2)<<"\n";}
			return GeneratedCone(combineLeftRight(a.getRayGenerators(mr2),b.getRayGenerators(mr2),mr2),temp,mr,mr2);
/*			return GeneratedCone(
					combineLeftRight(a.getRayGenerators(),b.getRayGenerators()),
					combineLeftRight(a.getLineGenerators(),b.getLineGenerators())
					);*/
		}

		/**
		 * Computes the projection to the ith coordinateHyperPlane and returns the projection as a cone in that space.
		 */
		GeneratedCone withIthCoordinateProjectedAway(int i)const
		{
			return GeneratedCone(getRayGenerators().withIthColumnRemoved(i),getLineGenerators().withIthColumnRemoved(i));
		}

		void ensureStateAsMinimum(int s)
		{
			if(state<1 && s>=1)
			{
//				std::cerr<<"GOING FROM 0 -> 1\n";
				findOrthogonalComplementAndDimension();
				state=1;
//				std::cerr<<"DONE GOING FROM 0 -> 1\n";
			}

			if(state<2 && s>=2)
			{
//				std::cerr<<"GOING FROM 1 -> 2\n";
				findBasisForLinealitySpaceAndItsDimension();
				state=2;
//				std::cerr<<"DONE GOING FROM 1 -> 2\n";
			}
			if(state<3 && s>=3)
			{
//				std::cerr<<"GOING FROM 2 -> 3\n";
				findExtremeRays();
				state=3;
//				std::cerr<<"DONE GOING FROM 2 -> 3\n";
			}
		}
/*	bool containsInLinealitySpace(Vector<mvtyp> const &v) //disabled, but may work
	{
		ensureStateAsMinimum(2);
		setSpecialColumnTransformed(v);
		int violatingIndex;
		return isSpecialColumnInLinealitySpace(&violatingIndex);
	}*/
		/**
		 * Decides whether v is contained in the cone represented by *this.
		 * If not and if separatingHyperplaneNormal is not null, then
		 * *separatingHyperplaneNormal is changed to be a separating normal.
		 * *separatingHyperplaneNormal must have been initialised
		 * to a vector with the right number of coordinates.
		 */
	bool contains(Vector<mvtyp> const &v, Vector<mvtyp> *separatingHyperplaneNormal=0)
	{
		assert(separatingHyperplaneNormal==0 || separatingHyperplaneNormal->size()==getHeight()-1);
		bool returnValue;
		int d=getHeight();
		assert(v.size()==d-1);
		assert(basisIndices[d-1]==d-1);
		assert(signAt(d-1,d-1));
		ensureStateAsMinimum(1); //Only when we have made the cone full dimensional i.e. have a basis of its span do we have a chance of deciding if v is in the cone.
		assert(basisIndices[d-1]==d-1);
		assert(signAt(d-1,d-1));

		setSpecialColumnTransformed(v);
		// It is now possible that v is not in the span of the cone in which case we return false
		int violatingIndex;
		if(isSpecialColumnOutsideSpan(&violatingIndex))
		{
			//std::cerr<<this->toString()<<"\n";
			std::cerr<<"vioindex"<<violatingIndex<<"\n";
			if(separatingHyperplaneNormal)
			{
				if(signAt(violatingIndex,d-1)<0)
					for(int i=0;i<separatingHyperplaneNormal->size();i++)
						(*separatingHyperplaneNormal)[i]=combinedMatrix[violatingIndex][i];
				else
					for(int i=0;i<separatingHyperplaneNormal->size();i++)
						(*separatingHyperplaneNormal)[i]=-combinedMatrix[violatingIndex][i];//could this overflow?
			}
			returnValue=false;
			goto restore;//THIS NEEDS TO BE FIXED: THE SEPARATING HYPERPLANE MUST BE ORIENTED
		}
		// Now we go towards column d-1
		{
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			assert(signAt(d-1,d-1));
			int TW=towards(d-1,lowerIndex,upperIndex,nonExtreme);
			returnValue=(TW==-1);
			if(!signAt(d-1,d-1))
			{
				std::cerr<<toString();
				assert(0);
			}
			if(TW>=0)
				if(separatingHyperplaneNormal)
				{
					std::cerr<<"TW"<<TW<<"\n";
					std::cerr<<"H"<<separatingHyperplaneNormal->toString()<<"\n"<<d<<"\n";

					if(signAt(TW,d-1)<0)
						for(int i=0;i<separatingHyperplaneNormal->size();i++)
						(*separatingHyperplaneNormal)[i]=combinedMatrix[TW][i];
					else
					for(int i=0;i<separatingHyperplaneNormal->size();i++)
						(*separatingHyperplaneNormal)[i]=-combinedMatrix[TW][i];
				}
		}
		restore:
		// We need to put the standard basis vector back in at column d-1. Note that last entry did not change.
		assert(CHECKCOL(d-1));
		for(int i=0;i<d-1;i++)
			combinedMatrix[i][d-1]=0;
		// We do not care about last row of matrix, since we are in mode where this should be ignored:
//		combinedMatrix[d-1][d-1]-determinantOfBasis;
//		if((combinedMatrix[d-1][d-1]-determinantOfBasis).isNonZero())
//		{
//			std::cerr<<toString();
//			assert(0);
//		}
//		std::cerr<<"CONTAINS-END\n";
		return returnValue;
	}

	bool knownToContainVector(Vector<mvtyp> const &v)const  //MAY OVERFLOW
	{
		assert(0);//The assertion below should always fail
		assert(v.size()==getAmbientDimension()-1);
		mvtyp sign=-determinantOfBasis.sign();
		assert(state>=1);//necessary?
		for(int i=0;i<getHeight()-1;i++)
		{
			mvtyp s=0;
			for(int j=0;j<v.size();j++)
			{
				assert(CHECKCOL(j));
				s+=v[j]*combinedMatrix[i][j];
			}
			if(inOrthogonalComplement[i])
			{
				if(!s.isZero())return false;
			}
			else if(!ignoreRows[i])
			{
				if((s*sign).isNegative())return false;
			}//In the third case there is nothing to check
		}
		return true;
	}
	bool knownToContainVectorCommaEpsilon(Vector<mvtyp> const &v)const  //MAY OVERFLOW
	{
		assert(v.size()==getAmbientDimension()-1);
		int sign=-determinantOfBasis.sign();
		assert(state>=1);//necessary?
		for(int i=0;i<getHeight()-1;i++)
		{
			mvtyp s=0;
			for(int j=0;j<v.size();j++)
			{
				assert(CHECKCOL(j));
				s+=v[j]*combinedMatrix[i][j];
			}
			if(inOrthogonalComplement[i])
			{
				if(!s.isZero() || signAt(v.size(),i))return false;
			}
			else if(!ignoreRows[i])
			{
				if((s.sign()==-sign) ||s.isZero() && (signAt(v.size(),i)==-sign))return false;
			}//In the third case there is nothing to check
		}
		return true;
	}
		int getAmbientDimension()const
		{
			return getHeight()-1;
		}
		int getDimension()
		{
			ensureStateAsMinimum(1);
			return dimension;
		}
		Vector<mvtyp> getGenericSupportingHyperplane(MR *mr=get_default_resource())
		{
			ensureStateAsMinimum(2);
			int d=getHeight();
			//CHECKCOL NEEDED on cols 0 to d-1
			return combinedMatrix.subRowVector(d-1,0,d-1,mr);
		}
		/*
		 * Returns a basis of the orthogonal complement as the rows of a matrix.
		 */
		Matrix<mvtyp> getOrthogonalComplement()
		{
			ensureStateAsMinimum(1);

			Matrix<mvtyp> ret(ambientDimension-dimension,ambientDimension);
			int I=0;
			for(int i=0;i<ambientDimension;i++)
				if(inOrthogonalComplement[i])
					//CHECKCOL NEEDED on cols 0 to ambiendDimension-1
					ret[I++]=combinedMatrix.submatrix(i,0,i+1,ambientDimension)[0];

			return ret;
		}


		int getDimensionOfLinealitySpace()
		{
			ensureStateAsMinimum(2);

			int ret=0;
			for(int i=0;i<getHeight()-1;i++)
				if(ignoreRows[i]&&!inOrthogonalComplement[i])ret++;

			return ret;
		}
		/*
		 * Returns a basis of the lineality space of the cone as the columns of a matrix.
		 */
		Matrix<mvtyp> getLinealitySpace(MR *mr=get_default_resource())
		{
			ensureStateAsMinimum(2);
			Matrix<mvtyp> ret(ambientDimension,getDimensionOfLinealitySpace(),mr);
			int I=0;
			for(int i=0;i<getHeight()-1;i++)
				if(ignoreRows[i]&&!inOrthogonalComplement[i])
				{
					int col=basisIndices[i]-getHeight();
					assert(col>=0);
					for(int j=0;j<ambientDimension;j++)
						ret.UNCHECKEDACCESS(j,I)=originalMatrix.UNCHECKEDACCESS(j,col);
					I++;
				}
			return ret;
		}
		Matrix<mvtyp> getOrthogonalComplementOfLinealitySpace(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			ensureStateAsMinimum(2);
			int d=getHeight();
			Matrix<mvtyp> ret(getAmbientDimension()-getDimensionOfLinealitySpace(),ambientDimension,mr);

			int I=0;
			for(int i=0;i<getHeight()-1;i++)
				if(!(ignoreRows[i]&&!inOrthogonalComplement[i]))
				{
					//CHECKCOL NEEDED on cols 0 to d-1
					ret[I++]=combinedMatrix.submatrix(i,0,i+1,d-1,mr2)[0].toVector(mr2);
				}
			assert(I==ret.getHeight());
			return ret;
		}
		int getNumberOfRays()
		{
			ensureStateAsMinimum(3);
			int ret=0;
			for(int i=getHeight();i<getWidth();i++)
				if(!nonExtreme[i])ret++;
			return ret;
		}
		/*
		 * Returns a Matrix whose columns are the rays of the cone. (Maybe the function should return rows instead since that would avoid a transpose at the end.)//FIX implementation
		 */
		Matrix<mvtyp> getRays(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			ensureStateAsMinimum(3);
			Matrix<mvtyp> ret(getNumberOfRays(),ambientDimension,mr2);
			int I=0;
			for(int i=getHeight();i<getWidth();i++)
				if(!nonExtreme[i])
					ret[I++]=originalMatrix.column(i-getHeight(),mr2).subvector(0,ambientDimension,mr2);
			return ret.transposed(mr);
		}
	private:
		/*
		 * Below follow routines for facet enumeration using reverse search.
		 */
		/*
		 * Assuming that the current simplicial cones is part of the regular lexicographic triangulation.
		 * For each index, the current simplicial cone has a facet which may be a facet of an adjacent
		 * simplicial cone in the triangulation. That adjacent simplicial cone involves one more vector.
		 * Only the "smallest" vector can form such simplicial cone. This is a routine for comparing the
		 * candidate columns a and b.
		 */
		bool compareCol(int violatedIndex, int a, int b)//true if b is better  //MAY OVERFLOW
		{
			assert(a!=b);
			int temp=basisIndices[violatedIndex];
			exchange(violatedIndex,b);
			bool ret=doesIthsLiftsHigherThanSimplicialCone(a);
			exchange(violatedIndex,a);
			bool ret2=!doesIthsLiftsHigherThanSimplicialCone(b);
			exchange(violatedIndex,temp);
			assert(ret==ret2);
			return ret;
		}
		/*
		 * Among the column vectors on the right side of the facet opposite of violatedIndex we
		 * find the one that compare smallest in the order given by the function above.
		 */
		int findOpposite(int violatedIndex,
				int lowerIndex,
				int upperIndex,
				pmrvector<bool> const &ignoreRows,
				pmrvector<bool> &nonExtreme)
		{
			int k=violatedIndex;
			vector<int> candidates(getWidth());
			int nCandidates=0;
			for(int i=lowerIndex;i<upperIndex;i++)
				if(!nonExtreme[i])
				{
					assert(CHECKCOL(i));
					if(combinedMatrix[k][i].sign()*determinantOfBasis.sign()==-1)
						{
							candidates[i]=true;
							nCandidates++;
						}
				}
//			std::cerr<<"CANDIDATES:"<<vectorToString(candidates)<<nCandidates<<"\n";
//			std::cerr<<"CANDIDATES:";
//			for(int i=0;i<candidates.size();i++)std::cerr<<candidates[i];std::cerr<<"\n";
			assert(nCandidates>0);
			int best;
			{
				for(best=lowerIndex;best<upperIndex;best++)if(candidates[best])break;
				for(int i=best+1;i<upperIndex;i++)
					if(candidates[i]&&compareCol(violatedIndex,best,i))best=i;
			}
			return best;
		}
		/*
		 * Determines whether the facet of the current simplicial cone defined by excluding the ith basis member is
		 * in-going in the reverse search tree.
		 */
		bool isInFacet(int i)
		{
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			for(int j=lowerIndex;j<upperIndex;j++)
				if(signAt(i,j))
				{
					assert(CHECKCOL(j));
					return combinedMatrix[i][j].isNegative()==determinantOfBasis.isNegative();
				}
			return false;
		}
		/*
		 * Determines whether the basis members (with the ith excluded) generate a hyperplane defining a facet of the
		 * cone represented by the class.
		 */
		bool definesFacetOfCone(int i)
		{
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			for(int j=lowerIndex;j<upperIndex;j++)
				if(!nonExtreme[j])
			{
				assert(CHECKCOL(j));
				if(combinedMatrix[i][j].sign()+determinantOfBasis.sign()==0)
					return false;
			}
			return true;
		}

		void checkThatIthsLiftsHigherThanSimplicialCone(int i)
		{
			assert(!inBasis[i]);
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
//			for(int j=lowerIndex;j<upperIndex;j++)
			for(int j=upperIndex-1;j>=lowerIndex;j--)
			{
				if(inBasis[j])
				{
					int J=0;
					for(;J<getHeight()-1;J++)if(basisIndices[J]==j)break;
					assert(J!=getHeight()-1);
					if(!ignoreRows[J])
					{
						assert(CHECKCOL(i));
						int t=(determinantOfBasis.sign()*combinedMatrix[J][i].sign());
						if(t>0)
						{
							std::cerr<<this->toString();
							std::cerr<<"Check fails at "<<i<<" with j="<<j<<" J="<<J<<"\n";
						}
						assert(t<=0);//checks that ith lifts higher
						if(t<0)break;
					}
				}
				if(j==i) //nothing to check since ith is lifted positively
				{
					break; //If this is reached we are ok
				}
			}
		}

		void checkThatAllLiftHigher()
		{
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			for(int i=lowerIndex;i<upperIndex;i++)
				if(!nonExtreme[i])
					if(!inBasis[i])checkThatIthsLiftsHigherThanSimplicialCone(i);
		}

		/*
		 * A facet normal of the current simplicial cone may appear for other cones in the triangulation with facets
		 * spanning the same hyperplane. This routine is a test for whether the normal should be returned at the
		 * current simplicial cone. It uses that the greedy simplicial cone of the vectors in the hyperplane
		 * always exists in a lexicographic triangulation. Therefore the normal is only output if the triangulation
		 * if the simplicial cone restricted to the hyperplane is the greedy simplicial cone.
		 */
		bool lexMaxFacetRepresentative(int i)
		{
			// Because of the lifting we know that the greedy simplex exists in the triangulation of the facet
//			std::cerr<<"i"<<i<<"\n";
//			return true;//CHANGE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			//assert(definesFacetOfCone(i));
			// We now need to check if there are other 0 entries in the ith row than the ones in our basis
			// If so, we need to check if there is a non-zero coordinate of that is not
			// If so it is not a lex min representative
			for(int I=0;I<getHeight()-1;I++)
				if(!inLineality[basisIndices[I]])
				if(I!=i)
			for(int j=lowerIndex;j<basisIndices[I];j++)
				if(!nonExtreme[j])
				if(!inBasis[j])
				{
					assert(CHECKCOL(j));
					if(combinedMatrix[i][j].isZero())
					{
						if(!combinedMatrix[I][j].isZero())
							return false;
/*						for(int k=I;k<getHeight()-1;k++)
							if(k!=i)
							if(!ignoreRows[k])
								if(!combinedMatrix[k][j].isZero())
//									if(basisIndices[k]>i)
								{std::cerr<<"I"<<I<<"j"<<j<<"k"<<k<<"\n";	return false;}
*/					}
				}
			return true;
		}
		/*
		 * Returns a string that describes the incident edges in the current vertex of the reverse search tree.
		 */
		string facetStatus()
		{
			stringstream s;
			for(int i=0;i<getHeight()-1;i++)
				if(!ignoreRows[i])
			{
				if(!definesFacetOfCone(i))
					s<<" "<<i<<":"<<(isInFacet(i)?"In ":"Out ");
				else
					s<<" "<<i<<":"<<"NoFlip ";
			}
			s<<"\n";
//			std::cerr<<"DETERMINANT"<<determinantOfBasis.toString()<<"\n";
			return s.str();
		}
		/*
		 * If we assume that the current simplicial cone is int the triangulation, this routine will tell whether the ith
		 * edge at the current simplicial cone in the flip graph is present in the reverse search tree.
		 */
		bool isInTree(int i)
		{
			bool ret=false;
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
//			std::cerr<<this->toString()<<"\n";
			auto temp=basisIndices;
			checkThatAllLiftHigher();
//			std::cerr<<"BASISINDICESBEFORE"<<toString2(basisIndices)<<"\n";
			{
				int newColumnIndex=findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme);
				exchange(i,newColumnIndex);
				checkThatAllLiftHigher();
			}
//			std::cerr<<facetStatus()<<"\n";
//			std::cerr<<"flipped:"<<toString();
			for(int II=lowerIndex;II<upperIndex;II++)
				if(inBasis[II])
				{				//int I=basisIndices[II];
					int I=0;while(basisIndices[I]!=II)I++;
			//for(int I=0;I<getHeight()-1;I++)  ///It seems running through the facets in this order will not work. Follow matrix columns instead.
				if(!ignoreRows[I])
				if(!definesFacetOfCone(I))
					if(!isInFacet(I))
					{
//						std::cerr<<"Setting return value"<<I<<i<<"\n";
						ret=(I==i);
						break;
					}
				}
			{
				int newColumnIndex=findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme);

				exchange(i,newColumnIndex);
				checkThatAllLiftHigher();
			}
//			std::cerr<<"BASISINDICESAFTER "<<toString2(basisIndices)<<"\n";
//			std::cerr<<this->toString()<<"\n";
			assert(basisIndices==temp);
			return ret;
		}
		/*
		 * Since we are using lifts (M^1,M^2,....) for M big, one simplicial cone that is guaranteed
		 * to be in the triangulation is the one obtained be the greedy algorithm (in the matroid sense)
		 * for finding a basis of the column space.
		 */
		void greedySimplex()
		{
//			std::cerr<<toString()<<"\n";

			int lowerIndex=getHeight();
			int upperIndex=getWidth();

			int pos=lowerIndex;
			int unused=getHeight()-1;
			vector<int> used(unused);
//			std::cerr<<unused<<"\n";
//			std::cerr<<this->toString();
			for(int i=0;i<used.size();i++)if(ignoreRows[i]){used[i]=true;unused--;}
			while(unused>0)
			{

	//			std::cerr<<"A\n"<<vectorToString(used)<<unused<<"\n";
				if(!nonExtreme[pos])
				{
		//			std::cerr<<"B"<<pos<<"\n";
					bool isNonZero=false;
					int i;
					assert(CHECKCOL(pos));
					for(i=0;i<used.size();i++)
						if(!used[i])if(combinedMatrix[i][pos].isNonZero()){isNonZero=true;break;}
			//		std::cerr<<"C\n";
					if(isNonZero)
					{
						exchange(i,pos);
//						std::cerr<<toString()<<"\n";
						used[i]=true;
						unused--;
	//					std::cerr<<"POS:"<<pos<<"\n";
					}
				}
				pos++;
			}
		}
		/*
		 * This routine lets the current simplicial cone be the root of the reverse search tree.
		 */
		void goToTop()
		{
//			std::cerr<<"GOTOTOP\n";
//			std::cout<<"GREEDY\n";
			greedySimplex();
			checkThatAllLiftHigher();
//			std::cout<<"DONEGREEDY\n";
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
			int i;
			do
			{
				i=0;
				for(;i<getHeight()-1;i++)
				{
					if(!ignoreRows[i])///ADDED Feb 20 2020. The reverse search code was not written to work in general with lineality space and orthorgonal complement. To make it work it seems that there should be tests like this around the code.
					if(!definesFacetOfCone(i))
						if(!isInFacet(i))
						{
//							std::cerr<<"At:"<<i<<"\n";
							exchange(i,findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme));
							break;
						}
				}
			}
			while(i!=getHeight()-1);
//			std::cerr<<"DONE GOTOTOP\n";
		}
		/*
		 * This routine traverses the reverse search subtree starting at the current cone.
		 * If the pointer is non-zero, the normals associated to the subtree are appended to
		 * the rows of the matrix. The width of the matrix must match the dimension of the ambient space.
		 */
		void recursiveTriangulation(Matrix<mvtyp> *collectionOfNormals)//  MAY OVERFLOW
		{
			checkThatAllLiftHigher();
//			static int i;std::cerr<<"REK:"<<i<<"\n";i++;
//			std::cerr<<this->toString()<<"FA"<<facetStatus()<<"\n";
			int lowerIndex=getHeight();
			int upperIndex=getWidth();
//			std::cerr<<facetStatus();
			for(int i=0;i<getHeight()-1;i++)
				if(!ignoreRows[i])
				if(definesFacetOfCone(i))
				{
					//CHECKCOL NEEDED
					temporaryNormal=combinedMatrix.submatrix(i,0,i+1,getHeight()-1)[0];
					if(determinantOfBasis.sign()<0)for(int i=0;i<temporaryNormal.size();i++)temporaryNormal[i]=-temporaryNormal[i];
				//	std::cerr<<"CHECKING:"<<temporaryNormal.toString()<<this->toString()<<"\n";
					if(lexMaxFacetRepresentative(i))
			{
					//	std::cerr<<"C:TRUE\n";

				temporaryNormal=combinedMatrix.submatrix(i,0,i+1,getHeight()-1)[0];
				if(determinantOfBasis.sign()<0)for(int i=0;i<temporaryNormal.size();i++)temporaryNormal[i]=-temporaryNormal[i];
				if(collectionOfNormals)collectionOfNormals->appendRow(temporaryNormal);
			}
				//	else
				//		std::cerr<<"C:FALSE\n";
				}
//			std::cerr<<this->toString()<<"FC"<<facetStatus()<<"\n";
			for(int i=0;i<getHeight()-1;i++)
				{
				if(!ignoreRows[i])
				if(!definesFacetOfCone(i))
					if(isInFacet(i))
					{
						if(isInTree(i))
						{
//							std::cerr<<"Facet status"<<facetStatus()<<"\n";
//							std::cerr<<"exchanging"<<i<<"\n";
							int oldEntry=basisIndices[i];
//							std::cerr<<"New"<<findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme)<<"Old"<<basisIndices[i]<<"\n";
							checkThatAllLiftHigher();
							exchange(i,findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme));
							checkThatAllLiftHigher();
							recursiveTriangulation(collectionOfNormals);
//							std::cerr<<"Facet status"<<facetStatus()<<"\n";
//							std::cerr<<"exchanging"<<i<<"back:\n";
//							std::cerr<<"New"<<findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme)<<"Old"<<basisIndices[i]<<"\n";
							assert(findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme)==oldEntry);
							checkThatAllLiftHigher();
							exchange(i,findOpposite(i,lowerIndex,upperIndex,ignoreRows,nonExtreme));
							checkThatAllLiftHigher();
	//						std::cerr<<this->toString()<<"FB"<<facetStatus()<<"\n";
						}
	//					std::cerr<<"FD"<<facetStatus()<<"\n";
					}
	//			std::cerr<<"FE"<<facetStatus()<<"\n";
				}
	//		std::cerr<<"At return\n"<<this->toString()<<"Facet status"<<facetStatus()<<"\n";
		}
	public:
		/**
		 * Returns the facet normals as the rows of a matrix.
		 */
		Matrix<mvtyp> getFacetNormals()//test routine for finding information about a cone given by inequalities.
		{
//			std::cout<<"GETFACETNORMALS1\n";
//			std::cout<<"orig:"<<this->originalMatrix<<"\n";
			ensureStateAsMinimum(3);
//			std::cout<<"GETFACETNORMALS2\n";

			/*
			 * The facets should be found by reverse search on the boundary facets of the cone.
			 * This will be dual to what lrs is doing.
			 * It seems that the easiest way of doing this is by using a pulling (reverse lex) triangulation
			 * of the extreme rays of the cone. Every maximal cone in this triangulation will contribute with
			 * at least one cone in the triangulation of the boundary. From such a normal vector can be extracted.
			 * Of course, since the cone may not have simplicial facets only, each normal vector may be computed more
			 * than once. This method seems to be equal to the lrs approach, although it seems that Avis considers
			 * lexicographic triangulations.
			 * In [Rambau, DeLoera, Santos], there is a different notion of lexicographic triangulation that mixes pushing and pulling.
			 * We could consider pulling for a single vector and pushing for the rest, if this turns out to give fewer simplices.
			 *
			 * ALTERNATIVES to the above would be to use the double description method (or beneath-beyond). However, such
			 * mehods seem to be non-recursive, memory consuming and not suitable for parallelisation
			 *  - It is unclear how to make the double description method memory-less.
			 */

			/*
			 * There are several different ways that the library could return the facets to the user.
			 * Besides the obvious incidence vs coordinate question, there is also the issue of whether
			 * an iterator or a call back should be used to collect the data. We choose the simple one
			 * where a matrix of rows being the facet normals of the cone is returned.
			 */
//			std::cout<<"GETFACETNORMALS3\n"<<toString();
			goToTop();
			checkThatAllLiftHigher();
//			std::cerr<<"RECURSION CAN START\n";
//			std::cout<<"GETFACETNORMALS4\n";
			Matrix<mvtyp> facetNormals(0,getHeight()-1);
//			std::cout<<"REKBEGIN\n";
			recursiveTriangulation(&facetNormals);
//			std::cout<<"REKEND\n";
//			std::cerr<<toString()<<"\n";
//			std::cerr<<"FACET NORMALS:\n";
//			std::cout<<facetNormals.toString();
			return facetNormals;
		}
		/*
		 * Returns the face of *this selected by w. w must be in the dual cone of *this.
		 */
		GeneratedCone face(const Vector<mvtyp> &w)const
		{
			// This is not efficient. getRayGenerators for example also calls submatrixColumns with its own lambda.
			auto R=getRayGenerators();
			auto L=getLineGenerators();
			return GeneratedCone(
					R.submatrixColumns([R,w](int i)->bool{return R.columnIDot(i,w).isZero();}),
					L.submatrixColumns([L,w](int i)->bool{return L.columnIDot(i,w).isZero();})
			);
		}
		std::string resourceString()const
		{
			return memoryResourceToString(combinedMatrix.data.get_allocator().resource());
		}
	};
	static_assert(std::is_move_constructible<GeneratedCone<CircuitTableInt32>>::value,"Not move constructible!" );
	static_assert(std::is_move_assignable<GeneratedCone<CircuitTableInt32>>::value,"Not move assignable!" );
//	static_assert(std::is_nothrow_move_assignable<GeneratedCone<CircuitTableInt32>>::value,"Not nothrow move assignable!" );
//	static_assert(std::is_nothrow_move_constructible<GeneratedCone<CircuitTableInt32>>::value,"Not nothrow move constructible!" );
//	static_assert(std::is_nothrow_destructible<GeneratedCone<CircuitTableInt32>>::value,"Not nothrow destructible!" );
	template<class mvtyp> class Cone{
	public:
		GeneratedCone<mvtyp> dualCone; // The dual cone is generated by inner facet normals.
		Cone(GeneratedCone<mvtyp> const &dualCone_, MR *mr=get_default_resource()):dualCone(dualCone_,mr){}
	public:
		~Cone(){/*std::cerr<<"Cone destructor\n";*/}
		Cone(MR *mr=get_default_resource())://should be autogenerated somehow
			dualCone(mr)
		{
		}
		/**
		 * The columns of M are inequalities ((most) corresponding to _inner_facet normals) of the cone.
		 */
		Cone(Cone const &a, MR *mr=get_default_resource())://should be autogenerated somehow
			dualCone(a.dualCone,mr)
		{
		}
		Cone(Cone &&a)=default;
		Cone &operator=(Cone &&a)=default;
		Cone &operator=(Cone const &a)=default;
		Cone(Cone &&a, MR *mr):Cone(mr){/*std::cerr<<"movecons\n";*/operator=(std::move(a));}
		explicit Cone(Matrix<mvtyp> const &M):
			dualCone(M)
		{
		}
		/*
		 * The inequalities are given as columns of M, while equations are given as columns of equations.
		 */
		Cone(Matrix<mvtyp> const &M, Matrix<mvtyp> const &equations, MR *mr=get_default_resource(), MR *mr2=get_default_resource())://equations supported through duplication
			dualCone(M,equations,mr,mr2)
		{
			assert(M.getHeight()==equations.getHeight());
		}
		// Template for conversion
		template <typename otherTyp>
		explicit Cone(Cone<otherTyp> const &a, MR *mr=get_default_resource()):
		    dualCone(a.dualCone,mr)
		{
		}

		  static Cone halfSpace(Vector<mvtyp> const &inequality)
		{
			return Cone(Matrix<mvtyp>::rowVectorMatrix(inequality).transposed(),Matrix<mvtyp>(inequality.size(),0));
		}
		static Cone hyperPlane(Vector<mvtyp> const &equation)
		{
			return Cone(Matrix<mvtyp>(equation.size(),0),Matrix<mvtyp>::rowVectorMatrix(equation).transposed());
		}

		void save(std::ostream &s)const
		{
			dualCone.save(s);
		}

		static Cone load(std::istream &s)
		{
			return Cone(GeneratedCone<mvtyp>::load(s));
		}

		int getAmbientDimension()const
		{
			return dualCone.getAmbientDimension();
		}

		bool contains(Vector<mvtyp> const &v)
		{
			return dualCone.isImplied(v);
		}
		bool containsInRelativeInterior(Vector<mvtyp> const &v)
		{
			return dualCone.isGenericSupportingHyperplane(v);
		}
		bool isImplied(Vector<mvtyp> const &v)
		{
			return dualCone.contains(v);
		}
		/*
		 * The name of the function suggests that the given inequality must be supporting. The current implementation does not require this - but maybe it should.
		 */
		bool supportingInequalityDefinesFacet(Vector<mvtyp> const &v)//v can only define a facet if it is supporting (meaning that the result depends on whether -v or v was passed). If v defines the whole cone, it is unclear what the returned value is.
		{
			return dualCone.containedVectorIsExtreme(v);
		}
		bool knownThatVectorCommaEpsilonIsImplied(Vector<mvtyp> const &v)const
		{
			return dualCone.knownToContainVectorCommaEpsilon(v);
		}
		bool knownToBeImplied(Vector<mvtyp> const &v)const
		{
			assert(0);//return was missing
			return dualCone.knownToContain(v);//method does not exist
		}
		/*
		 * Returns vector of bools - one bool for each inequality/equation specifying whether that inequality defines a facet. This is done in a way that only one inequality is marked out of several inequalities defining the same facet.
		 */
		pmrvector<bool> facetDefining()
		{
			return dualCone.markingOfExtremeRaysAmongAllGenerators();
		}
		/*
		 * The returned values of the following two functions together define the cone.
		 */
		/*
		 * Returns a set of equations of the cone as COLUMNS of a Matrix.
		 */
		Matrix<mvtyp> getEquations(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			return dualCone.getLineGenerators(mr);
		}
		/*
		 * Returns a set of inequalities of the cone as rows of a Matrix.
		 */
		Matrix<mvtyp> getInequalities(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			return dualCone.getRayGenerators(mr);
		}


		Matrix<mvtyp> getLinealitySpace()
		{
			return dualCone.getOrthogonalComplement();
		}
		Vector<mvtyp> getRelativeInteriorPoint(MR *mr=get_default_resource())
		{
			return dualCone.getGenericSupportingHyperplane(mr);
		}
		int getDimension()
		{
			return getAmbientDimension()-dualCone.getDimensionOfLinealitySpace();
		}
		int getDimensionOfLinealitySpace()
		{
			return getAmbientDimension()-dualCone.getDimension();
		}
		/**
		 * Returns a basis of the orthogonal complement of the Cone as columns of a matrix.
		 */
		Matrix<mvtyp> getOrthogonalComplement(MR *mr=get_default_resource())
		{
			return dualCone.getLinealitySpace(mr);
		}
		int getNumberOfFacets()
		{
			return dualCone.getNumberOfRays();
		}
		/**
		 * Returns the facet normals of the Cone as columns of a matrix.
		 */
		Matrix<mvtyp> getFacetNormals(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			return dualCone.getRays(mr,mr2);
		}
		/**
		 * Returns the rays as rows of a matrix.
		 */
		Matrix<mvtyp> getRays()
		{
			return dualCone.getFacetNormals();
		}
		/*
		 * Returns generators of the span of the cone as the row vectors of a matrix.
		 */
		Matrix<mvtyp> getSpan(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			return dualCone.getOrthogonalComplementOfLinealitySpace(mr,mr2);
		}
		std::string toString()
		{
			stringstream ret;
			ret	    <<"Cone:{\n";
			ret		<<"Ambient dimension:"<<getAmbientDimension()<<"\n";
//			std::cerr<<"A\n";
			ret		<<"Dimension:"<<getDimension();
//			std::cerr<<"B\n";
			ret		<<"Facet Normals: (as columns)\n"<<matrixToString(getFacetNormals())<<"\n"
					<<"Rays: (as columns)\n"<<matrixToString(getRays().transposed())<<"\n";
		//			<<"Product:\n"<<matrixToString(getFacetNormals()*getRays().transposed())<<"\n"
//			std::cerr<<"C\n";
//			ret		<<"LinealitySpace:{"<<getLinealitySpace().toString()<<"}\n";
//			std::cerr<<"D\n";
			ret		<<"OrthogonalComplement:"<<matrixToString(getOrthogonalComplement())<<"\n";
//			std::cerr<<"E\n";
			ret	<<"}\n";//<<"DUAL"<<dualCone.toString();
			return ret.str();
		}
		friend Cone intersection(Cone &a, Cone &b, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			Cone ret=Cone(sum(a.dualCone,b.dualCone,mr,mr2),mr);//IT SEEMS WE SHOULD SWAP RESOURCES HERE
//			std::cerr<<"DIM:"<<a.getDimension()<<b.getDimension()<<ret.getDimension()<<"\n\n";
//			std::cerr<<"INT:"<<a.toString()<<"\n";
//			std::cerr<<" AND "<<b.toString()<<"\n";
//			std::cerr<<"="<<ret.toString()<<"\n";

			return ret;
		}
		/**
		 * Returns the embedding of the intersection of the cone with the ith coordinate hyperplane into the ith coordinate hyperplane.
		 */
		Cone intersectedWithIthCoordinateHyperplaneEmbedded(int i)const
		{
			return Cone(dualCone.withIthCoordinateProjectedAway(i));
		}
		 /*
		  * Turns inequality i into a strict inequality by subtracting inequality j.
		  * For this to work the following must be satisfied:
		  *
		  */
		void regret(int i, int j)
		{
			dualCone.regret(i,j);
		}
		 /*
		  * Turns inequality i into an equation.
		  * For this to work the following must be satisfied:
		  *
		  */
		void regretMakeEquation(int i)
		{
			dualCone.regretMakeLine(i);
		}
		std::string resourceString()const
		{
			return dualCone.resourceString();
		}
		/*
		 * Returns the link of the cone at w. w must be in the cone represented by *this.
		 */
		Cone link(const Vector<mvtyp> &w)const
		{
			return Cone(dualCone.face(w));
		}
	};
	static_assert(std::is_move_constructible<Cone<CircuitTableInt32>>::value,"Not move constructible!" );
	static_assert(std::is_move_assignable<Cone<CircuitTableInt32>>::value,"Not move assignable!" );
//	static_assert(std::is_nothrow_move_constructible<Cone<CircuitTableInt32>>::value,"Not nothrow move constructible!" );
//	static_assert(std::is_nothrow_move_assignable<Cone<CircuitTableInt32>>::value,"Not nothrow move assignable!" );
//	static_assert(std::is_nothrow_destructible<Cone<CircuitTableInt32>>::value,"Not nothrow destructible!" );

//----------------------------------------------------------------------- seperate into halfopencone template
	class Subset{
		bool usesWord;
		vector<bool> V;//used if ground set is too big
		uint64_t v;
	public:
//			Subset(Subset const & b):v(b.v){}
		Subset():v(0),usesWord(true){}
		Subset(int sizeOfGroundSet):v(0),
				usesWord(sizeOfGroundSet<=64),
				V((sizeOfGroundSet<=64)?0:sizeOfGroundSet)
		{
		}
		bool operator<(Subset const &b)const{return usesWord?v>b.v:std::lexicographical_compare(b.V.begin(),b.V.end(),V.begin(),V.end());}///NOTE THAT THIS IS OPPOSITE TO MAKE SORTING EASIER
		string toString()const//not fully implemented
		{
			return bitset<8>(v).to_string();
		}
		bool isSupersetOf(Subset const &b)const{
//			std::cerr<<"This"<<toString()<<"is superset of"<<b.toString()<<":"<<!(b.v&~v)<<"\n";
			if(usesWord)return !(b.v&~v);
			for(int i=0;i<V.size();i++)
				if(b.V[i] & !V[i])
					return false;
			return true;
		}
		bool isSubsetOf(Subset const &b)const{
//			std::cerr<<"This"<<toString()<<"is subset of"<<b.toString()<<":"<<b.isSupersetOf(*this)<<"\n";
			return b.isSupersetOf(*this);
		}
		Subset intersect(Subset const &b)const
		{
			Subset ret=b;
			if(usesWord)
				ret.v=b.v&v;
			else
				for(int i=0;i<V.size();i++)ret.V[i]=ret.V[i]&V[i];
			return ret;
		}
		void insert(int i)
		{
			if(usesWord)
				v|=((uint64_t)1)<<i;
			else
				V[i]=true;
		}
		bool contains(int i)const
		{
			if(usesWord)
				return (((uint64_t)1)<<i)&v;
			return V[i];
		}
	};
	template<class mvtyp> class RayCollector{
		Matrix<mvtyp> linealitySpaceComplement;
		map<Vector<mvtyp>,int> indexMap;
	public:
		Matrix<mvtyp> rays;
	public:
		RayCollector(Matrix<mvtyp> const &linealitySpace_):
			linealitySpaceComplement(fromZMatrix<mvtyp>(kernel(toZMatrix(linealitySpace_)))),
			rays(0,linealitySpace_.getWidth())
		{
		}
		int lookup(Vector<mvtyp> const &v)
		{
//			std::cerr<<"Looking up"<<v.toString()<<"\n";
			auto temp=(linealitySpaceComplement*Matrix<mvtyp>::rowVectorMatrix(normalize(v)).transposed()).transposed()[0].toVector();
//			std::cerr<<"translates to "<<temp.toString()<<"\n";
//			int d=gcd(temp).v;
//			for(int i=0;i<temp.size();i++)temp[i].v=temp[i].v/d;//replace with call to normalize
//			temp=temp.normalized();
			temp=normalize(temp);
//			std::cerr<<"Normalized:"<<temp.toString()<<"\n";
			if(indexMap.find(temp)==indexMap.end())
			{
				indexMap[temp]=rays.getHeight();
				rays.appendRow(v);
			}
			return indexMap[temp];
		}
		Matrix<mvtyp> getRays()const
		{
			return rays;
		}
	};

static vector<int> intersection(vector<int> const &a, vector<int> const &b)
{
  vector<int> ret;
  set_intersection(a.begin(),a.end(),b.begin(),b.end(),back_inserter(ret));
  return ret;
}
static vector<int> setDifference(vector<int> const &a, vector<int> const &b)
{
  vector<int> ret;
  set_difference(a.begin(),a.end(),b.begin(),b.end(),back_inserter(ret));
  return ret;
}
/**
 * Given two sorted vectors a and b, this function decides whether the elements of a is a subset of the elements of b.
 */
static bool isSubsetOf(vector<int> const &a, vector<int> const &b)
{
	unsigned next=0;
    for(unsigned i=0;i<a.size();i++)
      {
    	while(1)
          {
            if(next>=b.size())return false;
            if(a[i]==b[next])break;
            next++;
          }
      }
    return true;
}


	template<class mvtyp> class Complex{
	public: // PreComplex needs to access privates
		int ambientDimension;
		mutable int dimension;
		Matrix<mvtyp> linealitySpace;
		Matrix<mvtyp> linealitySpaceComplement;
		Matrix<mvtyp> vertices;
		std::map<Vector<mvtyp>,int> indexMap;
		Vector<mvtyp> uniqueCoordinates(Vector<mvtyp> const &v){auto temp=normalize((linealitySpaceComplement*Matrix<mvtyp>::rowVectorMatrix(normalize(v)).transposed()).transposed()[0].toVector());return temp;}
	public:
		int lookUpIndex(Vector<mvtyp> const &v)
		{
			return indexMap[uniqueCoordinates(v)];
		}
		   class Cone
		  {
			   bool isKnownToBeNonMaximalFlag;
		  public:
			int multiplicity;
		    vector<int> indices;//always sorted
		    int dimension;
		    mutable bool markedAsDeleted;
		    Cone(std::set<int> const &indices_, int dimension_, Complex const &complex, bool knownToBeNonMaximal=false):
				dimension(dimension_),
				multiplicity(1),
				isKnownToBeNonMaximalFlag(knownToBeNonMaximal),
				markedAsDeleted(false)
		    {
		    	for(auto &i:indices_)indices.push_back(i);
		    }
		    void markAsDeleted()const
		    {
		    	markedAsDeleted=true;
		    }
		    std::set<int> indexSet()const
		    {
		      std::set<int> ret;
		      for(auto &i:indices)
		        ret.insert(i);

		      return ret;
		    }
		    string toString()const
		    {
		    	stringstream s;
		    	s<<"{";
		    	bool first=true;
		    	for(auto &i:indices)
		    		{
		    			if(!first)s<<",";
		    			s<<i;
		    			first=false;
		    		}
		    	s<<"}";
		    	return s.str();
		    }
		    bool operator<(const Cone & b)const
		    {
		    	return indices<b.indices;
		    }
		    bool isKnownToBeNonMaximal()const{return isKnownToBeNonMaximalFlag;}
		    bool isSubsetOf(Cone const &c)const
		    {
		    	return gfan::isSubsetOf(indices,c.indices);
		    }
//		    Vector<mvtyp> relativeInteriorPoint()const
//				{}
/*		    GeneratedCone<mvtyp> toGeneratedCone(Complex const & complex)const
				{
		    		Matrix<mvtyp> rays(complex.vertices.getWidth(),indices.size());
		    		for(int i=0;i<indices.size();i++)
		    			rays[i]=complex.vertices[indices[i]];
		    		return GeneratedCone<mvtyp>(rays.transposed(),complex.linealitySpace.transposed());
				}*/
		  };
		   typedef std::set<Cone> ConeContainer;
		   ConeContainer cones;
		   Complex(Matrix<mvtyp> const &rays, Matrix<mvtyp> const &linealitySpace_):
			   vertices(rays),
			   linealitySpace(linealitySpace_),
			   ambientDimension(rays.getWidth()),
			   dimension(-1),
			   linealitySpaceComplement(fromZMatrix<mvtyp>(kernel(toZMatrix(linealitySpace_))))
		   {
//			   std::cerr<<"LIN SPACE CoMP"<<linealitySpaceComplement<<"\n";
//			   std::cerr<<"Constructing with"<<rays.toString()<<linealitySpace_.toString()<<"\n";
			   for(int i=0;i<vertices.getHeight();i++)indexMap[uniqueCoordinates(vertices[i].toVector())]=i;
		   }
	  vector<vector<Cone const*> > rayIncidences;//for each ray, store one vector
	  void buildRayIncidences()
	  {
	    rayIncidences=vector<vector<Cone const*> >(vertices.getHeight());
	    for(auto &c:cones)
	      for(auto &a:c.indices)
		rayIncidences[a].push_back(&c);
	  }
	  vector<vector<int>> facesContainedInSubset(vector<int> const &a)const
		{
		  vector<vector<int>> ret;
		  for(auto &c:cones)
			  if(!c.markedAsDeleted)
				  if(isSubsetOf(c.indices,a))
					  ret.push_back(c.indices);
		  return ret;
		}
		bool insert(Cone const &c)
		{
	        if(c.dimension>dimension)dimension=c.dimension;
			return cones.insert(c).second;
//			std::cerr<<"Inserting "<<c.toString()<<"\n";
		}
		/**
		 * If faces have been deleted, it may be necessary to set recompute to true.
		 */
		int getMaxDim(bool recompute=false)const
		{
			if(recompute)
			{
				dimension=-1;
				for(auto &i:cones)
					if(!i.markedAsDeleted)
						if(i.dimension>dimension)dimension=i.dimension;
			}
		  return dimension;
		}
		int getMinDim()const//why is this needed?
		{
			int ret=100000;
			for(auto &i:cones)
				if(i.dimension<ret)ret=i.dimension;
			return ret;
		}
		bool isMaximal(Cone const &c)const
		{
		  if(c.isKnownToBeNonMaximal())return false;
		  if(c.dimension==dimension)return true;
///		  for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
//		    {
//		      Cone c2=c.permuted(*k,*this,false);
		      for(auto i=cones.begin();i!=cones.end();i++)
		        {
		          if(i->dimension>c.dimension)
		            if(c.isSubsetOf(*i) && !i->isSubsetOf(c))return false;
		        }
//		    }
		  return true;
		}
/*		bool intersectNotAtInfinity(Cone const &a, Cone const &b, int infinityCoordinate)
		{
			auto A=a.indexSet();
			auto B=b.indexSet();
			for(auto x:A)
				if(B.count(x))
					if(!vertices[x][infinityCoordinate].isZero())
						return true;
			return false;
		}*/
/*		void curveCleaning()
		{
			for(auto c:cones)
				if(isMaximal(c))
				{
					auto S=c.toGeneratedCone().getRecession(infinityCoordinate);
					for(auto f:ConeContainer)
						if(f.indices!=c.indices)
						{
							auto cap=intersection(c,f);
							auto relInt=cap.getRelativeInteriorPoint();
							auto F=f.toCone();
							S=S+F.link(relInt);
						}
					if(S.isPointed())
						std::cerr<<"This relatively open polyhedron can be ignored.\n";
					else
						std::cerr<<"This relatively open polyhedron cannot be ignored.\n";
				}
		}*/
		std::string toStringJustCones(int dimLow, int dimHigh, bool onlyMaximal, bool group, std::ostream *multiplicities=0, bool compressed=false, bool tPlaneSort=false)const
		{
			  std::stringstream ret;

			  ZVector additionalSortKeys(cones.size());

			  for(int d=dimLow;d<=dimHigh;d++)
			    {
			      int numberOfOrbitsOutput=0;
			      int numberOfOrbitsOfThisDimension=0;
			      bool newDimension=true;
			        {
			          int I=0;
			          for(typename ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++,I++)
			        	  if(!i->markedAsDeleted)
			                  if(i->dimension==d)
			                    {
			                      numberOfOrbitsOfThisDimension++;
			              if(!onlyMaximal || isMaximal(*i))
			                {
			                  numberOfOrbitsOutput++;
			                  bool isMax=true;/*isMaximal(*i);*/
			                  bool newOrbit=true;
			                  std::set<std::set<int> > temp;
/*			                    for(SymmetryGroup::ElementContainer::const_iterator k=sym.elements.begin();k!=sym.elements.end();k++)
			                      {
			                        Cone temp1=i->permuted(*k,*this,false);
			                        temp.insert(temp1.indexSet());
			                        if(compressed)break;
			                    }*/
			                  temp.insert(i->indexSet());
			                  for(std::set<std::set<int> >::const_iterator j=temp.begin();j!=temp.end();j++)
			                    {
			                      ret << "{";
			                      for(std::set<int>::const_iterator a=j->begin();a!=j->end();a++)
			                        {
			                          if(a!=j->begin())ret<<" ";
			                          ret << *a;
			                        }
			                      ret << "}";
			                      if(group)if(newOrbit)ret << "\t# New orbit";
			                      if(newDimension)ret << "\t# Dimension "<<d;
			                      ret <<std::endl;
			                      if(isMax)if(multiplicities)
			                        {
			                          *multiplicities << i->multiplicity;
			                          if(group)if(newOrbit)*multiplicities << "\t# New orbit";
			                          if(newDimension)*multiplicities << "\t# Dimension "<<d;
			                          *multiplicities << std::endl;
			                        }
			                      newOrbit=false;
			                      newDimension=false;
			                    }
			              }
			                    }
			        }
			    }

			  return ret.str();
		}
		Vector<int> fvector(bool boundedPart=false, bool infinityPart=false)const
		{
			auto dimHigh=getMaxDim();
			auto dimLow=getMinDim();
			assert(dimLow<=dimHigh);
			Vector<int> ret(dimHigh-dimLow+1);
			for(typename ConeContainer::const_iterator i=cones.begin();i!=cones.end();i++)
				if(!i->markedAsDeleted)
		    {
		      bool doAdd=true;
		      if(boundedPart)
		        {
		    	  doAdd=!boundedPart;
		          bool isBounded=true;
		          for(unsigned j=0;j<i->indices.size();j++)
		            if(vertices[i->indices[j]][0].sign()==0)isBounded=false;
		          doAdd=isBounded;
		        }
		      else
		    	  if(infinityPart)
		    	  {
		    		  bool isAtInfinity=true;
		    		  for(unsigned j=0;j<i->indices.size();j++)
		    			  if(vertices[i->indices[j]][0].sign()!=0)isAtInfinity=false;
		    		  doAdd=isAtInfinity;
		    	  }
		      if(i->dimension>dimHigh || i->dimension<dimLow)
		      {
		    	  std::cerr<<"In fvector computation for complex:\n";
		    	  std::cerr<<"A cone of dimension "<<i->dimension<<" was found but the complex was supposed to have cones of dimension "<<dimLow<<" to "<<dimHigh<<"\n";
		    	  assert(0);
		      }
		      if(doAdd)
		    	  ret[i->dimension-dimLow]++;
		    }
			return ret;
		}
		std::string toString(int flags=0)const
		{
		    PolymakeFile polymakeFile;
		    polymakeFile.create("NONAME","PolyhedralFan","PolyhedralFan",flags&FPF_xml);
		    polymakeFile.writeCardinalProperty("AMBIENT_DIM",ambientDimension);
//		    polymakeFile.writeCardinalProperty("DIM",getMaxDim());
		    polymakeFile.writeCardinalProperty("LINEALITY_DIM",linealitySpace.getHeight());
//		    if(flags&&FPF_boundedInfo)polymakeFile.writeMatrixProperty("RAYS",toZMatrix(vertices),true);//??????
		    if(flags&&FPF_rays)polymakeFile.writeMatrixProperty("RAYS",toZMatrix(vertices),true);
		    polymakeFile.writeCardinalProperty("N_RAYS",vertices.getHeight());

		    polymakeFile.writeMatrixProperty("LINEALITY_SPACE",toZMatrix(linealitySpace),ambientDimension);
		    polymakeFile.writeMatrixProperty("ORTH_LINEALITY_SPACE",kernel(toZMatrix(linealitySpace)),ambientDimension);

		    ZVector fvector=toZVector(this->fvector());
		    polymakeFile.writeCardinalVectorProperty("F_VECTOR",fvector);

//		    if(flags&FPF_boundedInfo)
//		      {
//		        ZVector fvectorBounded=this->fvector(true);
//		        polymakeFile.writeCardinalVectorProperty("F_VECTOR_BOUNDED",fvectorBounded);
//		      }
//		    polymakeFile.writeCardinalProperty("SIMPLICIAL",isSimplicial());
//		    polymakeFile.writeCardinalProperty("PURE",isPure());

		    if(flags&FPF_cones)polymakeFile.writeStringProperty("CONES",toStringJustCones(getMinDim(),getMaxDim(),false,flags&FPF_group, 0,false,flags&FPF_tPlaneSort));
		    if(flags&FPF_maximalCones)polymakeFile.writeStringProperty("MAXIMAL_CONES_OF_CLOSURE",toStringJustCones(getMinDim(),getMaxDim(),true,flags&FPF_group, 0,false,flags&FPF_tPlaneSort));
/*		    if(flags&FPF_conesCompressed)polymakeFile.writeStringProperty("CONES_ORBITS",toStringJustCones(getMinDim(),getMaxDim(),false,flags&FPF_group, 0,true,flags&FPF_tPlaneSort));
		    if((flags&FPF_conesCompressed) && (flags&FPF_maximalCones))polymakeFile.writeStringProperty("MAXIMAL_CONES_ORBITS",toStringJustCones(getMinDim(),getMaxDim(),true,flags&FPF_group, 0,true,flags&FPF_tPlaneSort));

		    if(!sym.isTrivial())
		      {
		        polymakeFile.writeMatrixProperty("SYMMETRY_GENERATORS",IntToZMatrix(sym.getGenerators()));
		      }
*/
		    std::stringstream s;
		    polymakeFile.writeStream(s);
		    return s.str();
		}
		void multipleVertexCheck()
		{
			for(int i=0;i<vertices.getHeight();i++)
				for(int j=0;j<i;j++)
					if(vertices[i].toVector()==vertices[j].toVector())
						{
							std::cerr<<"Row "<<i<<" and "<<j<<" are the same in vertex matrix:"<<vertices;
							assert(0);
						}
		}
		void substitute(vector<Cone*> oldCones, vector<pair<gfan::Matrix<mvtyp>,int>> &newCones)
		{
//			multipleVertexCheck();
			//look up which rays are not already in complex
			Matrix<mvtyp> newVertices(0,vertices.getWidth());
			for(auto &A:newCones)
				for(int i=0;i<A.first.getHeight();i++)
					if(indexMap.count(uniqueCoordinates(A.first[i].toVector()))==0)
						newVertices.appendRow(/*uniqueCoordinates*/(A.first[i].toVector()));
			newVertices.normalizeRows();
			newVertices.sortAndRemoveDuplicateRows();

			vertices=combineOnTop(vertices,newVertices);
//			std::cerr<<"New Vertices:\n"<<vertices.toString()<<"\n";
			for(int i=vertices.getHeight()-newVertices.getHeight();i<vertices.getHeight();i++)
				indexMap[uniqueCoordinates(vertices[i].toVector())]=i;

			for(auto r:oldCones)
				r->markAsDeleted();

			//std::erase_if(cones,[](auto const &c){return c.markedAsDeleted;});
			for (auto it{cones.begin()}, end{cones.end()}; it != end; ) {
			  if (it->markedAsDeleted){
			    it = cones.erase(it);
			  } else ++it;
			}
			for(auto &A:newCones)
				{
					std::set<int> indices;
//					std::cerr<<"RAYS OF CONE TO BE ADDED:\n"<<A.first.toString();
					for(int i=0;i<A.first.getHeight();i++)
						indices.insert(indexMap[uniqueCoordinates(A.first[i].toVector())]);
					auto temp=Cone(indices,A.second/*dimension*/,*this,false);
//					std::cerr<<"Inserting:"<<temp.toString()<<"\n";
					cones.insert(temp);
				}
			buildRayIncidences();
			// maybe we should recompute dimension of complex
//			multipleVertexCheck();
		}
	};
	/*
	 * There can be many rays in a fan, while individual cones may involve just a few rays.
	 * During face lattice extraction the rays involved in a cone are marked with a bit.
	 * To keep the number of required bits low, the each bit corresponds to a vector in a subcollection of rays.
	 * The purpose of the IndexTranslator is to setup such a subcollection of rays.
	 */
	template<class mvtyp> class IndexTranslator
	{
		Complex<mvtyp> &c;
		Matrix<mvtyp> rays;
	public:
		IndexTranslator(Complex<mvtyp> &c_, Matrix<mvtyp> const &rays_):
			c(c_),
			rays(rays_)
		{
		}
		Subset raysToSubset(Matrix<mvtyp> const &coneRays)
		{
			Subset ret(rays.getHeight());
			for(int i=0;i<coneRays.getHeight();i++)
				ret.insert(c.lookUpIndex(coneRays[i]));
			return ret;
		}
		typename Complex<mvtyp>::Cone subsetToCone(Subset const &s, int dim, bool knownToBeNonMaximal)const
		{
			set<int> indices;
			int S=0;
			for(int i=0;i<rays.getHeight();i++)
				if(s.contains(i))
				{
					S++;
					indices.insert(c.lookUpIndex(rays[i]));
				}
			assert(indices.size()==S);
			return typename Complex<mvtyp>::Cone(indices,dim,c,knownToBeNonMaximal);
		}
		//returns true if inserted and false if s existed already
		bool insert(Subset const &s, int dimension, bool knownToBeNonMaximal)
		{
			return c.insert(subsetToCone(s,dimension,knownToBeNonMaximal));
		}
	};
	static string toString(vector<Subset> const &V)
	{
		stringstream s;
		for(auto &a:V)s<<a.toString()<<"\n";
		return s.str();
	}

/**
 * To represent a HalfOpenCone, a Cone in higher dimension is constructed - the lifted cone.
 * The new coordinate is the _last_ coordinate in the new space.
 * For points in the lifted cone this coordinate is _negative_.
 * The represented half open cone is the projection of the intersection of the lifted cone
 * and the open negative halfspace.
 */
	template<class mvtyp> class HalfOpenCone{
		template<typename> friend class HalfOpenCone;
public:		Cone<mvtyp> lifted; //The strict inequality gets 1 as the additional coordinate, meaning that the halfopen cone is the projection of the intersection of the lifted cone with the strictly negative halfspace.
private:
		int specialColumnIndex;//column containing (0,0,0,0,1).
		HalfOpenCone(Cone<mvtyp> const &lifted_, MR *mr=get_default_resource()):lifted(lifted_,mr){/*std::cerr<<"copy\n";*/specialColumnIndex=-1;}
	public:
//			~HalfOpenCone()NOEXCEPT_{/*std::cerr<<"HalfOpenCone destructor\n";*/}
		HalfOpenCone(MR *mr=get_default_resource()):lifted(mr){/*std::cerr<<"copy\n";*/specialColumnIndex=-1;}
		HalfOpenCone(const HalfOpenCone &c, MR *mr):lifted(c.lifted,mr),specialColumnIndex(c.specialColumnIndex){/*std::cerr<<"copy\n";*/}
//		HalfOpenCone(const HalfOpenCone &c):lifted(c.lifted),specialColumnIndex(c.specialColumnIndex){/*std::cerr<<"copy\n";*/}
//		HalfOpenCone(HalfOpenCone<mvtyp> &&c)NOEXCEPT_:lifted(std::move(c.lifted)),specialColumnIndex(c.specialColumnIndex){/*std::cerr<<"move\n";*/}
//		HalfOpenCone(HalfOpenCone<mvtyp> &&c, MR *mr)NOEXCEPT_:lifted(std::move(c.lifted),((StackResource*)mr)/*->dump()*/),specialColumnIndex(c.specialColumnIndex){/*std::cerr<<"move construct\n";*/}
//		HalfOpenCone &operator=(HalfOpenCone<mvtyp> &&c)NOEXCEPT_{/*std::cerr<<"move assign\n";*/lifted=std::move(c.lifted);specialColumnIndex=c.specialColumnIndex;return *this;}
//		HalfOpenCone &operator=(HalfOpenCone<mvtyp> const &c){/*std::cerr<<"assign\n";*/lifted=c.lifted;specialColumnIndex=c.specialColumnIndex;return *this;}
//		HalfOpenCone()//only here so that std::vector can be used on this.???
//		{
//		}
		// We append an additional inequality to make a better correspondence between facets of the two cones
/*
 * inequalities:
 * NNSS 0
 * NNSS 0
 * 0011-1
 * equations:
 * EE
 * EE
 * 00
 */
		/** The (in)equalities are columns of the given matrices. */
		HalfOpenCone(Matrix<mvtyp> const &nonStrict, Matrix<mvtyp> const &equations, Matrix<mvtyp> const &strict, MR *mr=get_default_resource(), MR *mr2=get_default_resource()):
		lifted(
#if 1
				mr
#else
				combineLeftRight(//HERE
						combineOnTop(nonStrict,Matrix<mvtyp>::filled(1,nonStrict.getWidth(),mvtyp(0)),mr2),
						combineLeftRight(
								combineOnTop(strict,Matrix<mvtyp>::filled(1,strict.getWidth(),mvtyp(1)),mr2),
								combineOnTop(Matrix<mvtyp>(nonStrict.getHeight(),1,mr2),Matrix<mvtyp>::filled(1,1,mvtyp(-1)),mr2),
								mr2),
							mr2),
				combineOnTop(equations,Matrix<mvtyp>::filled(1,equations.getWidth(),mvtyp(0)),mr2),mr//,mr2//TO HERE
#endif
				),
				specialColumnIndex(nonStrict.getWidth()+strict.getWidth())
		{
#if 1
			Matrix<mvtyp> nonStrictM(nonStrict.getHeight()+1,nonStrict.getWidth()+strict.getWidth()+1,mr2);
			nonStrictM.setSubMatrix(0,0,nonStrict.getHeight(),nonStrict.getWidth(),nonStrict);
			nonStrictM.setSubMatrix(0,nonStrict.getWidth(),strict.getHeight(),nonStrict.getWidth()+strict.getWidth(),strict);
			for(int i=0;i<strict.getWidth();i++)nonStrictM[nonStrict.getHeight()][nonStrict.getWidth()+i]=mvtyp(1);
			nonStrictM[nonStrict.getHeight()][nonStrict.getWidth()+strict.getWidth()]=mvtyp(-1);
			Matrix<mvtyp> equationsM(equations.getHeight()+1,equations.getWidth(),mr2);
			equationsM.setSubMatrix(0,0,equations.getHeight(),equations.getWidth(),equations);

			lifted=Cone<mvtyp>(nonStrictM,equationsM,mr,mr2);
#endif
//			std::cerr<<"HalfOpenCone construtor on"<<" "<<nonStrict.addressOfData()<<equations.addressOfData()<<" "<<strict.addressOfData()<<"\n";
			assert(nonStrict.getHeight()==strict.getHeight());
//			std::cerr<<"location of lifted data"<<lifted.dualCone.combinedMatrix.addressOfData()<<"\n";
//			std::cerr<<"Constructor:"<<
//					combineOnTop(combineLeftRight(nonStrict,Matrix<mvtyp>::filled(nonStrict.getHeight(),1,mvtyp(0))),
//							combineLeftRight(strict,Matrix<mvtyp>::filled(strict.getHeight(),1,mvtyp(1)))).toString()<<"\n";
		}

		// Template for conversion
		template <typename otherTyp>
		explicit HalfOpenCone(HalfOpenCone<otherTyp> const &a, MR *mr=get_default_resource()):
		    lifted(a.lifted,mr),
			specialColumnIndex(a.specialColumnIndex)
		{
//			std::cerr<<"HalfOpenCone conversion\n";
		}


		static HalfOpenCone openHalfSpace(Vector<mvtyp> const &normal, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			return HalfOpenCone(
					Matrix<mvtyp>(normal.size(),0,mr2),
					Matrix<mvtyp>(normal.size(),0,mr2),
					Matrix<mvtyp>::rowVectorMatrix(normal,mr2).transposed(mr2),
					mr,
					mr2);
		}
		static HalfOpenCone hyperplane(Vector<mvtyp> const &normal, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			return HalfOpenCone(
					Matrix<mvtyp>(normal.size(),0,mr2),
					Matrix<mvtyp>::rowVectorMatrix(normal,mr2).transposed(mr2),
					Matrix<mvtyp>(normal.size(),0,mr2),
					mr,
					mr2);
		}

		void save(std::ostream &s)const
		{
			lifted.save(s);
		}

		static HalfOpenCone load(std::istream &s)
		{
			return HalfOpenCone(Cone<mvtyp>::load(s));
		}
		/*
		 * If the cone was constructed by giving nonStrict, equations and strict, this method returns
		 * a vector of bools indexing columns of nonStrict equation and strict and marking those that define facets.
		 * THIS DOCUMENTATION SEEMS WRONG.
		 */
		pmrvector<bool> facetDefining()
		{
			return lifted.facetDefining();
		}
		bool isEmpty(MR *mr=get_default_resource()){
//			if(lifted.dualCone.getState()>=2)
//				lifted.dualCone.containsInLinealitySpace(Vector<mvtyp>::standardVector(lifted.getAmbientDimension(),lifted.getAmbientDimension()-1,mr));
			return lifted.isImplied(Vector<mvtyp>::standardVector(lifted.getAmbientDimension(),lifted.getAmbientDimension()-1,mr));
		}
		Cone<mvtyp> closure(MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			assert(!isEmpty(mr));
			int n=lifted.getAmbientDimension()-1;
//			std::cerr<<"BLA"<<(lifted.getFacetNormals().toString()+lifted.getOrthogonalComplement().toString())<<"BLA\n";

/*			std::cerr<<"taking closure of"<<lifted.toString()<<
					lifted.getFacetNormals().submatrixColumns(0,n)<<
					lifted.getOrthogonalComplement().submatrixColumns(0,n)<<"\n";
*/
			return Cone<mvtyp>(
					lifted.getFacetNormals(mr2,mr).submatrixRows(0,n,mr2),//HERE
					lifted.getOrthogonalComplement(mr2).submatrixRows(0,n,mr2),mr,mr2//TO HERE
			);
		}
		bool containsOrigin()
		{
			return lifted.contains(-Vector<mvtyp>::standardVector(lifted.getAmbientDimension(),lifted.getAmbientDimension()-1));
		}
		/**
		 * If a HalfOpenCone is nonempty, then a given vector v is contained in its closure iff (v,0) is contained in the lift.
		 */
		bool closureOfNonEmptyConeContains(Vector<mvtyp> const &v, MR *mr=get_default_resource())
		{
			return lifted.contains(concatenation(v,Vector<mvtyp>(1,mr),mr));
		}
		/**
		 * The projection of a relative interior point of the lifted cone is a relative interior point of the cone.
		 * It is of course a mistake to ask for a relative interior point of an empty set.
		 */
		Vector<mvtyp> getRelativeInteriorPoint(MR *mr=get_default_resource(),MR *mr2=get_default_resource())
		{
			return lifted.getRelativeInteriorPoint(mr2).subvector(0,lifted.getAmbientDimension()-1,mr);
		}

		friend HalfOpenCone intersection2(HalfOpenCone &a, HalfOpenCone &b, MR *mr=get_default_resource(), MR *mr2=get_default_resource())
		{
			HalfOpenCone ret=HalfOpenCone(intersection(a.lifted,b.lifted,mr,mr2),mr);//IT SEEMS RESOURCES SHOULD BE SWAPPED
//			ret.lifted.dualCone.ensureStateAsMinimum(3);
//			std::cerr<<"INTERSECTION2 "<<a.containsOrigin()<<b.containsOrigin()<<ret.containsOrigin()<<"\n";
			return ret;
		}
		int getDimension(MR *mr=get_default_resource()){
			if(isEmpty(mr))return -1;
			return lifted.getDimension()-1;
		}
		int getAmbientDimension()const
		{
			return lifted.getAmbientDimension()-1;
		}
		std::string toString()
		{
			return std::string(
					"<halfopencone\n")+
					(
							isEmpty()?
									std::string("EMPTY")
			:
									std::string()+(containsOrigin()?"ZERO CONTAINED":"ZERO NOT CONTAINED")+"\n--------CLOSURE:\n"+closure().toString()+"\n---------LIFTED:\n"+lifted.toString()
					)
					+std::string("\nhalfopencone>");
		}
		/*
		 * The name of the function suggests that the given inequality must be valid. The current implementation does not require this - but maybe it should.
		 */
		bool validInequalityDefinesExistingFacet(Vector<mvtyp> const &v, MR *mr2=get_default_resource())
		{
			bool returnValue=lifted.supportingInequalityDefinesFacet(concatenation(v,Vector<mvtyp>(1),mr2));//Is this really the correct test?

//			std::cerr<<"definesExistingFacetOn:"<<v.toString()<<"and"<<toString()<<"return value"<<returnValue<<"\n";
			return returnValue;
		}
		bool knownToBeContainedInOpenHalfSpace(Vector<mvtyp> const &normal)const
		{
			return lifted.knownThatVectorCommaEpsilonIsImplied(normal);
		}
		bool knownToBeContainedInClosedHalfSpace(Vector<mvtyp> const &normal)const
		{
			assert(0);//DEAD CODE
			return lifted.knownToBeImplied(concatenation(normal,Vector<mvtyp>(1)));
		}
		bool knownToNotIntersectOpenHalfSpace(Vector<mvtyp> const &normal)const
		{
			assert(0);//DEAD CODE
			return knownToBeContainedInClosedHalfSpace(-normal);
		}
		bool knownToNotIntersectClosedHalfSpace(Vector<mvtyp> const &normal)const
		{
			return knownToBeContainedInOpenHalfSpace(-normal);
		}
		bool knownToNotIntersectHyperplane(Vector<mvtyp> const &normal)const
		{
			return knownToBeContainedInOpenHalfSpace(normal) || knownToBeContainedInOpenHalfSpace(-normal);
		}
//		bool knownToBeContainedInHyperplane It should be possible to make this exact if the lineality space of the dual cone is known.


		 /*
		  * Turns a ith nonstrict inequality into an equation.
		  * For this to work the following must be satisfied:
		  *
		  */
		void makeStrictInequality(int i)
		{
			assert(specialColumnIndex!=-1);
			lifted.regret(i,specialColumnIndex);
		}
		 /*
		  * Turns the ith nonstrict inequality into an equation.
		  * For this to work the following must be satisfied:
		  *
		  */
		void makeEquation(int i)
		{
			lifted.regretMakeEquation(i);
		}
		Vector<int> fVector()//this function is slow
		{
			int n=lifted.getAmbientDimension()-1;
			if(isEmpty())return Vector<int>(n+1);
			auto C=closure();

			auto D=C.getFacetNormals();
			for(int i=0;i<D.getWidth();i++)
				if(this->validInequalityDefinesExistingFacet(D.column(i)))
				{
					auto H=hyperplane(D.column(i));
					auto A=intersection2(*this,H);
					auto O=openHalfSpace(D.column(i));
					auto B=intersection2(*this,O);
					return A.fVector()+B.fVector();
				}
			return Vector<int>::standardVector(n+1,lifted.getDimension()-1);
		}
		/* Finding facets in halfopencone:
		 * Assume extreme rays have been computed
		 * Assume candidate facet normals given
		 * Assume excluding inequalities given
		 * Compute indices to extreme rays
		 * Iterate over candidat facets:
		 * 		find subset of supported rays at facet
		 * 		rank computation could show if this defines facet
		 * 		check if one excluding inequality excludes all supported rays. If so, ignore this facet
		 * 	Alternative to rank computation is inclusion tests. Note that only non-excluded faces needs to be tested
		 *
		 */


//		void forEachRelativelyOpenFace(std::function<bool (HalfOpenCone &)> func)
//		{
//		}

		static string toString(vector<Subset> const &V)///appears other  places too
		{
			stringstream s;
			for(auto &a:V)s<<a.toString()<<"\n";
			return s.str();
		}
		void addFacesToComplex(/*Complex<mvtyp> &complex*/IndexTranslator<mvtyp> &complex, Subset const &rayIndices, vector<Subset> const &inequalitySupports, vector<Subset> &excludingSubsets, int dimension, bool knownToBeNonMaximal=false, bool onlyThis=false)
		{
//			std::cerr<<"called with\nRays:"<<rayIndices.toString()<<"\nIneq:\n"<<toString(inequalitySupports)<<"Excl:\n"<<toString(excludingSubsets)<<"\nDimension:"<<dimension<<"\n";
			if(not complex.insert(rayIndices,dimension,knownToBeNonMaximal) or onlyThis) return;//Note that in this case lower non-excluded faces might not get inserted. Therefore there are two modes for using the Complex. Either only "disjoint halfopen cones" are inserted or only cones with no excluding faces.
			vector<Subset> facetCandidates;
			for(auto &f:inequalitySupports)
				if(!rayIndices.isSubsetOf(f))//Ignore inequalities that do not remove anything
				{
					auto face=f.intersect(rayIndices);//This looks like a mistake: why exclude when finding facets? No, this is ok, because if a facet that shows that a face is non-maximal is excluded, then so is the face
					if(!std::any_of(excludingSubsets.begin(),excludingSubsets.end(),[&](Subset &a){return a.isSupersetOf(face);}))
						facetCandidates.push_back(face);
				}
			sort(facetCandidates.begin(),facetCandidates.end());//make sure that supersets appear last
//			std::cerr<<"FacetCandidates:\n"<<toString(facetCandidates)<<"\n";
			//check for inclusions
			vector<Subset> trueFacets;
			for(auto i=facetCandidates.begin();i!=facetCandidates.end();i++)
				if(!std::any_of(facetCandidates.begin(),i,[&](Subset &a){return i->isSubsetOf(a);}))
					trueFacets.push_back(*i);
//			std::cerr<<"TrueFacets:\n"<<toString(trueFacets)<<"\n";
			for(auto &f:trueFacets)
			{
			  addFacesToComplex(complex,f,trueFacets,excludingSubsets,dimension-1,true);
				excludingSubsets.push_back(f);//Exclude this facet in other calls
			}
			excludingSubsets.resize(excludingSubsets.size()-trueFacets.size());//undo changes to excludingSubsets
		}
		bool representedPolyhedronIsBounded(int specialCoordinate=0)
		{
			if(isEmpty())return true;
			auto closed=this->closure();
			if(0)
			{
				auto rays=closed.getRays();
				auto lines=closed.getLinealitySpace();
				for(int i=0;i<rays.getHeight();i++)if(rays[i][0].isZero())return false;
				for(int i=0;i<lines.getHeight();i++)if(lines[i][0].isZero())return false;
				return true;
			}
			auto temp1=Cone<mvtyp>::halfSpace(
					Vector<mvtyp>::standardVector(
							closed.getAmbientDimension(),
							specialCoordinate
							)
					);
			auto temp=intersection(
					closed,
					temp1);
			return temp.getDimension()==0;
		}
		bool representedPolyhedronsRecessionConeIsContainedInNegativeHalfSpace(int specialCoordinate=0)//except the origin
		{
			auto closed=this->closure();
			auto rays=closed.getRays();
			auto lines=closed.getLinealitySpace();

			for(int i=0;i<rays.getHeight();i++)
				if(rays[i][0].isZero())
					{
						Vector<mvtyp> v=rays[i].toVector();
						if(!v.sum().isNegative())return false;
					}
			for(int i=0;i<lines.getHeight();i++)
				if(lines[i][0].isZero())
					if(!lines[i].toVector().sum().isZero())return false;
			return true;
		}
		/*
		 * This method should maybe be deleted.
		 * It definitely can be improved, since the recession cone of a halfopen polyhedron is a halfopen cone.
		 */
		Cone<mvtyp> representedPolyhedronsRecessionConeInHigherDim(int specialCoordinate=0)
		{
			assert(!isEmpty());
			auto closed=closure();
			auto H=Cone<mvtyp>::hyperPlane(
					Vector<mvtyp>::standardVector(
							closed.getAmbientDimension(),
							specialCoordinate
							)
					);
			return intersection(
					closed,
					H);
		}
	};
	static_assert(std::is_move_constructible<HalfOpenCone<CircuitTableInt32>>::value,"Not move constructible!" );
	static_assert(std::is_move_assignable<HalfOpenCone<CircuitTableInt32>>::value,"Not move assignable!" );
//	static_assert(std::is_nothrow_move_constructible<HalfOpenCone<CircuitTableInt32>>::value,"Not nothrow move constructible!" );
//	static_assert(std::is_nothrow_move_assignable<HalfOpenCone<CircuitTableInt32>>::value,"Not nothrow move assignable!" );
//	static_assert(std::is_nothrow_destructible<HalfOpenCone<CircuitTableInt32>>::value,"Not nothrow destructible!" );

	template<typename mvtyp>
	  //static string 
	  Complex<mvtyp>
	  extractFaceComplex(vector<HalfOpenCone<mvtyp>> &cones, RayCollector<mvtyp> &collector, Matrix<mvtyp> const &linealitySpace, int fpf=FPF_cones|FPF_maximalCones, bool takeClosure=false, bool onlyMax=false, int logLevel=0)//make non-static
	{
		Complex<mvtyp> C(collector.rays,linealitySpace);
		if(logLevel>=1)std::cerr<<"EXTRACTING FACE COMPLEX. NCONES="<<cones.size()<<"\n";
		int count = 0;
		for(auto &c:cones)
		{
			if(logLevel>=1)if (count++ % 1000 == 0) std::cerr << cones.size() - count+1 << " cones left.\n";
			auto closedCone=c.closure();
			int dimension=closedCone.getDimension();
			auto inequalities=closedCone.getFacetNormals().transposed();
			auto localRays=closedCone.getRays();

//			std::cerr<<"DIMENSION:"<<dimension<<"\n";
//			std::cerr<<"DIMENSION:"<<closedCone.getdimension<<"\n";
//			std::cerr<<"FACETNORMALS:"<<inequalities.toString();
//			std::cerr<<"LOCALRAYS"<<localRays<<"\n";

/*			if(localRays.getHeight()>64){
				std::cerr<<"SKIPPING CONE WITH MORE THAN 64 RAYS\n";
				continue;
			}*/
			IndexTranslator<mvtyp> translator(C,localRays);
			Subset rayIndices(localRays.getHeight());
			for(int i=0;i<localRays.getHeight();i++)rayIndices.insert(i);
			vector<Subset> excludingSubsets;
			auto exclusionsWithExtraCoordinates=c.lifted.getFacetNormals().transposed();
			if(!takeClosure)
				for(int i=0;i<exclusionsWithExtraCoordinates.getHeight();i++)
					if(exclusionsWithExtraCoordinates[i][exclusionsWithExtraCoordinates.getWidth()-1].isPositive())
					{
						Subset s(localRays.getHeight());
						for(int j=0;j<localRays.getHeight();j++)
						{
							if(localRays[j].dot(exclusionsWithExtraCoordinates.submatrixColumns(0,exclusionsWithExtraCoordinates.getWidth()-1)[i]).isZero())s.insert(j);
						}
						excludingSubsets.push_back(s);
					}
			vector<Subset> inequalitiesSupports;
			for(int i=0;i<inequalities.getHeight();i++)
			{
				Subset s(localRays.getHeight());
				for(int j=0;j<localRays.getHeight();j++)
				{
					if(localRays[j].dot(inequalities[i]).isZero())s.insert(j);
				}
				inequalitiesSupports.push_back(s);
			}
			c.addFacesToComplex(translator,rayIndices,inequalitiesSupports,excludingSubsets,dimension,false/*not known to be max*/,onlyMax/*add only c*/);
		}
		if (not onlyMax) {
//		  std::cerr<<"Fvector from complex"<<C.fvector().toString()<<"\n";//This vector should actually be part of the returned string.
//		  std::cerr<<"Infinity Fvector from complex (only makes sense if first coordinate is non-negative/non-positive)"<<C.fvector(false,true).toString()<<"\n";
		}
		return C;
		//return C.toString(fpf);
	}
	template<typename mvtyp>
	static string extractFaceComplex(vector<HalfOpenCone<mvtyp>> &cones)//make non-static
	{
		for(auto &F:cones)
			if(!F.isEmpty())
			{
				auto closure=F.closure();
				Matrix<mvtyp> linealitySpace=closure.getLinealitySpace();
				RayCollector<mvtyp> collector(linealitySpace);
				for(auto &c:cones)
					if(!c.isEmpty())
					{
						auto rays=c.closure().getRays();
						for(int i=0;i<rays.getHeight();i++)
							collector.lookup(normalize(rays[i].toVector()));//duplicate code?
					}
				return extractFaceComplex(cones,collector,linealitySpace).toString(FPF_cones|FPF_maximalCones);
			}
		return string("Empty\n");
	}
	template<typename typ>
	static string extractFaceComplex(HalfOpenCone<typ> &cone)//make non-static
	{
		vector<HalfOpenCone<typ>> C(1,cone);
		return extractFaceComplex(C);
	}

#if 0
	void allocatortest()
	 {
		 auto mr=&stackResource;
		 auto mr2=&stackResource2;

		 Matrix<CircuitTableInt32> nonstrict(3,3,mr);
		 Matrix<CircuitTableInt32> strict(3,3,mr);
		 HalfOpenCone<CircuitTableInt32> cone(strict.transposed(),nonstrict.transposed(),nonstrict.transposed(),mr);
//		 assert(0);
		 vector<HalfOpenCone<CircuitTableInt32>> ret;
		 std::cerr<<"At start\n";
		 mr2->dump();
		 for(int i=0;i<10;i++)
		 {
			 std::cerr<<"Iteration"<<i<<"\n";
//			 ResourceInvariant a(mr2);
			 {
				 mr2->dump();
				 std::cerr<<"Reserving\n";
			 ret.reserve(ret.size()+1);

			 mr2->dump();
			 std::cerr<<"Building A\n";
			 auto A=HalfOpenCone<CircuitTableInt32>(nonstrict.transposed(),Matrix<CircuitTableInt32>(0,nonstrict.getWidth()).transposed(),strict.transposed(),mr2,mr);
			 std::cerr<<"Done building A\n";
			 mr2->dump();

			 ret.emplace_back(std::move(intersection2(A,cone,mr2)),mr);//FAILS IF ENABLED
			 }
//			 a.assertIsInvariant();
		 }
//		 mr->dump();
		 mr2->dump();
		 assert(0);
	 }
#endif

	template<class mvtyp>
	std::vector<HalfOpenCone<mvtyp>> loadHalfOpenConeVector(std::istream &s)
	{
		auto a=parseSequence(s);
		std::vector<HalfOpenCone<mvtyp>> ret;
		ret.reserve(a.size());
		for(auto &s:a)
		{
			stringstream S(s);
			ret.push_back(HalfOpenCone<mvtyp>::load(S));
		}
		return ret;
	}
	template<class mvtyp>
	void save(std::vector<HalfOpenCone<mvtyp>> const &v,std::ostream &s)
	{
		s<<"(";
		for(int i=0;i<v.size();i++)
		{
			if(i)
				s<<",\n";
			v[i].save(s);
		}
		s<<")\n";
	}
};


#endif /* GFANLIB_TABLEAU_H_ */
