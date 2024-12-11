/*
 * lib_zmatrix.h
 *
 *  Created on: Sep 28, 2010
 *      Author: anders
 */

#ifndef LIB_ZMATRIX_H_
#define LIB_ZMATRIX_H_

#include <vector>
#include <algorithm>
#include <functional>
#include "gfanlib_vector.h"

#include "gfanlib_memoryresource.h"

#include<type_traits>
#include<iostream>
#include<iomanip>

// These definitions come from gist.github.com  maddouri/has_member.hpp
#define define_has_member(member_name) \
template<typename T> \
class has_member_##member_name \
{ \
	typedef char yes_type; \
	typedef long no_type; \
	template<typename U>static yes_type test(decltype(&U::member_name)); \
	template<typename U>static no_type test(...); \
public: \
	static constexpr bool value=sizeof(test<T>(0))==sizeof(yes_type); \
};
define_has_member(POD)
define_has_member(POD2)
#define has_member(class_,member_name) has_member_##member_name<class_>::value

namespace gfan{

//taken from wikibooks.org:


/*#define GENERATE_HAS_MEMBER_TYPE(Type)
template <class T>
class HasMemberType_##Type
{
private:
/*	using Yes=char[2];
	using No=char[1];
	struct Fallback {struct Type{};};
	struct Derived:T,Fallback{};
	template<class U>
	static No& test (typename U::Type*);
	template<typename U>
	static Yes& test (U*);
public:
	static constexpr bool RESULT=sizeof(test<Derived>(nullptr))==sizeof(Yes);
};
*/
//template <typename... values> class Matrix;
//template <class typ, auto> class Matrix;

//template <class typ, class POD=typ::POD>

// Here POD is supposed to be an easy-to-initialise, easy-to-move/copy version of typ (if possible).
//	template <class typ, class POD>

	template<bool selector, typename FirstType, typename SecondType>
	struct IfElse
	{
		typedef SecondType type;
	};
	template<typename FirstType, typename SecondType>
	struct IfElse<true,FirstType,SecondType>
	{
		typedef FirstType type;
	};

/*	template <typename typ,
			typename POD=typename IfElse<
				has_member(typ,POD2),
				std::enable_if<has_member(typ,POD2),typename typ::POD>,
				std::enable_if<true,typ>
>::type::type
	>*/
//	template <typename typ,auto> class Matrix;

	constexpr bool hasPod(...){return false;}
	constexpr bool hasPod(int *){return false;}
	constexpr bool hasPod(class CircuitTableInt32 *){return true;}


//	template<typename typ,typename POD=std::conditional_t<
//			false/*hasPod((typ*)(0))*/,
//			std::enable_if<false,typename typ::POD>,
//			typ>>
//	template <typename typ,typename POD>

	namespace aux{	template<typename T>
		struct element_type{
			typedef typename T::POD/*::element_type*/ type;
		};}
	namespace aux{	template<typename T>
		struct element_type2{
			typedef  T/*::element_type*/ type;
		};}
	template <typename typ>
class Matrix{
		typedef typename std::conditional<hasPod((typ*)(0)),aux::element_type<typ>,aux::element_type2<typ>>::type func_;
		typedef typename func_::type POD;
//	template <class typ, class POD/*=typ*/> class Matrix{
  int width,height;
//  std::vector<Vector<typ> > rows;
//  std::vector<typ> data;
public:  pmrvector<POD> data;
public:
  // rowIterator;
 // std::vector<Vector<typ> >::iterator rowsBegin(){return rows.begin();}
//  std::vector<Vector<typ> >::iterator rowsEnd(){return rows.end();}
  inline int getHeight()const{return height;};
  inline int getWidth()const{return width;};
  Matrix(Matrix &&a)noexcept:width(a.width),height(a.height),data(std::move(a.data)){};
  Matrix &operator=(Matrix &&a)noexcept{/*std::cerr<<"Matrix move assignment\n";*/width=a.width;height=a.height;data=std::move(a.data);return *this;};
  Matrix &operator=(Matrix const &a)=default;
  Matrix(Matrix &&a, MR *mr)noexcept:data(mr){operator=(std::move(a));}
  Matrix(const Matrix &a, MR *mr=get_default_resource()):width(a.getWidth()),height(a.getHeight()),data(a.data,mr){
  }
  Matrix(int height_, int width_):width(width_),height(height_),data(width_*height_){
    assert(height>=0);
    assert(width>=0);
  };//(signed long int)0
  Matrix(int height_, int width_, MR *mr):width(width_),height(height_),data(width_*height_, POD{},mr){
    assert(height>=0);
    assert(width>=0);
  };
  ~Matrix(){
  };
  Matrix(MR *mr=get_default_resource()):width(0),height(0),data(0,POD{}/*(signed long int)0*/,mr)
  {
//	  assertIsZero();
//	  std::cerr<<"Matrix constructor\n";
  };
  // Template for conversion
  template <typename otherTyp>
  explicit Matrix(Matrix<otherTyp> const &a, MR *mr=get_default_resource()):
  	  Matrix(a.getHeight(),a.getWidth(),mr)
  {
	  for(int i=0;i<a.getHeight();i++)
		  for(int j=0;j<a.getWidth();j++)
			  (*this)[i][j]=typ(a[i][j]);
  }
  typ* addressOfData()const
  {
	  if(height*width==0)return 0;
	  return (typ*)&(data[0]);
  }
  void reserveRows(int nrows)
  {
	  data.reserve(nrows*width);
  }
  static Matrix rowVectorMatrix(Vector<typ> const &v, MR *mr=get_default_resource())
  {
    Matrix ret(1,v.size(),mr);
    for(int i=0;i<v.size();i++)ret[0][i]=v[i];
    return ret;//MOVE??
  }
  Vector<typ> column(int i, MR *mr=get_default_resource())const
    {
      assert(i>=0);
      assert(i<getWidth());
      Vector<typ> ret(getHeight(),mr);
      for(int j=0;j<getHeight();j++)ret[j]=(*this)[j][i];
      return ret;
    }
  typ columnIDot(int i, Vector<typ> const &v)const
  {
	  assert(v.size()==getHeight());
	  typ ret;
	  for(int j=0;j<v.size();j++)ret+=UNCHECKEDACCESS(j,i)*v[j];
	  return ret;
  }
  Matrix transposed(MR *mr=get_default_resource())const
    {
      Matrix ret(getWidth(),getHeight(),mr);
      for(int i=0;i<getWidth();i++)
          for(int j=0;j<getHeight();j++)
        	  ret[i][j]=(*this)[j][i];
      return ret;
    }
  static Matrix identity(int n)
    {
      Matrix m(n,n);
      for(int i=0;i<n;i++)m[i][i]=typ(1);
      return m;
    }
  static Matrix filled(int h, int w, typ v)
    {
      Matrix m(h,w);
      for(int i=0;i<h;i++)
          for(int j=0;j<w;j++)
        	  m[i][j]=v;
      return m;
    }
  void append(Matrix const &m)
    {
      if(m.getWidth()!=width)
	{
	  std::cerr<<"this:"<<height<<"x"<<width<<"\n";
	  std::cerr<<"m:"<<m.height<<"x"<<m.width<<"\n";
	}
      assert(m.getWidth()==width);
	  data.resize((height+m.height)*width);
	  int oldHeight=height;
      height+=m.height;
      for(int i=0;i<m.height;i++)
        {
          for(int j=0;j<m.width;j++)
        	  (*this)[i+oldHeight][j]=m[i][j];
        }
    }
  void appendRow(Vector<typ> const &v)
    {
	  if(v.size()!=width)
	  {
		  std::cerr<<"Appending row of size "<<v.size()<<" to a matrix of size "<<height<<"x"<<width<<"\n";
	  }
	  assert(v.size()==width);
#if 0
	  data.resize((height+1)*width); // quadratic complexity in height if repeated
	  for(int j=0;j<width;j++)
		  (*this)[height][j]=v[j];
#else
	  for(int j=0;j<width;j++)data.push_back(v[j]); // amortised linear complexity in height if repeated
#endif
	  height++;
    }
  void eraseLastRow()
  {
    assert(height>0);
    data.resize((height-1)*width);
    height--;
  }
/*  IntegerVector vectormultiply(IntegerVector const &v)const//slow
    {
      assert(v.size()==width);
      IntegerVector ret(height);
      for(int i=0;i<height;i++)
        ret[i]=dot(rows[i],v);
      return ret;
    }*/
  /**
   * Decides if v is in the kernel of the matrix.
   */
/*  bool inKernel(IntegerVector const &v)const
    {
      assert(v.size()==width);
      for(int i=0;i<height;i++)
          if(dotLong(rows[i],v)!=0)return false;
      return true;
    }
*/
  friend Matrix operator*(const typ &s, const Matrix& q)
    {
      Matrix p=q;
      for(int i=0;i<q.height;i++)
          for(int j=0;j<q.width;j++)
        	  p[i][j]=s*(q[i][j]);
      return p;
    }
  friend Matrix operator*(const Matrix& a, const Matrix& b)
    {
	  if(a.width!=b.height)
		  std::cerr<<"Trying to multiply "<<a.getHeight()<<"x"<<a.getWidth()<<" with "<<b.getHeight()<<"x"<<b.getWidth()<<" matrix!\n";
      assert(a.width==b.height);
      Matrix ret(a.height,b.width);
      for(int i=0;i<a.height;i++)
      {
    	  for(int k=0;k<b.width;k++)
    	  {
    		  typ sum(0);
    		  for(int j=0;j<a.width;j++)
    			  sum+=a[i][j]*b[j][k];
    		  ret[i][k]=sum;
    	  }
      }
      return ret;
    }
  /*  template<class T>
    Matrix<T>(const Matrix<T>& c):v(c.size()){
    for(int i=0;i<size();i++)v[i]=typ(c[i]);}
  */
  friend Matrix operator-(const Matrix &b)
  {
    Matrix ret(b.height,b.width);
    for(int i=0;i<b.height;i++)ret[i]=-b[i];
    return ret;
  }


  void copyEntriesRestrict(typ const * __restrict__ src, typ * __restrict__ dest, int n)
  {
	  for(int i=0;i<n;i++)dest[i]=src[i];
  }
  void setSubMatrix(int startRow, int startColumn, int endRow, int endColumn, Matrix const &m, int rowOffsetInM=0, int columnOffsetInM=0)
  {
	  assert(0<=startRow && startRow<=getHeight());
	  assert(0<=startColumn && startColumn<=getWidth());
	  assert(0<=endRow && endRow<=getHeight());
	  assert(0<=endColumn && endColumn<=getWidth());
	  assert(endRow-startRow<=m.getHeight());
	  assert(endColumn-startColumn<=m.getWidth());
	  if(endColumn!=startColumn)  // <- test is needed because we may otherwise be indexing out of bound
	  for(int i=startRow;i<endRow;i++)
	  {
		  auto r=m[i-startRow+rowOffsetInM];
		  auto R=(*this)[i];
		  copyEntriesRestrict(&(r[columnOffsetInM]),&(R[startColumn]),endColumn-startColumn);
//		  for(int a=startColumn;a<endColumn;a++)
//			  R.UNCHECKEDACCESS(a)=r.UNCHECKEDACCESS(a-startColumn+columnOffsetInM);
	  }
  }
  /**
     Returns the specified submatrix. The endRow and endColumn are not included.
   */
  Matrix submatrix(int startRow, int startColumn, int endRow, int endColumn, MR *mr=get_default_resource())const
  {
    assert(startRow>=0);
    assert(startColumn>=0);
    assert(endRow>=startRow);
    assert(endColumn>=startColumn);
    assert(endRow<=height);
    assert(endColumn<=width);
    Matrix ret(endRow-startRow,endColumn-startColumn,mr);
    ret.setSubMatrix(0,0,endRow-startRow,endColumn-startColumn,*this,startRow,startColumn);
//    for(int i=startRow;i<endRow;i++)
//      for(int j=startColumn;j<endColumn;j++)
//        ret[i-startRow][j-startColumn]=(*this)[i][j];
    return ret;
  }

  Matrix submatrixColumns(int startColumn, int endColumn, MR *mr=get_default_resource())const
  {
	  return submatrix(0,startColumn,getHeight(),endColumn,mr);
  }

  Matrix submatrixColumns(std::function<bool (int)> func, MR *mr=get_default_resource())const
  {
	  int width=0;for(int i=0;i<getWidth();i++)if(func(i))width++;
	  Matrix ret(getHeight(),width,mr);
	  int I=0;
	  for(int i=0;i<getWidth();i++)
		  if(func(i))
		  {
			  for(int j=0;j<getHeight();j++)
				  ret.UNCHECKEDACCESS(j,I)=(*this).UNCHECKEDACCESS(j,i);
			  I++;
		  }
	  return ret;
  }

  Matrix withIthColumnRemoved(int i)const
  {
	  return submatrixColumns([i](int a)->bool{return a!=i;});
  }

  Matrix submatrixRows(int startRow, int endRow, MR *mr=get_default_resource())const
  {
	  return submatrix(startRow,0,endRow,getWidth(),mr);
  }

  Vector<typ> subRowVector(int startRow, int startColumn, int endColumn, MR *mr=get_default_resource())const
  {
	  Vector<typ> ret(endColumn-startColumn,mr);
	  for(int i=0;i<ret.size();i++)ret[i]=(*this)[startRow][startColumn+i];
	  return ret;
  }

  class RowRef;
  class const_RowRef{
    int rowNumM;
    Matrix const &matrix;
    friend class RowRef;
  public:
  inline const_RowRef(const Matrix  &matrix_, int rowNum_)__attribute__((always_inline)):
    rowNumM(rowNum_*matrix_.width),
      matrix(matrix_)
      {
      }
  inline typ const &operator[](int j)const __attribute__((always_inline))
    {
	assert(j>=0);
	assert(j<matrix.width);
	return static_cast<typ const&>(matrix.data[rowNumM+j]);
    }
  inline typ const &UNCHECKEDACCESS(int j)const __attribute__((always_inline))
    {
	return static_cast<typ const&>(matrix.data[rowNumM+j]);
    }
    const Vector<typ> toVector(MR *mr=get_default_resource())const
    {
      Vector<typ> ret(matrix.width,mr);
      for(int j=0;j<matrix.width;j++)
    	  ret[j]=matrix.data[rowNumM+j];
      return ret;
    }
    operator Vector<typ>()const
		{
			return toVector();
		}
    bool operator==(Vector<typ> const &b)const
		{
			return toVector()==b;
		}
/*    typ dot(Vector<typ> const &b)const
		{
			return dot(toVector(),b);
		}*/
    Vector<typ> operator-()const
    {
    	return -toVector();
    }
  };
  class RowRef{
    int rowNumM;
    Matrix &matrix;
  public:
  inline RowRef(Matrix &matrix_, int rowNum_):
    rowNumM(rowNum_*matrix_.width),
      matrix(matrix_)
      {
      }
    inline typ &operator[](int j) __attribute__((always_inline))
      {
    	assert(j>=0);
    	assert(j<matrix.width);
    	return static_cast<typ&>(matrix.data[rowNumM+j]);
      }
    inline typ &UNCHECKEDACCESS(int j)
      {
    	return static_cast<typ&>(matrix.data[rowNumM+j]);
      }
    RowRef &operator=(Vector<typ> const &v)
    {
        assert(v.size()==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]=v[j];

    	return *this;
    }
    RowRef &operator=(RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
/*    RowRef &operator+=(Vector<typ> const &v)
    {
        assert(v.size()==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]+=v.v[j];

    	return *this;
    }*/
    RowRef &operator+=(RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]+=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    RowRef &operator+=(const_RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]+=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    RowRef &operator=(const_RowRef const &v)
    {
        assert(v.matrix.width==matrix.width);
        for(int j=0;j<matrix.width;j++)
        	matrix.data[rowNumM+j]=v.matrix.data[v.rowNumM+j];

    	return *this;
    }
    const Vector<typ> toVector(MR *mr=get_default_resource())const
    {
      Vector<typ> ret(matrix.width,mr);
      for(int j=0;j<matrix.width;j++)
    	  ret[j]=matrix.data[rowNumM+j];
      return ret;
    }
    operator Vector<typ>()const
		{
			return toVector();
		}
/*    typ dot(Vector<typ> const &b)const
		{
			return dot(toVector(),b);
		}*/
    typ dot(RowRef const& b)const
    {
    	assert(matrix.width==b.matrix.width);
    	typ ret((int)0);
    	for(int j=0;j<matrix.width;j++)
    		ret+=matrix.data[rowNumM+j]*b.matrix.data[b.rowNumM+j];
    	return ret;
    }
    bool isZero()const
        {
          for(int j=0;j<matrix.width;j++)if(!(matrix.data[rowNumM+j].isZero()))return false;
          return true;
        }
  };


  inline RowRef operator[](int i) __attribute__((always_inline))
  {
    assert(i>=0);
    assert(i<height);
    return RowRef(*this,i);
  }
  inline const_RowRef operator[](int i)const __attribute__((always_inline))
  {
    assert(i>=0);
    assert(i<height);
    return const_RowRef(*this,i);
  }


  const typ& UNCHECKEDACCESS(int i,int j)const __attribute__((always_inline)){
/*	    assert(i>=0);
	    assert(i<height);
	    assert(j>=0);
	    assert(j<width);*/
	  return static_cast<typ const&>(data[i*width+j]);}
  typ& UNCHECKEDACCESS(int i,int j) __attribute__((always_inline)){
/*	    assert(i>=0);
	    assert(i<height);
	    assert(j>=0);
	    assert(j<width);*/
	    return static_cast<typ&>(data[i*width+j]);}



  bool operator<(const Matrix & b)const
  {
    if(getWidth()<b.getWidth())return true;
    if(b.getWidth()<getWidth())return false;
    if(getHeight()<b.getHeight())return true;
    if(b.getHeight()<getHeight())return false;

    for(int i=0;i<getHeight();i++)
      {
        if((*this)[i].toVector()<b[i].toVector())return true;
        if(b[i].toVector()<(*this)[i].toVector())return false;
      }
    return false;
  }
  /**
     Adds a times the i th row to the j th row.
  */
  void madd(int i, typ a, int j)
  {
    assert(i!=j);
    assert(i>=0 && i<height);
    assert(j>=0 && j<height);

    if(!a.isZero())
    for(int k=0;k<width;k++)
      if(!(*this)[i][k].isZero())
    	  (*this)[j][k].madd((*this)[i][k],a);
  }

  friend std::ostream &operator<<(std::ostream &f, Matrix const &a){
    f<<"{";
    for(int i=0;i<a.getHeight();i++)
      {
        if(i)f<<","<<std::endl;
        f<<a[i].toVector();
      }
    f<<"}"<<std::endl;
    return f;
  }


  static Matrix readMatrix(std::istream &f, int width)
  {
	  auto v=parseSequence(f);
	  int height=v.size();
	  Matrix ret(height,width);
	  int i=0;
	  for(auto c:v)
	  {
	    //std::cerr<<"\""<<c<<"\""<<c.size()<<"\n";
		  if(width>0)										// We cannot tell the difference between a 0x0 and a 1x0 matrix!
			  if(c.size()==0)continue;
		  Vector<typ> v;
		  std::stringstream t(c);
		  t>>v;
		  ret[i++]=v;
	  }
	  return ret;
  }

  std::string toString()const
  {
	  std::stringstream f;
	  f<<*this;
	  return f.str();
  }

  /**
     Swaps the i th and the j th row.
   */
  void swapRows(int i, int j)
  {
    for(int a=0;a<width;a++)std::swap((*this)[i][a],(*this)[j][a]);
  }
  /**
     This method is used for iterating through the pivots in a matrix
     in row echelon form. To find the first pivot put i=-1 and
     j=-1 and call this routine. The (i,j) th entry of the matrix is a
     pivot. Call the routine again to get the next pivot. When no more
     pivots are found the routine returns false.
  */
  bool nextPivot(int &i, int &j)const
  {
    i++;
    if(i>=height)return false;
    while(++j<width)
      {
        if(!(*this)[i][j].isZero()) return true;
      }
    return false;
  }
  /**
     Returns the indices of the columns containing a pivot.
     The returned list is sorted.
     The matrix must be in row echelon form.
   */
  std::vector<int> pivotColumns()const
  {
    std::vector<int> ret;
    int pivotI=-1;
    int pivotJ=-1;
    while(nextPivot(pivotI,pivotJ))ret.push_back(pivotJ);
    return ret;
  }
  /**
     Returns the indices of the columns not containing a pivot.
     The returned list is sorted.
     The matrix must be in row echelon form.
   */
  std::vector<int> nonPivotColumns()const
  {
    std::vector<int> ret;
    int pivotI=-1;
    int pivotJ=-1;
    int firstPossiblePivot=0;
    while(nextPivot(pivotI,pivotJ))
      {
        for(int j=firstPossiblePivot;j<pivotJ;j++)
          ret.push_back(j);
        firstPossiblePivot=pivotJ+1;
      }
    for(int j=firstPossiblePivot;j<getWidth();j++)
      ret.push_back(j);

    return ret;
  }
  /**
     This routine removes the zero rows of the matrix.
   */
  void removeZeroRows()
  {
    int nonZeros=0;
    for(int i=0;i<height;i++)if(!(*this)[i].isZero())nonZeros++;
    if(nonZeros==height)return;

    Matrix b(nonZeros,width);

    int j=0;
    for(int i=0;i<height;i++)
      {
        if(!(*this)[i].isZero())
          {
            b[j]=(*this)[i];
            j++;
          }
      }
    *this=b;
  }
  /**
     Returns the index of a row whose index is at least currentRow and
     which has a non-zero element on the column th entry. If no such
     row exists then -1 is returned. This routine is used in the Gauss
     reduction. To make the reduction more efficient the routine
     chooses its row with as few non-zero entries as possible.
   */
  int findRowIndex(int column, int currentRow)const
  {
    int best=-1;
    int bestNumberOfNonZero=0;
    for(int i=currentRow;i<height;i++)
      if(!(*this)[i][column].isZero())
        {
          int nz=0;
          for(int k=column+1;k<width;k++)
            if(!(*this)[i][k].isZero())nz++;
          if(best==-1)
            {
              best=i;
              bestNumberOfNonZero=nz;
            }
          else if(nz<bestNumberOfNonZero)
            {
              best=i;
              bestNumberOfNonZero=nz;
            }
        }
    return best;
  }
  /**
     Performs a Gauss reduction and returns the number of row swaps (and negative scalings)
     done. The result is a matrix in row echelon form. The pivots may
     not be all 1.  In terms of Groebner bases, what is computed is a
     minimal (not necessarily reduced) Groebner basis of the linear
     ideal generated by the rows.  The number of swaps is need if one
     wants to compute the determinant afterwards. In this case it is
     also a good idea to set the flag returnIfZeroDeterminant which
     make the routine terminate before completion if it discovers that
     the determinant is zero.
  */
  int reduce(bool returnIfZeroDeterminant=false, bool integral=false, bool makePivotsOne=false)
  {
    assert(integral || typ::isField());
    assert(!makePivotsOne || !integral);

    int retSwaps=0;
    int currentRow=0;

    for(int i=0;i<width;i++)
      {
        int s=findRowIndex(i,currentRow);

        if(s!=-1)
          {
            if(s!=currentRow)
              {
                swapRows(currentRow,s);
                retSwaps++;
              }
            if(makePivotsOne)
              {//THE PIVOT SHOULD BE SET TO ONE IF INTEGRAL IS FALSE
                if((*this)[currentRow][i].sign()>=0)retSwaps++;
                typ inverse=typ(1)/(*this)[currentRow][i];
                //                if(!rows[currentRow][i].isOne())
                  {
                    for(int k=0;k<width;k++)
                      if(!(*this)[currentRow][k].isZero())
                        (*this)[currentRow][k]*=inverse;
                  }
              }
            for(int j=currentRow+1;j<height;j++)
              if(integral)
                {
                  if(!(*this)[j][i].isZero())
                    {
                      typ s;
                      typ t;

                      typ g=typ::gcd((*this)[currentRow][i],(*this)[j][i],s,t);
                      typ u=-(*this)[j][i]/g;
                      typ v=(*this)[currentRow][i]/g;
                        /* We want the (s,t) vector to be as small as possible.
                         * We are allowed to adjust by multiples of (u,v).
                         * The following computes the correct multiplier (in most cases).
                         */
/*                        {
                          FieldElement multiplier=(s*u+t*v)*((u*u+v*v).inverse());
                          double d=mpq_get_d(*(multiplier.getGmpRationalTemporaryPointer()));
                          multiplier=multiplier.getField()->zHomomorphism(-(((int)(d+10000.5))-10000));
                          s.madd(multiplier,u);
                          t.madd(multiplier,v);
                        }*/
                        for(int k=0;k<width;k++)
                          {
                            typ A=(*this)[currentRow][k];
                            typ B=(*this)[j][k];

                            (*this)[currentRow][k]=s*A+t*B;
                            (*this)[j][k]=u*A+v*B;
                          }
                      }
                  }
                else
                  {
                    if(!(*this)[j][i].isZero())
                      madd(currentRow,-(*this)[j][i]/(*this)[currentRow][i],j);
                  }
              currentRow++;
            }
          else
            if(returnIfZeroDeterminant)return -1;
        }

      return retSwaps;
    }
  /**
     Computes a reduced row echelon form from a row echelon form. In
     Groebner basis terms this is the same as tranforming a minimal
     Groebner basis to a reduced one except that we do not force
     pivots to be 1 unless the scalePivotsToOne parameter is set.
   */
  int REformToRREform(bool scalePivotsToOne=false)
  {
    int ret=0;
    int pivotI=-1;
    int pivotJ=-1;
    while(nextPivot(pivotI,pivotJ))
      {
    	if(scalePivotsToOne)
          (*this)[pivotI]=(*this)[pivotI].toVector()/(*this)[pivotI][pivotJ];
        for(int i=0;i<pivotI;i++)
          if(!(*this)[i][pivotJ].isZero())
            madd(pivotI,-(*this)[i][pivotJ]/(*this)[pivotI][pivotJ],i);
        }
    return ret;
  }
  /**
     This function may be called if the FieldMatrix is in Row Echelon
     Form. The input is a FieldVector which is rewritten modulo the
     rows of the matrix. The result is unique and is the same as the
     normal form of the input vector modulo the Groebner basis of the
     linear ideal generated by the rows of the matrix.
  */
  Vector<typ> canonicalize(Vector<typ> v)const
  {
    assert(typ::isField());
    assert((int)v.size()==getWidth());

    int pivotI=-1;
    int pivotJ=-1;

    while(nextPivot(pivotI,pivotJ))
      if(!v[pivotJ].isZero())
      {
        typ s=-v[pivotJ]/(*this)[pivotI][pivotJ];

        for(int k=0;k<width;k++)
          if(!(*this)[pivotI][k].isZero())
            v[k].madd((*this)[pivotI][k],s);
      }
    return v;
  }
  /**
     Calls reduce() and constructs matrix whose rows forms a basis of
     the kernel of the linear map defined by the original matrix. The
     return value is the new matrix.
   */
  Matrix reduceAndComputeKernel()
  {
    Matrix ret(width-reduceAndComputeRank(),width);

    REformToRREform();

    int k=0;
    int pivotI=-1;
    int pivotJ=-1;
    bool pivotExists=nextPivot(pivotI,pivotJ);
    for(int j=0;j<width;j++)
      {
        if(pivotExists && (pivotJ==j))
          {
            pivotExists=nextPivot(pivotI,pivotJ);
            continue;
          }
        int pivot2I=-1;
        int pivot2J=-1;
        while(nextPivot(pivot2I,pivot2J))
          {
            ret[k][pivot2J]=(*this)[pivot2I][j]/(*this)[pivot2I][pivot2J];
          }
        ret[k][j]=typ(-1);
        k++;
      }
    return ret;
  }
  /**
     Assumes that the matrix has a kernel of dimension 1.
     Calls reduce() and returns a non-zero vector in the kernel.
     If the matrix is an (n-1)x(n) matrix then the returned vector has
     the property that if it was appended as a row to the original matrix
     then the determinant of that matrix would be positive. Of course
     this property, as described here, only makes sense for ordered fields.
     Only allowed for fields at the moment.
   */
  Vector<typ> reduceAndComputeVectorInKernel()
  {
    assert(typ::isField());
    // TODO: (optimization) if the field is ordered, then it is better to just keep track of signs rather than
    // multiplying by sign*diagonalProduct*lastEntry at the end.
    int nswaps=this->reduce();
    typ sign=typ(1-2*(nswaps&1));
    int rank=reduceAndComputeRank();
    assert(rank+1==width);

    REformToRREform();

    Vector<typ> ret(width);

    typ diagonalProduct(1);
    {
      int pivot2I=-1;
      int pivot2J=-1;
      while(nextPivot(pivot2I,pivot2J))
        {
          diagonalProduct*=(*this)[pivot2I][pivot2J];
        }
    }
    {
      int j=nonPivotColumns().front();
      int pivot2I=-1;
      int pivot2J=-1;
      ret[j]=typ(-1);
      // Pretend that we are appending ret to the matrix, and reducing this
      // new row by the previous ones. The last entry of the resulting matrix
      // is lastEntry.
      typ lastEntry=ret[j];
      while(nextPivot(pivot2I,pivot2J))
        {
          ret[pivot2J]=(*this)[pivot2I][j]/(*this)[pivot2I][pivot2J];
          lastEntry-=ret[pivot2J]*ret[pivot2J];
        }
      ret=(sign*(diagonalProduct*lastEntry))*ret;
    }

    return ret;
  }

  /**
     Calls reduce() and returns the rank of the matrix.
   */
  int reduceAndComputeRank()
  {
    reduce(false,!typ::isField(),false);
    int ret=0;
    int pivotI=-1;
    int pivotJ=-1;
    while(nextPivot(pivotI,pivotJ))ret++;
    return ret;
  }
  /**
   * Sort the rows of the matrix.
   */
  struct rowComparer{
    bool operator()(std::pair<Matrix*,int> i, std::pair<Matrix*,int> j) {return ((*i.first)[i.second].toVector()<(*j.first)[j.second].toVector());}
  } theRowComparer;
  void sortRows()
  {
	  std::vector<std::pair<Matrix*,int> > v;
	  for(int i=0;i<height;i++)v.push_back(std::pair<Matrix*,int>(this,i));
	  std::sort(v.begin(),v.end(),theRowComparer);
	  Matrix result(height,width);
	  for(int i=0;i<height;i++)
		  result[i]=(*this)[v[i].second].toVector();
	  data=result.data;
  }
  /**
   * Sort the rows of the matrix and remove duplicate rows.
   */
  void sortAndRemoveDuplicateRows()
  {
    sortRows();
    if(getHeight()==0)return;
    Matrix B(0,getWidth());
    B.appendRow((*this)[0]);
    for(int i=1;i<getHeight();i++)
      if((*this)[i].toVector()!=(*this)[i-1].toVector())B.appendRow((*this)[i].toVector());
    *this=B;
  }
  /**
     Takes two matrices with the same number of columns and construct
     a new matrix which has the rows of the matrix top on the top and
     the rows of the matrix bottom at the bottom. The return value is
     the constructed matrix.
   */
  friend Matrix combineOnTop(Matrix const &top, Matrix const &bottom, MR *mr=get_default_resource())
  {
    assert(bottom.getWidth()==top.getWidth());
    Matrix ret(top.getHeight()+bottom.getHeight(),top.getWidth(),mr);
    for(int i=0;i<top.getHeight();i++)ret[i]=top[i];
    for(int i=0;i<bottom.getHeight();i++)ret[i+top.getHeight()]=bottom[i];

    return ret;
  }
  /**
     Takes two matrices with the same number of rows and construct
     a new matrix which has the columns of the matrix left on the left and
     the columns of the matrix right on the right. The return value is
     the constructed matrix.
   */
  friend Matrix combineLeftRight(Matrix const &left, Matrix const &right, MR *mr=get_default_resource())
  {
    assert(left.getHeight()==right.getHeight());
    Matrix ret(left.getHeight(),left.getWidth()+right.getWidth(),mr);
/*    for(int i=0;i<left.getHeight();i++)
      {
//        for(int j=0;j<left.getWidth();j++)ret.UNCHECKEDACCESS(i,j)=left.UNCHECKEDACCESS(i,j);
//        for(int j=0;j<right.getWidth();j++)ret.UNCHECKEDACCESS(i,j+left.getWidth())=right.UNCHECKEDACCESS(i,j);
        for(int j=0;j<left.getWidth();j++)ret[i][j]=left[i][j];
        for(int j=0;j<right.getWidth();j++)ret[i][j+left.getWidth()]=right[i][j];
      }*/
    ret.setSubMatrix(0,0,left.getHeight(),left.getWidth(),left);
    ret.setSubMatrix(0,left.getWidth(),left.getHeight(),left.getWidth()+right.getWidth(),right);
    return ret;//MOVE??
  }
#if 1
	  /*
	    Scale each row positively to make the row primitive.
	   */
	  void normalizeRows()
	  {
	    for(int i=0;i<getHeight();i++)
	      {
		Vector<typ> a=(*this)[i];
		(*this)[i]=a.normalized();
	      }
	  }
	  /*
	    Return a normalised version of the matrix, where all rows have been made primitive by positive scaling.
	   */
	  Matrix normalizedRows()const
	  {
		  Matrix ret=*this;
		  ret.normalizeRows();
		  return ret;
	  }
#endif
};

	static std::string setColor(int c)
	{
		std::stringstream t;t<<(c&7);
		std::string s(" [");
		s[0]=0x1B;
		return s+"3"+t.str()+"m";
	}
	static std::string setBGColor(int c)
	{
		std::stringstream t;t<<(c&7);
		std::string s(" [");
		s[0]=0x1B;
		return s+"4"+t.str()+"m";
	}
	static std::string resetColors()
	{
		std::string s(" [");
		s[0]=0x1B;
		return s+"0m";
	}
	template <class matrixType> std::string matrixToString(matrixType const &m, matrixType const *colors=0)
	{
		int h=m.getHeight();
		int w=m.getWidth();
		std::vector<std::vector<std::string> > theStrings;
		for(int i=0;i<h;i++)
		{
			std::vector<std::string> row;
			for(int j=0;j<w;j++)row.push_back(m[i][j].toString());
			theStrings.push_back(row);
		}
		std::vector<int> columnWidths(w);
		for(int j=0;j<w;j++)
		{
			int width=1;
			for(int i=0;i<h;i++)if(width<theStrings[i][j].size()+1)width=theStrings[i][j].size();
			columnWidths[j]=width;
		}
		std::stringstream s;
		s<<"h:"<<m.getHeight()<<"w:"<<m.getWidth()<<"\n";
		for(int i=0;i<h;i++)
		{
			for(int j=0;j<w;j++)
			{
				if(colors)
					s<<setColor(((*colors)[i][j].toInt64())&7)<<setBGColor((((*colors)[i][j].toInt64())>>4)&7);
				if(j)s<<" ";
				for(int k=theStrings[i][j].size();k<columnWidths[j];k++)s<<" ";
				s<<theStrings[i][j];
			}
			if(colors)
				s<<resetColors();
			s<<"\n";
		}
		return s.str();
	}

//template<class typ>
//class Matrix<typ,typename std::enable_if<!has_member(typ,POD2),typ>::type >{};

//->decltype(typ.POD)
//template<class typ> class  Matrix<typ,std::enable_if<std::>typ::POD> {};
//template<class typ> class  Matrix<typ,std::enable_if<std::>typ::POD> {};

//template<class typ>
//class Matrix<typ,typ>{};

//template<class typ>
//class Matrix<typ,typename typ::POD>;

// If typ has POD defined we use it. Otherwise, we just use typ itself.

//template<typename typ>
//class Matrix<typ,std::conditional_t<has_member(typ,POD2),typename typ::POD, typ>>{};
//	hasPod((typ*)(0))
//template<typename typ>
//class Matrix<typ,std::conditional_t<
//				hasPod((typ*)(0)),
//				typename typ::POD,
//				typ>>{};

//template<class typ>
//class Matrix<typ,typename std::enable_if<has_member(typ,POD2),typename typ::POD> >{};
//template<class typ>
//class Matrix<typ,typ>;//{} Matrix;

typedef Matrix<Integer> ZMatrix;
typedef Matrix<Rational> QMatrix;
typedef Matrix<int> IntMatrix;

inline QMatrix ZToQMatrix(ZMatrix const &m)
{
  QMatrix ret(m.getHeight(),m.getWidth());
  for(int i=0;i<m.getHeight();i++)ret[i]=ZToQVector(m[i]);
  return ret;
}

inline ZMatrix QToZMatrixPrimitive(QMatrix const &m)
{
  ZMatrix ret(m.getHeight(),m.getWidth());
  for(int i=0;i<m.getHeight();i++)ret[i]=QToZVectorPrimitive(m[i]);
  return ret;
}


inline IntMatrix ZToIntMatrix(ZMatrix const &m)
{
  IntMatrix ret(m.getHeight(),m.getWidth());
  for(int i=0;i<m.getHeight();i++)ret[i]=ZToIntVector(m[i]);
  return ret;
}


inline ZMatrix IntToZMatrix(IntMatrix const &m)
{
  ZMatrix ret(m.getHeight(),m.getWidth());
  for(int i=0;i<m.getHeight();i++)ret[i]=IntToZVector(m[i]);
  return ret;
}

inline QMatrix canonicalizeSubspace(QMatrix const &m)
{
  QMatrix temp=m;
  temp.reduce();
  temp.REformToRREform();
  temp.removeZeroRows();
  return temp;
}

inline ZMatrix canonicalizeSubspace(ZMatrix const &m)
{
  return QToZMatrixPrimitive(canonicalizeSubspace(ZToQMatrix(m)));
}


inline QMatrix kernel(QMatrix const &m)
{
  QMatrix temp=m;
  return temp.reduceAndComputeKernel();
}

inline ZMatrix kernel(ZMatrix const &m)
{
  return QToZMatrixPrimitive(kernel(ZToQMatrix(m)));
}

static_assert(std::is_move_constructible<Matrix<int>>::value,"Not move constructible!" );
static_assert(std::is_move_assignable<Matrix<int>>::value,"Not move assignable!" );
static_assert(std::is_nothrow_move_assignable<Matrix<int>>::value,"Not nothrow move assignable!" );
static_assert(std::is_nothrow_move_constructible<Matrix<int>>::value,"Not nothrow move constructible!" );
static_assert(std::is_nothrow_destructible<Matrix<int>>::value,"Not nothrow destructible!" );
}


#endif /* LIB_ZMATRIX_H_ */
