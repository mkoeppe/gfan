/*
 * lib_zvector.h
 *
 *  Created on: Sep 6, 2010
 *      Author: anders
 */

#ifndef GFANLIB_ZVECTOR_H_INCLUDED
#define GFANLIB_ZVECTOR_H_INCLUDED

#include <vector>
#include <list>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "gfanlib_memoryresource.h"

#include "gfanlib_z.h"
#include "gfanlib_q.h"

namespace gfan{

static std::vector<std::string> parseSequence(std::istream &s)
{
	s>>std::noskipws;
	std::vector<std::string> ret;
	std::vector<char> stack;
	int state=0;
	std::string currentString;
	while(!s.eof())
	{
		char c;s>>c;
//				std::cerr<<"PARSED::"<<c<<"\n";

		if(c=='(' || c=='{')
		{
			stack.push_back(c);
			if(stack.size()>1)currentString.push_back(c);
		}
		else
		if(c==',' && stack.size()==1)
		{
		  //			ret.push_back(currentString.str());
		  //			currentString=std::stringstream("");
		  std::string t;
		  std::swap(t,currentString);
		  ret.emplace_back(std::move(t));
		}
		else
		if(c==')' || c=='}')
		{
		  if(stack.size()>1)currentString.push_back(c);
//			assert(stack.size());
			stack.pop_back();
			if(stack.size()==0)
			{
			  ret.push_back(std::move(currentString));
				break;
			}
		}
		else
		  currentString.push_back(c);
	}
	return ret;
}


inline void outOfRange(int i, int n)
{
  std::cerr<<"Index out of range. i="<<i<<" n="<<n<<std::endl;
  assert(0);
}

template <class typ> class Vector{
private:
  pmrvector<typ> v;
public:
  //--------------
  // Constructors
  //--------------
  Vector(Vector &&a)=default;
  Vector(Vector &&a, MR *mr):v(std::move(v),mr){};
  Vector &operator=(Vector &&a)noexcept{v=std::move(a.v);return *this;};
  Vector &operator=(Vector const &a)=default;
  Vector(const Vector &a, MR *mr=get_default_resource()):v(a.v,mr){
  }
  Vector(int n, MR *mr=get_default_resource()):v(n,typ{}/*typ()*/,mr){
    assert(n>=0);
  };
/*  Vektor(const typ *p, int n):v(n){
    assert(n>=0);
    support=0;for(int i=0;i<n;i++)v[i]=p[i];
  };*/
  ~Vector(){
  };
  Vector(){
  };
  // Template for conversion
  template<typename> friend class Vector;
  template <typename otherTyp>
  explicit Vector(Vector<otherTyp> const &a, MR *mr=get_default_resource()):
  	  Vector(a.v.size(),mr)
  {
	  for(int i=0;i<a.v.size();i++)
		  (*this)[i]=typ(a[i]);
  }


  static Vector standardVector(int n, int i, MR *mr=get_default_resource())
    {
      Vector v(n,mr);
      v[i]=typ(1);
      return v;
    }
  static Vector allOnes(int n)
    {
      Vector v(n);
      for(int i=0;i<n;i++)
        v[i]=typ(1);
      return v;
    }
  //-------------
  // Conversions
  //-------------
/*  template<class T> Vektor(const Vektor<T>& c)
          :v(c.size())
        {
          for(int i=0;i<size();i++)v[i]=typ(c[i]);
        }
*/
  //--------
  // Access
  //--------
  typ& operator[](int n)__attribute__((always_inline))
    {
      if(!(n>=0 && n<(int)v.size()))outOfRange(n,v.size());
      return (v[n]);
    }
  const typ& operator[](int n)const __attribute__((always_inline)){assert(n>=0 && n<(int)v.size());return (v[n]);}
  typ& UNCHECKEDACCESS(int n)__attribute__((always_inline)){return (v[n]);}
  const typ& UNCHECKEDACCESS(int n)const __attribute__((always_inline)){return (v[n]);}

  //-------------
  // STL related
  //-------------
  unsigned int size()const{return v.size();};
  void resize(int n){v.resize(n,typ());};
  void grow(int i){if(size()<i)resize(i);}
  void push_back(typ a)
  {
    v.push_back(a);
  }
  void sort()
    {
      std::sort(v.begin(),v.end());
    }
  bool nextPermutation()
    {
      return std::next_permutation(v.begin(),v.end());
    }
  bool operator<(const Vector & b)const
    {
      if(size()<b.size())return true;
      if(size()>b.size())return false;
      for(unsigned i=0;i<size();i++)
        {
          if(v[i]<b[i])return true;
          if(b[i]<v[i])return false;
        }
      return false;
    }

  //-----------------
  // Arithmetic fast
  //-----------------
  typ sum()const{typ f;for(typename pmrvector<typ>::const_iterator i=v.begin();i!=v.end();i++)f+=*i;return f;};
  Vector& operator+=(const Vector& q){assert(size()==q.size());typename pmrvector<typ>::const_iterator j=q.v.begin();for(typename pmrvector<typ>::iterator i=v.begin();i!=v.end();++i,++j)(*i)+=(*j);return *this;}
  Vector& operator-=(const Vector& q){assert(size()==q.size());typename pmrvector<typ>::const_iterator j=q.v.begin();for(typename pmrvector<typ>::iterator i=v.begin();i!=v.end();++i,++j)(*i)-=(*j);return *this;}
  Vector& operator/=(const Vector& q){assert(size()==q.size());typename pmrvector<typ>::const_iterator j=q.v.begin();for(typename pmrvector<typ>::iterator i=v.begin();i!=v.end();++i,++j)(*i)/=(*j);return *this;}
  inline friend typ dot(const Vector& p, const Vector& q){assert(p.size()==q.size());typ s;typename pmrvector<typ>::const_iterator j=q.v.begin();for(typename pmrvector<typ>::const_iterator i=p.v.begin();i!=p.v.end();++i,++j)s+=(*i)*(*j);return s;}
//  inline friend int64 dotLong(const Vektor& p, const Vektor& q){assert(p.size()==q.size());int64 s=0;for(int i=0;i<p.size();i++)s+=(int64)p[i]*(int64)q[i];return s;}
  bool operator==(const Vector & q)const{if(size()!=q.size())return false;typename pmrvector<typ>::const_iterator j=q.v.begin();for(typename pmrvector<typ>::const_iterator i=v.begin();i!=v.end();++i,++j)if(!((*i)==(*j)))return false;return true;}
  bool operator!=(const Vector & q)const {return !(operator==(q));}
  bool isZero() const
    {
      for(typename pmrvector<typ>::const_iterator i=v.begin();i!=v.end();++i)if(!(*i).isZero())return false;
      return true;
    }
  bool isPositive() const
    {
      for(typename pmrvector<typ>::const_iterator i=v.begin();i!=v.end();++i)if((*i).sign()<=0)return false;
      return true;
    }
  bool isNonNegative() const
    {
      for(typename pmrvector<typ>::const_iterator i=v.begin();i!=v.end();++i)if(i->sign()<0)return false;
      return true;
    }
/*  int max()const
  {
    int ret=-0x7fffffff; //not completely correct, but kind of works for 64bit
    for(int i=0;i<v.size();i++)if(ret<v[i])ret=v[i];
    return ret;
  }
  int min()const
  {
    int ret=0x7fffffff;
    for(int i=0;i<v.size();i++)if(ret>v[i])ret=v[i];
    return ret;
  }
*/
  friend bool dependent(const Vector& p, const Vector& q)
    {
          /*
      typ pp=dot(p,p);
      typ qq=dot(q,q);
      typ pq=dot(p,q);
      return pq*pq==pp*qq;
*/
          unsigned n=p.size();
          assert(n==q.size());
          unsigned i;
          for(i=0;i<n;i++)
          {
            if(!p.v[i].isZero())break;
          }
          if(i==n)return true;
          if(q.v[i].isZero())return q.isZero();
          typ a=p.v[i];
          typ b=q.v[i];
          for(unsigned j=0;j<n;j++)
            if(a*q.v[j]!=b*p.v[j])return false;
          return true;
    }

  //-----------------
  // Arithmetic slow
  //-----------------
  inline friend Vector operator*(typ s, const Vector& q){Vector p=q;for(unsigned i=0;i<q.size();i++)p[i]*=s;return p;}
//  inline friend Vektor operator/(const Vektor& q, typ s){Vektor p=q;for(int i=0;i<q.size();i++)p[i]/=s;return p;}
/*  inline friend Vector operator*(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1.v[i]*=q.v[i];return p1;}
  inline friend Vektor operator+(const Vektor& p, const Vektor& q){assert(p.size()==q.size());Vektor p1=p;for(int i=0;i<p.size();i++)p1[i]+=q[i];return p1;}
*/
  inline friend Vector operator/(const Vector& p, typ const &s){Vector ret(p.size());for(unsigned i=0;i<p.size();i++)ret[i]=p[i]/s;return ret;}
  inline friend Vector operator+(const Vector& p, const Vector& q){assert(p.size()==q.size());Vector p1=p;for(unsigned i=0;i<p.size();i++)p1[i]+=q[i];return p1;}
  inline friend Vector operator-(const Vector& p, const Vector& q){assert(p.size()==q.size());Vector p1=p;for(unsigned i=0;i<p.size();i++)p1[i]-=q[i];return p1;}
  inline friend Vector operatorMinus(const Vector& p, const Vector& q, MR *mr=get_default_resource()){assert(p.size()==q.size());Vector p1(p,mr);for(unsigned i=0;i<p.size();i++)p1[i]-=q[i];return p1;}
  inline friend Vector max(const Vector& p, const Vector& q){assert(p.size()==q.size());Vector p1=p;for(unsigned i=0;i<p.size();i++)if(p1[i]<q[i])p1[i]=q[i];return p1;}
  inline friend Vector min(const Vector& p, const Vector& q){assert(p.size()==q.size());Vector p1=p;for(unsigned i=0;i<p.size();i++)if(p1[i]>q[i])p1[i]=q[i];return p1;}

  friend Vector operator-(const Vector &b)
  {
    Vector ret(b.size());
    for(unsigned i=0;i<b.size();i++)ret[i]=-b[i];
    return ret;
  }
  //------------------
  // Monomial related
  //------------------
/*  int divides(const Vektor& q) const
    {
      assert(size()==q.size());
      int n=v.size();
      for(int i=0;i<n;i++)
        {
          if(v[i]>0)if(q.v[i]<v[i])return 0;
        }
      return 1;
    }
  inline friend bool relativelyPrime(const Vektor& p, const Vektor& q)
    {
      assert(p.size()==q.size());
      int n=p.size();
      for(int t=0;t<n;t++)if((p[t]>0)&&(q[t]>0)) return false;
      return true;
    }
  Vektor supportVector()const
    {
      Vektor r(v.size());
      for(int i=0;i<size();i++)
        r[i]=(v[i]!=0);
      return r;
    }
*/
  //------------------------------
  // Subvectors and concatenation
  //------------------------------
  Vector subvector(int begin, int end, MR *mr=get_default_resource())const
    {
      assert(begin>=0);
      assert(end<=(int)size());
      assert(end>=begin);
      Vector ret(end-begin,mr);
      for(int i=0;i<end-begin;i++)
        ret[i]=v[begin+i];
      return ret;
    }
/*  Vektor subvector(list<int> const &chosenIndices)const
  {
    Vektor ret(chosenIndices.size());
    int I=0;
    for(list<int>::const_iterator i=chosenIndices.begin();i!=chosenIndices.end();i++,I++)ret[I]=v[*i];
    return ret;
  }
*/
  friend Vector concatenation(Vector const &a, Vector const &b, MR *mr=get_default_resource())
  {
    Vector ret(a.size()+b.size(),mr);
    for(int i=0;i<a.size();i++)ret[i]=a[i];
    for(int i=0;i<b.size();i++)ret[i+a.size()]=b[i];
    return ret;
  }
  //----------------------
  // Zero/nonZero entries
  //----------------------
/*  int indexOfLargestNonzeroEntry()const
  {
    int ret=-1;
    for(int i=0;i<v.size();i++)
      {
        if(v[i])ret=i;
      }
    return ret;
  }
  Vektor supportIndices()const
  {
    Vektor ret(0);
    for(int i=0;i<v.size();i++)
      if(v[i]!=0)ret.push_back(i);
    return ret;
  }
  Vektor supportAsZeroOneVector()const
  {
    Vektor ret(v.size());
    for(int i=0;i<v.size();i++)ret[i]=bool(v[i]);
    return ret;
  }
  void calcsupport(void)
    {
      support=0;
      for(int i=0;i<v.size();i++)support=(support<<1)|(((v[i]>0)==true)&1);
    }
*/
  friend std::ostream &operator<<(std::ostream &f, Vector const &a)
  {
    f<<"(";
    for(typename pmrvector<typ>::const_iterator i=a.v.begin();i!=a.v.end();i++)
    {
      if(i!=a.v.begin()) f<<",";
      f<<*i;
    }
    return f<<")";
  }

  friend std::istream &operator>>(std::istream &f, Vector &a)
  {
	  auto v=parseSequence(f);
	  a.resize(v.size());
	  int i=0;
	  for(auto c:v)
		  a[i++]=typ(c);
	  return f;
  }

  std::string toString()const
  {
	  std::stringstream f;
	  f<<*this;
	  return f.str();
  }

  typ gcd()const
  {
//    typ temp1,temp2;
    typ ret(int(0));
    //    std::cerr<<"gcd1\n";
    for(unsigned i=0;i<size();i++)
      ret=typ::gcd(ret,v[i]/*,temp1,temp2*/);
    //std::cerr<<"gcd2\n";
    //    std::cerr<<"gcd"<<ret<<"\n";
    return ret;
  }
  Vector normalized()const
  {
    //    std::cerr<<"Before Normalizing"<<toString()<<"\n";
      //    std::cerr<<"norm1\n";
    assert(!typ::isField());
    //std::cerr<<"norm2\n";
    auto g=gcd();
    if(!g.isZero())
      return (*this)/gcd();
    return *this;
  }
};

typedef Vector<Integer> ZVector;
typedef Vector<Rational> QVector;
typedef Vector<int> IntVector;

inline QVector ZToQVector(ZVector const &v)
{
  QVector ret(v.size());
  for(unsigned i=0;i<v.size();i++)ret[i]=Rational(v[i]);
  return ret;
}


inline IntVector ZToIntVector(ZVector const &v)
{
  IntVector ret(v.size());
  for(unsigned i=0;i<v.size();i++)ret[i]=v[i].toInt();
  return ret;
}


inline ZVector IntToZVector(IntVector const &v)
{
  ZVector ret(v.size());
  for(unsigned i=0;i<v.size();i++)ret[i]=Integer(v[i]);
  return ret;
}

inline ZVector QToZVectorPrimitive(QVector const &v)
{
    int n=v.size();
    ZVector ret(n);

    mpz_t lcm;
    mpz_t gcd;
    mpz_init_set_ui(lcm, 1);
    mpz_init_set_ui(gcd, 0);

    mpq_t a;
    mpq_init(a);
    for(int j=0;j<n;j++)
      {
        v[j].setGmp(a);
        if(mpz_cmp_si(mpq_denref(a),1)!=0)
          mpz_lcm(lcm,lcm,mpq_denref(a));
        if(mpz_sgn(mpq_numref(a))!=0)
          mpz_gcd(gcd,gcd,mpq_numref(a));
      }
    mpq_clear(a);
    if(mpz_sgn(gcd)!=0)//v is non-zero
      {
        if((mpz_cmp_si(lcm,1)==0)&&(mpz_cmp_si(gcd,1)==0)) //gcd=lcm=1
          {
            mpq_t a;
            mpq_init(a);

            for(int i=0;i<n;i++)
              {
                v[i].setGmp(a);
                ret[i]=Integer(mpq_numref(a));
              }
            mpq_clear(a);
          }
        else
          {
            mpq_t a;
            mpq_init(a);
            mpz_t tempA;
            mpz_t tempB;
            mpz_init(tempA);
            mpz_init(tempB);
            for(int i=0;i<n;i++)
              {
                v[i].setGmp(a);
                mpz_set(tempA,mpq_denref(a));
                mpz_set(tempB,mpq_numref(a));
                mpz_mul(tempA,gcd,tempA);
                mpz_mul(tempB,lcm,tempB);
                mpz_divexact(tempA,tempB,tempA);
                ret[i]=Integer(tempA);
              }
            mpz_clear(tempB);
            mpz_clear(tempA);
            mpq_clear(a);
          }
      }
    mpz_clear(gcd);
    mpz_clear(lcm);

    return ret;
}

static_assert(std::is_move_constructible<Vector<int>>::value,"Not move constructible!" );
static_assert(std::is_move_assignable<Vector<int>>::value,"Not move assignable!" );
static_assert(std::is_nothrow_move_assignable<Vector<int>>::value,"Not nothrow move assignable!" );
static_assert(std::is_nothrow_move_constructible<Vector<int>>::value,"Not nothrow move constructible!" );
static_assert(std::is_nothrow_destructible<Vector<int>>::value,"Not nothrow destructible!" );

}

#endif /* LIB_ZVECTOR_H_ */
