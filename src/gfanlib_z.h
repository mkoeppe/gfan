/*
 * lib_z.h
 *
 *  Created on: Sep 28, 2010
 *      Author: anders
 */

#ifndef LIB_Z_H_
#define LIB_Z_H_

#include <string.h>
#include <ostream>
#include <iostream>
#define OLD 1
#if OLD
#include "gmp.h"

namespace gfan{

  class IntGMP
  {
    mpz_t value;
  public:
    typedef IntGMP POD;//Not actually a POD type. But here we can specify a POD type if one exists.
    static bool isField()
    {
      return false;
    }
    IntGMP()
    {
      mpz_init(value);
    }
    IntGMP(signed long int value_)
    {
      mpz_init(value);
      mpz_set_si(value,value_);
    }
    IntGMP(IntGMP const & value_)
    {
      mpz_init_set(value,value_.value);
    }
    IntGMP(mpz_t const value_)
    {
      mpz_init_set(value,value_);
    }
    IntGMP(std::string const&s) {
      mpz_init(value);
      mpz_set_str(value, s.c_str(), 10);
    }

    ~IntGMP()
    {
      mpz_clear(value);
    }
    IntGMP& operator=(const IntGMP& a)
      {
        const IntGMP *A=(const IntGMP*)&a;
        if (this != A) {
          mpz_clear(value);
          mpz_init_set(value, a.value);
        }
        return *this;
      }
    bool isZero()const{
      return mpz_sgn(value)==0;
    }
    friend std::ostream &operator<<(std::ostream &f, IntGMP const &a)
    {
      void (*freefunc)(void *, size_t);
      mp_get_memory_functions(0,0,&freefunc);
      char *str=mpz_get_str(0,10,a.value);
      f<<str;
      freefunc(str,strlen(str)+1);
      return f;
    }
    IntGMP& operator+=(const IntGMP& a)
      {
        mpz_add(value,value,a.value);
        return *this;
      }
    IntGMP& operator-=(const IntGMP& a)
      {
        mpz_sub(value,value,a.value);
        return *this;
      }
    IntGMP& operator*=(const IntGMP& a)
      {
        mpz_mul(value,value,a.value);
        return *this;
      }
    IntGMP& operator/=(const IntGMP& a)
      {
        mpz_div(value,value,a.value);
        return *this;
      }
    friend IntGMP operator-(const IntGMP &b)
    {
      IntGMP ret;
      ret-=b;
      return ret;
    }
    IntGMP operator+(const IntGMP &a)const
    {
      IntGMP ret(*this);
      ret+=a;
      return ret;
    }
    IntGMP operator-(const IntGMP &a)const
    {
      IntGMP ret(*this);
      ret-=a;
      return ret;
    }
    IntGMP operator*(const IntGMP &a)const
    {
      IntGMP ret(*this);
      ret*=a;
      return ret;
    }
    IntGMP operator/(const IntGMP &a)const
    {
      IntGMP ret(*this);
      ret/=a;
      return ret;
    }
    void madd(const IntGMP &a,const IntGMP &b)
      {
        mpz_t temp;
        mpz_init(temp);
        mpz_mul(temp,a.value,b.value);
        mpz_add(value,value,temp);
        mpz_clear(temp);
      }
    bool operator<(const IntGMP &a)const
    {
      return mpz_cmp(value,a.value)<0;
    }
    bool operator==(const IntGMP &a)const
    {
      return mpz_cmp(value,a.value)==0;
    }
    bool operator!=(const IntGMP &a)const
    {
      return mpz_cmp(value,a.value)!=0;
    }
    int sign()const
    {
      return mpz_sgn(value);
    }
    static IntGMP gcd(IntGMP const &a, IntGMP const &b, IntGMP &s, IntGMP &t)
    {
      mpz_t r;
      mpz_init(r);
      mpz_gcdext(r,s.value,t.value,a.value,b.value);
      IntGMP ret(r);
      mpz_clear(r);
      return ret;
    }
    /**
     * Assigns the value to z. z must have been initialized as a gmp variable.
     */
    void setGmp(mpz_t z)const
    {
      mpz_set(z,value);
    }
    /**
     * Returns a value which is useful for computing hash functions.
     */
    signed long int hashValue()const
    {
      return mpz_get_si(value);
    }
    bool fitsInInt()const
    {
      mpz_t v;
      mpz_init(v);
      this->setGmp(v);
      bool ret=(mpz_fits_sint_p(v)!=0);
      mpz_clear(v);
      return ret;
    }
    int64_t toInt()const//works for 64 bit integers
    {
      mpz_t v;
      mpz_init(v);
      this->setGmp(v);
      int ret=0;
  //    if(mpz_fits_sint_p(v))
        ret=mpz_get_si(v);
  //    else
  //      ok=false;
      mpz_clear(v);
      return ret;
    }
  };

/*
 * Integer is supposed to be a more efficient version of IntGMP. However, the current implementation is not portable.
 * Therefore the following if can be used to enable/disable the efficient version.
 */
#if 0 // enable more efficient integer here
  typedef IntGMP Integer;
#else
  class Integer2
{
  /* Integer2 represents an integer by either an GMP integer or an immediate value.
   *
   * The value needs to be interpreted in two different ways depending on whether the pointer value._mp_d is odd or even.
   * We assume that odd pointers are not actually used for storing limbs - although that is not guaranteed on all architectures.
   * If the pointer is even, then Integer uses gmp to store its value in an mpz_t. If the pointer is odd, then a signed 32 bit
   * value is stored in value._mp_alloc.
   *
   * Because the above conventions may change in the future, access to the machine integer representation has been encapsulated
   * in the private methods below.
   */
  mpz_t value;
  bool hasLimbs()const
  {
	  return !(((long int)(value->_mp_d))&1);
  }
  void initSetInt32(std::int32_t v)//should only be used in constructor or after mpz_clear has been called
  {
	  value->_mp_alloc=v;
	  *((long int*)&(value->_mp_d))=1;
  }
  void setInt32(std::int32_t v)
  {
	  assert(!hasLimbs());
	  value->_mp_alloc=v;
  }
  std::int32_t getInt32()const
  {
	  assert(!hasLimbs());
	  return int32_t(value[0]._mp_alloc);
  }
  void extend()
  {
	  assert(!hasLimbs());
	  mpz_init_set_si(value,getInt32());
  }
  bool fitsInInt32(signed long int v)
  {
	  return (v>=-65536*(long)0x00008000)&&(v<=0x7fffffff);
  }
  void tryShrinking()
  {
	  if(hasLimbs())
		  if(mpz_fits_slong_p(value))
			  if(fitsInInt32(mpz_get_si(value)))
			  {
				  auto v=mpz_get_si(value);
				  mpz_clear(value);
				  initSetInt32(v);
			  }
  }
  void dump()const
  {
	  std::cerr<<value->_mp_d<<" "<<value->_mp_alloc<<"\n";
  }
public:
  typedef Integer2 POD;//Not actually a POD type. But here we can specify a POD type if one exists.
  static bool isField()
  {
    return false;
  }
  Integer2()
  {
	  initSetInt32(0);
//    mpz_init(value);
  }
  Integer2(signed long int value_)
  {
	  if(fitsInInt32(value_))
		  initSetInt32(value_);
	  else
		  mpz_init_set_si(value,value_);
  }
  Integer2(Integer2 const & value_)
  {
	  if(value_.hasLimbs())
		  mpz_init_set(value,value_.value);
	  else
		  value[0]=value_.value[0];
  }
  Integer2(mpz_t const value_)
  {
	  if(mpz_fits_sint_p(value_))
		  initSetInt32(mpz_get_si(value_));
	  else
		  mpz_init_set(value,value_);
  }
  Integer2(std::string const&s) {
    mpz_init(value);
    mpz_set_str(value, s.c_str(), 10);
    tryShrinking();
  }

  ~Integer2()
  {
	  if(hasLimbs())
		  mpz_clear(value);
  }
  Integer2& operator=(const Integer2& a)
    {
      const Integer2 *A=(const Integer2*)&a;
      if (this != A) {
    	  if(hasLimbs())
    		  if(a.hasLimbs())
				  mpz_set(value,a.value);
    		  else
    		  {
    			  mpz_clear(value);
    			  initSetInt32(a.getInt32());
    		  }
    	  else
    		  if(a.hasLimbs())
    			  mpz_init_set(value, a.value);
    		  else
    			  setInt32(a.getInt32());
      }
      return *this;
    }
  bool isZero()const{
//	  dump();
	  if(hasLimbs())
		  return mpz_sgn(value)==0;
	  return getInt32()==0;
  }
  friend std::ostream &operator<<(std::ostream &f, Integer2 const &a)
  {
	  if(a.hasLimbs())
	  {
		void (*freefunc)(void *, size_t);
		mp_get_memory_functions(0,0,&freefunc);
		char *str=mpz_get_str(0,10,a.value);
		f<<str;
		freefunc(str,strlen(str)+1);
	  }
	  else
		  f<<a.getInt32();
    return f;
  }
  Integer2& operator+=(const Integer2& a)
    {
	  if(hasLimbs())
		  if(a.hasLimbs())
			  mpz_add(value,value,a.value);
		  else
		  {
			  auto v=a.getInt32();
			  if(v>0)
				  mpz_add_ui(value,value,v);
			  else
				  mpz_sub_ui(value,value,((unsigned long int)0xffffffff)&((unsigned long int)-v));
		  }
	  else
		  if(a.hasLimbs())
		  {
			  mpz_init_set_si(value,getInt32());
			  mpz_add(value,value,a.value);
		  }
		  else
		  {
			  auto r=((long int)getInt32())+(long int)a.getInt32();
			  if(fitsInInt32(r))
				  setInt32(r);
			  else
				  mpz_init_set_si(value,r);
		  }
	  tryShrinking();
      return *this;
    }
  Integer2& operator-=(const Integer2& a)
    {
	  if(hasLimbs())
		  if(a.hasLimbs())
			  mpz_sub(value,value,a.value);
		  else
		  {
			  auto v=a.getInt32();
			  if(v>0)
				  mpz_sub_ui(value,value,v);
			  else
				  mpz_add_ui(value,value,((unsigned long int)0xffffffffL)&((unsigned long int)-v));
		  }
	  else
		  if(a.hasLimbs())
		  {
			  mpz_init_set_si(value,getInt32());
			  mpz_sub(value,value,a.value);
		  }
		  else
		  {
			  auto r=((long int)getInt32())-(long int)a.getInt32();
			  if(fitsInInt32(r))
				  setInt32(r);
			  else
				  mpz_init_set_si(value,r);
		  }
	  tryShrinking();
      return *this;
    }
  Integer2& operator*=(const Integer2& a)
    {
	  if(hasLimbs())
		  if(a.hasLimbs())
			  mpz_mul(value,value,a.value);
		  else
			  mpz_mul_si(value,value,a.getInt32());
	  else
		  if(a.hasLimbs())
		  {
			  mpz_init_set_si(value,getInt32());
			  mpz_mul(value,value,a.value);
		  }
		  else
		  {
			  auto r=((long int)getInt32())*(long int)a.getInt32();
			  if(fitsInInt32(r))
				  setInt32(r);
			  else
				  mpz_init_set_si(value,r);
		  }
	  tryShrinking();
	  return *this;
    }
  Integer2& operator/=(const Integer2& a)
	{
	  if(hasLimbs())
		  if(a.hasLimbs())
			  mpz_fdiv_q(value,value,a.value);
		  else
		  {
			  int v=a.getInt32();
			  if(v>0)
				  mpz_fdiv_q_ui(value,value,a.getInt32());
			  else
			  {
				  mpz_fdiv_q_ui(value,value,((unsigned long int)0xffffffffL)&((unsigned long int)-v));
				  mpz_neg(value,value);
			  }
		  }
	  else
		  if(a.hasLimbs())
		  {
			  mpz_init_set_si(value,getInt32());
			  mpz_fdiv_q(value,value,a.value);
		  }
		  else
			  setInt32(getInt32()/a.getInt32());
	  tryShrinking();
      return *this;
    }
  Integer2& exactlyDivideBy(const Integer2 &a)
  {
	  assert(0);//should be optimised
	  if(hasLimbs())
		  if(a.hasLimbs())
			  mpz_divexact(value,value,a.value);
		  else
		  {
			  int v=a.getInt32();
			  if(v>0)
				  mpz_divexact_ui(value,value,a.getInt32());
			  else
			  {
				  mpz_divexact_ui(value,value,((unsigned long int)0xffffffffL)&((unsigned long int)-v));
				  mpz_neg(value,value);
			  }
		  }
	  else
		  if(a.hasLimbs())
		  {
			  mpz_init_set_si(value,getInt32());
			  mpz_divexact(value,value,a.value);
		  }
		  else
			  setInt32(getInt32()/a.getInt32());
	  tryShrinking();
      return *this;
  }
  friend Integer2 operator-(const Integer2 &b)
  {
    Integer2 ret;
    ret-=b;
    return ret;
  }
  Integer2 operator+(const Integer2 &a)const
  {
    Integer2 ret(*this);
    ret+=a;
    return ret;
  }
  Integer2 operator-(const Integer2 &a)const
  {
    Integer2 ret(*this);
    ret-=a;
    return ret;
  }
  Integer2 operator*(const Integer2 &a)const
  {
    Integer2 ret(*this);
    ret*=a;
    return ret;
  }
  Integer2 operator/(const Integer2 &a)const
  {
    Integer2 ret(*this);
    ret/=a;
    return ret;
  }
  Integer2 exactQuotientBy(const Integer2 &a)const
  {
    Integer2 ret(*this);
    ret.exactlyDivideBy(a);
    return ret;
  }
  void madd(const Integer2 &a,const Integer2 &b)
    {
/*      mpz_t temp;
      mpz_init(temp);
      mpz_mul(temp,a.value,b.value);
      mpz_add(value,value,temp);
      mpz_clear(temp);*/

	  if((!hasLimbs())&&(a.hasLimbs()||b.hasLimbs()))
		  mpz_init_set_si(value,getInt32());

	  if(a.hasLimbs())
		  if(b.hasLimbs())
			  mpz_addmul(value,a.value,b.value);
		  else
			  if(b.getInt32()>0)
				  mpz_addmul_ui(value,a.value,b.getInt32());
			  else
				  mpz_submul_ui(value,a.value,((unsigned long int)0xffffffffL)&((unsigned long int)-b.getInt32()));
	  else
		  if(b.hasLimbs())
			  if(a.getInt32()>0)
				  mpz_addmul_ui(value,b.value,a.getInt32());
			  else
				  mpz_submul_ui(value,b.value,((unsigned long int)0xffffffffL)&((unsigned long int)-a.getInt32()));
		  else
		  {
			  if(hasLimbs())
			  {
				  auto prod=((long int)a.getInt32())*(long int)b.getInt32();
				  if(prod>=0)
					  mpz_add_ui(value,value,prod);
				  else
					  mpz_sub_ui(value,value,-prod);
			  }
			  else
			  {
				  auto r=getInt32()+((long int)a.getInt32())*(long int)b.getInt32();
				  if(fitsInInt32(r))
					  setInt32(r);
				  else
					  mpz_init_set_si(value,r);
			  }
		  }
	  tryShrinking();
    }
  bool operator<(const Integer2 &a)const
  {
	  if(hasLimbs())
		  if(a.hasLimbs())
			  return mpz_cmp(value,a.value)<0;
		  else
			  return mpz_cmp_si(value,a.getInt32())<0;
	  else
		  if(a.hasLimbs())
			  return -mpz_cmp_si(a.value,getInt32())<0;
		  else
			  return getInt32()<a.getInt32();
  }
  bool operator==(const Integer2 &a)const
  {
	  if(hasLimbs())
		  if(a.hasLimbs())
			  return mpz_cmp(value,a.value)==0;
		  else
			  return mpz_cmp_si(value,a.getInt32())==0;
	  else
		  if(a.hasLimbs())
			  return mpz_cmp_si(a.value,getInt32())==0;
		  else
			  return getInt32()==a.getInt32();
  }
  bool operator!=(const Integer2 &a)const
  {
	  return !(*this==a);
  }
  int sign()const
  {
	  if(hasLimbs())
		  return mpz_sgn(value);
	  else
		  if(getInt32()<0)return -1;
		  else
			  if(getInt32()>0)return 1;
			  else return 0;
  }
  static Integer2 gcd(Integer2 const &a, Integer2 const &b, Integer2 &s, Integer2 &t)
  {
	  auto A(a);
	  auto B(b);
	  if(!A.hasLimbs())A.extend();
	  if(!B.hasLimbs())B.extend();
	  if(!s.hasLimbs())s.extend();
	  if(!t.hasLimbs())t.extend();

    mpz_t r;
    mpz_init(r);
    mpz_gcdext(r,s.value,t.value,A.value,B.value);
    Integer2 ret(r);
    mpz_clear(r);

    	s.tryShrinking();
    	t.tryShrinking();

    return ret;
  }
  static Integer2 gcd(Integer2 const &a, Integer2 const &b)
  {
	  auto A(a);
	  auto B(b);
	  if(!A.hasLimbs())A.extend();
	  if(!B.hasLimbs())B.extend();

    mpz_t r;
    mpz_init(r);
    mpz_gcd(r,A.value,B.value);
    Integer2 ret(r);
    mpz_clear(r);

    return ret;
  }
  /**
   * Assigns the value to z. z must have been initialized as a gmp variable.
   */
  void setGmp(mpz_t z)const
  {
	  if(hasLimbs())
		  mpz_set(z,value);
	  else
		  mpz_set_si(z,getInt32());
  }
  /**
   * Returns a value which is useful for computing hash functions.
   */
  signed long int hashValue()const
  {
	  if(hasLimbs())
		  return mpz_get_si(value);
	  return getInt32();
  }
  bool fitsInInt()const
  {
	  if(hasLimbs())
	  {
		mpz_t v;
		mpz_init(v);
		this->setGmp(v);
		bool ret=(mpz_fits_sint_p(v)!=0);
		mpz_clear(v);
		return ret;
	  }
	  return true;
  }
  int64_t toInt()const//works for 64 bit integers
  {
	  if(hasLimbs())
	  {
		mpz_t v;
		mpz_init(v);
		this->setGmp(v);
		int ret=0;
	//    if(mpz_fits_sint_p(v))
		  ret=mpz_get_si(v);
	//    else
	//      ok=false;
		mpz_clear(v);
		return ret;
	  }
	  return getInt32();
  }
};


#if 0 // enable debug integer here
  class IntegerBoth
{
	  IntGMP A;
	  Integer2 B;
	  bool test()const
	  {
		  mpz_t a,b;
		  mpz_init(a);
		  mpz_init(b);
		  A.setGmp(a);
		  B.setGmp(b);
		  bool ret=(mpz_cmp(a,b)==0);
		  mpz_clear(b);
		  mpz_clear(a);
		  return ret;
	  }
public:
  typedef IntegerBoth POD;//Not actually a POD type. But here we can specify a POD type if one exists.
  static bool isField()
  {
    return false;
  }
  IntegerBoth():
	  A(),B()
  {
	  assert(test());
  }
  IntegerBoth(signed long int value_):
	  A(value_),
	  B(value_)
  {
	  assert(test());
  }
  IntegerBoth(IntegerBoth const & value_):
	  A(value_.A),
	  B(value_.B)
  {
	  assert(test());
  }
  IntegerBoth(mpz_t const value_):
	  A(value_),
	  B(value_)
  {
	  assert(test());
  }
  IntegerBoth(std::string const&s):
	  A(s),
	  B(s)
  {
	  assert(test());
  }

  bool isZero()const{
	  assert(test());
	  assert(A.isZero()==B.isZero());
	  return A.isZero();
  }
  friend std::ostream &operator<<(std::ostream &f, IntegerBoth const &a)
  {
	  assert(a.test());
	  return f<<a.B;
  }
  IntegerBoth& operator+=(const IntegerBoth& a)
    {
	  A+=a.A;
	  B+=a.B;
	  assert(test());
      return *this;
    }
  IntegerBoth& operator-=(const IntegerBoth& a)
    {
	  A-=a.A;
	  B-=a.B;
	  assert(test());
      return *this;
    }
  IntegerBoth& operator*=(const IntegerBoth& a)
    {
	  A*=a.A;
	  B*=a.B;
	  assert(test());
	  return *this;
    }
  IntegerBoth& operator/=(const IntegerBoth& a)
	{
	  assert(a.test());
	  assert(test());
//	  std::cerr<<*this<<"\n";
//	  std::cerr<<a<<"\n";
	  A/=a.A;
	  B/=a.B;
//	  std::cerr<<*this<<"\n";
	  assert(test());
      return *this;
    }
  friend IntegerBoth operator-(const IntegerBoth &b)
  {
    IntegerBoth ret;
    assert(ret.test());
    ret-=b;
    assert(ret.test());
    return ret;
  }
  IntegerBoth operator+(const IntegerBoth &a)const
  {
    IntegerBoth ret(*this);
    assert(ret.test());
    ret+=a;
    assert(ret.test());
    return ret;
  }
  IntegerBoth operator-(const IntegerBoth &a)const
  {
    IntegerBoth ret(*this);
    assert(ret.test());
    ret-=a;
    assert(ret.test());
    return ret;
  }
  IntegerBoth operator*(const IntegerBoth &a)const
  {
    IntegerBoth ret(*this);
    assert(ret.test());
    ret*=a;
    assert(ret.test());
    return ret;
  }
  IntegerBoth operator/(const IntegerBoth &a)const
  {
    IntegerBoth ret(*this);
    assert(ret.test());
    ret/=a;
    assert(ret.test());
    return ret;
  }
  void madd(const IntegerBoth &a,const IntegerBoth &b)
    {
	    assert(test());
	    assert(a.test());
	    assert(b.test());
	  A.madd(a.A,b.A);
	  B.madd(a.B,b.B);
	  assert(test());
    }
  bool operator<(const IntegerBoth &a)const
  {
	  assert((A<a.A)==(B<a.B));
	  return A<a.A;
  }
  bool operator==(const IntegerBoth &a)const
  {
	  assert((A==a.A)==(B==a.B));
	  return A==a.A;
  }
  bool operator!=(const IntegerBoth &a)const
  {
	  return !(*this==a);
  }
  int sign()const
  {
	  assert(A.sign()==B.sign());
	  return A.sign();
  }
  static IntegerBoth gcd(IntegerBoth const &a, IntegerBoth const &b, IntegerBoth &s, IntegerBoth &t)
  {
	  IntegerBoth ret;
	    assert(a.test());
	    assert(b.test());
	  ret.A=IntGMP::gcd(a.A,b.A,s.A,t.A);
	  ret.B=Integer2::gcd(a.B,b.B,s.B,t.B);
	  assert(ret.test());
	  return ret;
  }
  /**
   * Assigns the value to z. z must have been initialized as a gmp variable.
   */
  void setGmp(mpz_t z)const
  {
	  A.setGmp(z);
	  {
		  mpz_t z2;
		  mpz_init(z2);
		  B.setGmp(z2);
		  assert(mpz_cmp(z,z2)==0);
		  mpz_clear(z2);
	  }
  }
  /**
   * Returns a value which is useful for computing hash functions.
   */
  signed long int hashValue()const
  {
	  assert(A.hashValue()==B.hashValue());
	  return A.hashValue();
  }
  bool fitsInInt()const
  {
	  assert(A.fitsInInt()==B.fitsInInt());
	  return A.fitsInInt();
  }
  int64_t toInt()const//works for 64 bit integers
  {
	  assert(A.toInt()==B.toInt());
	  return A.toInt();
  }
};

typedef IntegerBoth Integer;
#else
typedef Integer2 Integer;
#endif
#endif

}

#else
namespace gfan{
  typedef signed long int word;//processor type for integers
  typedef int64_t LimbWord;//memory word. Assumed to be at least as big as word
  typedef uint64_t uLimbWord;//memory word. Assumed to be at least as big as word
  const int limbSizeInBytes=8;
#include <stdint.h>
  struct spec64malloc{
    int64_t data;
    bool hasLimbs()
    {
      return data&1;
    }
    word fitsMask()
    {
      return 0xc000000000000000;
    }
    void assign(word value)//assuming data fits
    {
      data=value<<1;
    }
    void alloc(int nlimbs)//one limb is added since that will tell the number of remaining limbs
    {
      int64_t temp=(int64_t)malloc((nlimbs+1)*limbSizeInBytes);
      assert(temp);
      data=temp+1;
    }
    LimbWord *limbs()//assuming that limbs exist
    {
      return (LimbWord*)(data-1);
    }
    word nlimbs()//assuming limbs exist
    {
      return (word)*limbs();
    }
    void copy(spec64malloc b)
    {
      if(hasLimbs())
        {
          word n2=b.nlimbs()+1;
          int64_t temp=(int64_t)malloc((n2)*limbSizeInBytes);
          assert(temp);

          data=temp+1;
          memcpy((LimbWord*)temp,limbs(),n2*limbSizeInBytes);
        }
      else
        {
          data=b.data;
        }
    }
    void doFree()//assuming that limbs exist
    {
      free((void*)(data-1));
    }
    word valueToWord()//assume no limbs and that word is big enough to contain int64_t
    {
      return data>>1;
    }
  };
  template <struct spec> struct IntegerTemplate : public spec
  {
  private:
    bool fits(word v)
    {
      return !((value_&fitsMask())^((value_<<1)&fitsMask));
    }
  public:
    static bool isField()
    {
      return false;
    }
    IntegerTemplate()
    {
      spec.data=0;
    }
    void assignWordNoFree()
    {
      if(fits(value_))
        {
          assign(value);
        }
      else
        {
          alloc(1);
          limbs()[0]=1;
          limbs()[1]=value_;
        }
    }
    IntegerTemplate(word value_)
    {
      assignWordNoFree(value_);
    }
    IntegerTemplate(IntegerTemplate const & value_)
    {
      if(value_.hasLimbs())
        {
          copy(value_);
        }
      else
        data=value.data;
    }
/*    Integer(mpz_t value_)
    {
      mpz_init_set(value,value_);
    }*/
    ~IntegerTemplate()
    {
      if(hasLimbs())doFree();
    }
    IntegerTemplate& operator=(const IntegerTemplate& a)
      {
        if(this!=&a)
          {
            if(hasLimps())doFree();
            copy(a);
          }
        return *this;
      }
    bool isZero()const{
      return data==0;
    }
    friend std::ostream &operator<<(std::ostream &f, IntegerTemplate const &a)
    {
      if(hasLimps())
        {
          LimpWord *p=limbs();
          int l=*p++;
          for(int i=0;i<l;i++)
            f<<":"<<p[i];
        }
      else
        {
          f<<valueToWord();
        }
      return f;
    }
    LimbWord signExtension(LimbWord a)
    {
      return 0-(a<0);
    }
    void addWordToLimbs(word v)
    {
      int n=nlimbs();
      LimbWord V=v;
      LimbWord *p=limbs()+1;
      LimbWord s=signExtension(V);
      for(int i=0;i<n;i++)
        {
          LimbWord r=V+*p;
          bool carry=(uLimbWord)r<(uLimbWord)V;
          V=carry+s;
        }

    }


    IntegerTemplate& operator+=(const IntegerTemplate& a)
      {
        if(hasLimbs()||a.hasLimbs())
          {
            if(!hasLimbs())
              {

              }
            else
              {
              }

          }
        else
          {
            word C=valueToWord();
            word A=a.valueToWord();
            word D=A+C;
            assignWordNoFree(D);
          }
        return *this;
      }

  };

  typedef IntegerTemplate<spec64malloc> Integer;
};
#endif

#endif /* LIB_Z_H_ */

