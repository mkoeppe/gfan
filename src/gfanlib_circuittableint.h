/*
 * gfanlib_circuittabletypeint.h
 *
 *  Created on: Apr 10, 2016
 *      Author: anders
 */

#ifndef GFANLIB_CIRCUITTABLEINT_H_
#define GFANLIB_CIRCUITTABLEINT_H_

#include <cstdint>
#include <exception>
#include <sstream>
#include <limits>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include "gfanlib_frequencytable.h"

namespace gfan{


  /* We would like to use std::make_unsigned<> but that does not necessarily
     work for 128 bit integers. Therefore we define: */
  template<typename> struct MyMakeUnsigned;
  template <> struct MyMakeUnsigned<int>{typedef unsigned int type;};
  template <> struct MyMakeUnsigned<long int>{typedef unsigned long int type;};
  template <> struct MyMakeUnsigned<__int128>{typedef unsigned __int128 type;};

  class MVMachineIntegerOverflow: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Exception: Overflow (or possible future overflow) in 32/64 bit integers in tropical homotopy.";
  }
};

extern MVMachineIntegerOverflow MVMachineIntegerOverflow;

// The following could be moved to circuittableinteger.h, but there is no object file
class IntegerConversionException: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Exception: Integer conversion exception.";
  }
};

extern IntegerConversionException IntegerConversionException;

// If there is support for operator<< on 128 bit values, then the following 6 procedures can be simplified.

static std::string toStr(__int128_t b)
{
	std::stringstream s;
	s<<"["<<std::setfill('0')<<std::setw(16)<<std::hex<<(unsigned long int)(b>>64)<<":"<<std::setfill('0')<<std::setw(16)<<std::hex<<(unsigned long int)b<<"]";
	return s.str();
}

static std::string toStr(__uint128_t b)
{
	std::stringstream s;
	s<<"["<<std::setfill('0')<<std::setw(16)<<std::hex<<(unsigned long int)(b>>64)<<":"<<std::setfill('0')<<std::setw(16)<<std::hex<<(unsigned long int)b<<"]";
	return s.str();
}

static std::string toStr(__int64_t b)
{
	std::stringstream s;
	s<<b;
	return s.str();
}

static std::string toStr(__uint64_t b)
{
	std::stringstream s;
	s<<b;
	return s.str();
}

static std::string toStr(__int32_t b)
{
	std::stringstream s;
	s<<b;
	return s.str();
}

static std::string toStr(__uint32_t b)
{
	std::stringstream s;
	s<<b;
	return s.str();
}

class my256s{
public:
	__int128_t lo,hi;
	my256s(__int128_t lo_,__int128_t hi_):
		lo(lo_),
		hi(hi_)
	{
	}
	my256s(__int128_t v):
		lo(v),
		hi(0)
	{
		if(v<0)hi=-1;
	}
	my256s operator+(my256s b)
	{
		__uint128_t newLo=lo+b.lo;
		bool carry=newLo<__uint128_t(lo);
		my256s r(newLo,hi+b.hi+carry);
//		std::cerr<<"sum"<<toStr(*this)<<toStr(b)<<"equals"<<toStr(r)<<"\n";
		return r;
	}
	friend std::string toStr(my256s const &m)
//	string toStr()const
	{
		return ::gfan::toStr(m.hi)+::gfan::toStr(m.lo);
	}
	my256s operator/(__int128_t b)const
	{
		std::cerr<<"Dividing"<<toStr(*this)<<"by"<<::gfan::toStr(b)<<"\n";
		assert(0);
		return *this;
	}
	my256s operator<<(int i)const
	{
		if(i==128)return my256s(0,lo);
		if(!(i>=0 && i<=127))std::cerr<<"i"<<i<<"\n";
		assert(i>=0 && i<=127);
		__uint128_t carry=i?(__uint128_t(lo))>>(128-i):0;
		return my256s(lo<<i,(hi<<i)+carry);
	}
	my256s operator>>(int i)const
	{
		assert(i>=0 && i<=127);
//		__uint128_t carry=(__uint128_t(hi))<<(128-i);
		__uint128_t carry=i?(__uint128_t(hi))<<(128-i):0;// a shift of 128 is undefined according to standard
//		std::cerr<<"carry"<<toStr(carry)<<"\n";
		my256s r((__uint128_t(lo)>>i)+carry,hi>>i);

//		std::cerr<<"Shift"<<toStr(*this)<<">>"<<i<<"="<<toStr(r)<<"\n";
		return r;
	}
	bool operator<=(my256s b)const
		{
			if(hi<b.hi)return true;
			if(hi>b.hi)return false;
			return ((__uint128_t)lo)<=((__uint128_t)b.lo);
		}
	bool operator<(my256s b)const
		{
			if(hi<b.hi)return true;
			if(hi>b.hi)return false;
			return ((__uint128_t)lo)<((__uint128_t)b.lo);
		}
	bool operator<=(__int128_t b)const
		{
		return *this<=my256s(b);
		}
	my256s operator-()const
	{
		my256s temp(~lo,~hi);
		return temp+my256s(1,0);
	}
	explicit operator __int128()const
	{
		return lo;
	}
};

class my256u{
public:
	__int128_t lo,hi;
	my256u(my256s const &a):
		lo(a.lo),
		hi(a.hi)
	{
	}
	friend std::string toStr(my256u const &m)
	{
		return ::gfan::toStr(m.hi)+::gfan::toStr(m.lo);
	}
	bool operator<=(my256s const &b)//full implementation not required as this simulates a cast to unsigned and a comparison.
	{
		if((__uint128_t)hi<(__uint128_t)b.hi)return true;
		if((__uint128_t)hi>(__uint128_t)b.hi)return false;
		return ((__uint128_t)lo)<=((__uint128_t)b.lo);
	}
	bool operator<(my256s const &b)//full implementation not required as this simulates a cast to unsigned and a comparison.
	{
		if((__uint128_t)hi<(__uint128_t)b.hi)return true;
		if((__uint128_t)hi>(__uint128_t)b.hi)return false;
		return ((__uint128_t)lo)<((__uint128_t)b.lo);
	}
	my256u operator<<(int i)const
	{
		assert(i>=0 && i<=127);
		__uint128_t carry=i?(__uint128_t(lo))>>(128-i):0;
		return my256s(lo<<i,__int128_t((hi<<i)+carry));
	}
};

template <> struct MyMakeUnsigned<my256s>{typedef my256u type;};

static long int extMul(int a, int b)
{
	return ((long int)a)*((long int)b);
}
static __int128_t extMul(long int a, long int b)
{
	return ((__int128_t)a)*((__int128_t)b);
}

static __uint128_t unsignedProd64(uint64_t x,uint64_t y)
{
	return __uint128_t(x)*__uint128_t(y);
}

static my256u unsignedProd128(__uint128_t x,__uint128_t y)
{
	my256s a(unsignedProd64(x,y),0);
	my256s b(unsignedProd64(x>>64,y),0);
	my256s c(unsignedProd64(x,y>>64),0);
	my256s d(unsignedProd64(x>>64,y>>64),0);
	return a+(b<<64)+(c<<64)+(d<<128);
}

static my256s extMul(__int128_t a, __int128_t b)
{
//	std::cerr<<"mul"<<toStr(a)<<"*"<<toStr(b)<<"\n";
	auto temp=unsignedProd128(a,b);
	my256s r(temp.lo,temp.hi);
	if(a<0)r=r+my256s(0,-b);
	if(b<0)r=r+my256s(0,-a);
//	std::cerr<<"result"<<toStr(r)<<"\n";
	return my256s(r);
}

/*
 * The philosophy here is that if this class overflows, then the computation needs to be restarted. Therefore
 * all overflows must be caught.
 */
/*
 * Remark: Really we would like to have just one one type for CircuitTableInt32POD (with constructors).
 * Unfortunately if the class is not POD, gcc inserts n null pointer check in the constructor of the class
 * when constructing vector<CircuitTableInt32POD> of size n. These checks are not optimised away. On the other
 * hand vector will initialise POD to zero, leading to more efficient code. There seems to be no reason for
 * why constructors force this behaviour. Hopefully it will be fixed in the next C++ standard.
 * For now we define CircuitTableInt32POD as a subclass of CircuitTableInt32PODPOD exploiting slicing.
 */
template<typename word, typename longword> class CircuitTableIntPOD{
 public:
	class Divisor{
	public:
		word v;
		int shift;
		word multiplicativeInverse;
		Divisor(CircuitTableIntPOD const &a)// for exact division
		{ // A good place to read about these tricks seems to be the book "Hacker's Delight" by Warren.
			v=a.v;
			shift=0;
			word t=v;
			assert(t);
			while(!(t&1)){t>>=	1;shift++;}
			word inverse=t;
			while(t*inverse!=1)inverse*=2-t*inverse;
			multiplicativeInverse=inverse;

//			std::cerr<<"DIVOBJ:v"<<toStr(v)<<"s"<<shift<<"inv"<<toStr(multiplicativeInverse)<<"\n";
		}
	};
	class Double{
	public:
		longword v;
		Double():v(0){};
		Double(longword a):v(a){};
		Double &operator+=(Double a){v+=a.v;return *this;}
		Double &operator-=(Double a){v-=a.v;return *this;}
		CircuitTableIntPOD castToSingle()const;
		bool isZero()const{return v==0;}
		bool isNegative()const{return v<0;}
		Double operator-()const{Double ret;ret.v=-v;return ret;}
		friend Double operator-(Double const &a, Double const &b){return Double(a.v-b.v);}
		friend Double operator+(Double const &a, Double const &b){return Double(a.v+b.v);}
		Double &addWithOverflowCheck(Double const &a)
		{
			longword bound=longword{1}<<(std::numeric_limits<longword>::digits-2);
			if(a.v<0 || v<0 || a.v>=bound || v>=bound) throw MVMachineIntegerOverflow;
			v+=a.v;
			return *this;
		}
		std::string toString()const{std::stringstream s;s<<v;return s.str();}
	};
 private:
public:
	word v;
	friend CircuitTableIntPOD operator+(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){CircuitTableIntPOD ret;ret.v=a.v+b.v;return ret;}
	friend CircuitTableIntPOD operator-(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){CircuitTableIntPOD ret;ret.v=a.v-b.v;return ret;}
	friend CircuitTableIntPOD operator*(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){CircuitTableIntPOD ret;ret.v=a.v*b.v;return ret;}
	friend CircuitTableIntPOD operator/(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){CircuitTableIntPOD ret;ret.v=a.v/b.v;return ret;}//This is used very few times. Should we require this be an exact division?
 public:
	static const word halfBound{(word{1}<<(std::numeric_limits<word>::digits/2-2))-1};
	// In the code products of CircuitTableIntPOD objects are summed. To avoid overflows one of the factors must "fit in half" and the "number of summands" may not exceed a certain number. The bounds are specified in the following functions:
	bool fitsInHalf()const{return v>-halfBound && v<halfBound;}// Is it better to allow only non-negative numbers?
	static bool isSuitableSummationNumber(int numberOfSummands){return numberOfSummands<halfBound;}
//	CircuitTableIntPOD(CircuitTableIntPOD&& other)noexcept:v(other.v){}
// 	class CircuitTableIntPOD &operator=(CircuitTableIntPOD a){v=a.v;return *this;}
  	class CircuitTableIntPOD &operator=(word a){v=a;return *this;}
	CircuitTableIntPOD operator-()const{CircuitTableIntPOD ret;ret.v=-v;return ret;}
	CircuitTableIntPOD &operator-=(CircuitTableIntPOD a){v-=a.v;return *this;}
	CircuitTableIntPOD &operator+=(CircuitTableIntPOD a){v+=a.v;return *this;}
	CircuitTableIntPOD &operator*=(CircuitTableIntPOD a){v*=a.v;return *this;}
	friend bool operator<(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){return a.v<b.v;}
	friend bool operator<=(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){return a.v<=b.v;}
	friend bool operator==(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){return a.v==b.v;}
	bool isZero()const{return v==0;}
	bool isOne()const{return v==1;}
	bool isNonZero()const{return v!=0;}
	bool isNegative()const{return v<0;}
	bool isPositive()const{return v>0;}
	int sign()const{return isNegative()?-1:isPositive();}
	friend int determinantSign1(CircuitTableIntPOD const &a, CircuitTableIntPOD const &c, CircuitTableIntPOD const &b, CircuitTableIntPOD const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		longword r=((longword)a.v)*((longword)c.v)-((longword)b.v)*((longword)d.v);
		if(r>0)return 1;
		if(r<0)return -1;
		return 0;
	}
	friend int determinantSign3(CircuitTableIntPOD const &a, CircuitTableIntPOD const &c, CircuitTableIntPOD const &b, CircuitTableIntPOD const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		auto X=extMul(a.v,c.v);
		auto Y=extMul(b.v,d.v);
		if(Y<X)return 1;
		if(X<Y)return -1;
		return 0;
	}
	/**
	 * The following does not compute with full precision. The idea of passing a double is just that for example a sign change before calling this function could cause overflow, but with a double this in not a problem. The computation can be handled in double precission.
	 */
	friend int determinantSign2(CircuitTableIntPOD::Double const &a, CircuitTableIntPOD const &c, CircuitTableIntPOD::Double const &b, CircuitTableIntPOD const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		longword r1=(a.v)*((longword)c.v);
		longword r2=(b.v)*((longword)d.v);
		if(r1>r2)return 1;
		if(r1<r2)return -1;
		return 0;
	}
	std::string toString()const{
	//	return std::to_string((int64_t)v);/*cast seems to be needed (128 bit will not work)*/
		return toStr(v);
	}
	int64_t toInt64()const{return v;}//WHAT SHOULD HAPPEN TO THIS FUNCTION?
  friend std::ostream &operator<<(std::ostream &f, CircuitTableIntPOD const &a){f<</*(int)*/a.toString();return f;}
	Double extend()const{Double ret;ret.v=v;return ret;}
	CircuitTableIntPOD &maddWithOverflowChecking(CircuitTableIntPOD const &a, CircuitTableIntPOD const&b){Double t=this->extend();t+=extendedMultiplication(a,b);*this=t.castToSingle();return *this;}
	CircuitTableIntPOD &msubWithOverflowChecking(CircuitTableIntPOD const &a, CircuitTableIntPOD const&b){Double s=this->extend();s-=extendedMultiplication(a,b);*this=s.castToSingle();return *this;}
	CircuitTableIntPOD &subWithOverflowChecking(CircuitTableIntPOD const &a){Double t=this->extend();t-=a.extend();*this=t.castToSingle();return *this;}
	CircuitTableIntPOD &addWithOverflowChecking(CircuitTableIntPOD const &a){Double t=this->extend();t+=a.extend();*this=t.castToSingle();return *this;}
	CircuitTableIntPOD negatedWithOverflowChecking()const{return (-extend()).castToSingle();}
	friend Double extendedMultiplication(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){Double ret;ret.v=((longword)a.v)*((longword)b.v);return ret;}
	friend Double extendedMultiplication(CircuitTableIntPOD const &a, int32_t b){Double ret;ret.v=((longword)a.v)*((longword)b);return ret;}//to be removed?
	friend CircuitTableIntPOD min(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){return (a.v>=b.v)?b:a;}
	friend CircuitTableIntPOD max(CircuitTableIntPOD const &a, CircuitTableIntPOD const &b){return (a.v<=b.v)?b:a;}
	friend CircuitTableIntPOD negabs(CircuitTableIntPOD const &a){return min(a,-a);}
	friend CircuitTableIntPOD abs(CircuitTableIntPOD const &a){return max(a,a.negatedWithOverflowChecking());}
	friend CircuitTableIntPOD dotDiv(CircuitTableIntPOD const &s, CircuitTableIntPOD const &a, CircuitTableIntPOD const &t, CircuitTableIntPOD const &b, CircuitTableIntPOD::Divisor const &q)
	{
		CircuitTableIntPOD ret;
		ret.v=((((word)((/*(uint64_t)*/(((longword)s.v)*((longword)a.v)+((longword)t.v)*((longword)b.v)))>>q.shift)))*q.multiplicativeInverse);
		return ret;
	}
	friend void dotDivAssign(CircuitTableIntPOD &s, CircuitTableIntPOD const &a, CircuitTableIntPOD const &t, CircuitTableIntPOD const &b, CircuitTableIntPOD::Divisor const &q)//as above but assigning to first argument. This helps the vectorizer.
	{
		s.v=((((word)((/*(uint64_t)*/(((longword)s.v)*((longword)a.v)+((longword)t.v)*((longword)b.v)))>>q.shift)))*q.multiplicativeInverse);
	}
	static word MIN(word a, word b)
	{
		return (a<b)?a:b;
	}
	static word MAX(word a, word b)
	{
		return (a>b)?a:b;
	}
	static CircuitTableIntPOD computeNegativeBound(CircuitTableIntPOD * __restrict__ Ai, int w)
	{
		CircuitTableIntPOD M{};assert(M.v==0);//POD types are zero initialised if value initialisation is used. The braces force value initialisation  and not default initialisation, which would be not do anything.
		CircuitTableIntPOD m{};assert(m.v==0);
		for(int j=0;j<w;j++)
		{
			m.v=MIN(m.v,Ai[j].v);
			M.v=MAX(M.v,Ai[j].v);
		}
//		std::cerr<<"BOUND COMP:"<<min(m,-M).toString()<<"\n";
		return min(m,-M);
	}
	static CircuitTableIntPOD quickNoShiftBounded(CircuitTableIntPOD * __restrict__ a, CircuitTableIntPOD * __restrict__ b, CircuitTableIntPOD s, CircuitTableIntPOD t, CircuitTableIntPOD::Divisor denominatorDivisor,int c)
	{
		CircuitTableIntPOD *aa=a;
		CircuitTableIntPOD *bb=b;
		CircuitTableIntPOD max{};assert(max.v==0);
	    CircuitTableIntPOD min{};assert(min.v==0);
	    CircuitTableIntPOD ret{};assert(ret.v==0);
	    for(int i=0;i<c;i++)
	      {
	    	aa[i].v=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v+
	    			(t.v*denominatorDivisor.multiplicativeInverse)*bb[i].v);
	    	min.v=MIN(min.v,aa[i].v);
	    	max.v=MAX(max.v,aa[i].v);
	      }
	    if(-max<=min)min=-max;
	    ret=min;
	    return ret;
	}
	static CircuitTableIntPOD dotDivVector(CircuitTableIntPOD * __restrict__ aa, CircuitTableIntPOD * __restrict__ bb, CircuitTableIntPOD s, CircuitTableIntPOD t, CircuitTableIntPOD::Divisor denominatorDivisor,int c, CircuitTableIntPOD boundA, CircuitTableIntPOD boundB)__attribute__ ((always_inline))
			//assigning to first argument. The return negative of the return value is an upper bound on the absolute values of the produced entries
	{
//	      std::cerr<<"ENTDOTDIV\n";
//		{static int i;i++;if(i>100000000){i=0;std::stringstream S;S<<"c"<<c<<"s"<<s<<"t"<<t<<"v"<<denominatorDivisor.v<<std::string(s.isNonZero()?"\t *":"")<<"\n";std::cerr<<S.str();}}

		while(((s.v|t.v)&1)==0 && denominatorDivisor.shift>0){s.v>>=1;t.v>>=1;denominatorDivisor.shift--;denominatorDivisor.v>>=1;}

		CircuitTableIntPOD max{};assert(max.v==0);
		CircuitTableIntPOD min{};assert(min.v==0);

	  // This is a bound on the absolute value of the sums of products before dividing by the denominator
	  // The calculation can overflow, but since we are casting to unsigned long, this is not a problem.
		typename MyMakeUnsigned<longword>::type positiveResultBoundTimesD=extMul(negabs(s).v,boundA.v)+extMul(negabs(t).v,boundB.v);

	  /* The first case is the case where we know that the intermediate results fit in 32bit before shifting down.
	   * In other words, the shifted-in bits of the results are equal.
	   * In that case the intermediate result can be computed with 32bits.
	   */
	  /*
	   * This is a discussion of what happens if word refers to a 64 bit integer on avx:
	   * In this case the loop will not be vectorised (although it may happen on avx2).
	   * This is because avx only has 128 bit vector instructions for multiplying 32 bit
	   * ints and not 256 bit vector instructions. One therefore needs 3 32 bit vector
	   * multiplications for doing 4 64*64->64 bit multiplications. Then it is equally
	   * efficient to just do 4 usual 64*64->64 bit multiplications without vector
	   * instructions. In most cases, however, either the scalar, the bound or both will
	   * fit in 32 bit. In such cases only one or two of the three 32 bit vector
	   * multiplications will be needed - leading to much faster code. Therefore these
	   * should be treated as separate cases in the future. This discussion ignores the
	   * inverses that would usually not fit in 32 bit! So in fact, it might be better to
	   * check if the result is known to fit in 32 bits. In that case only the lowest 32
	   * bit need to be computed - the remaining are obtained via sign extension. It seems
	   * that avx does not have instructions for shrinking vectors from 64 bit entries
	   * to 32 bit entries. And it would make little sense to waste memory in this way.
	   * Maybe the entire matrix should be 32 bit then.
	   *
	   * Note: On avx2 8 32 bit multiplications can be done with instruction vpmulld.
	   * Note: On avx512vl and avx512dq 64 bit multiplications are vectorised with
	   *       instruction vpmullq (up to 8 multiplications per instruction on avx512dq).
	   */
	  /*
	   * The old bound below seems to be wrong. The following example illustrates the problem:
	   * Division of 300 as a 16 bit number by 20 with the result being 8 bit.
	   * shift=2, inverse of 5 is 205. 300=$12c
	   * ($2c>>2)*205=11*205=2255=$8cf
	   * Here the 6 lower bits are correct, but c is not.
	   * Of course this is not going to pass the original test because we test  300<32*5 which is false.
	   * But the power digits-2=5 seems arbitrary anyway because digits is 7 and not 8.
	   * It is probably possible to construct an example where the code fails.
	   */
//		if(positiveResultBoundTimesD<(((longword{1}<<(std::numeric_limits<word>::digits-2))*denominatorDivisor.v)>>denominatorDivisor.shift))//check this carefully ORIGINAL
		if(positiveResultBoundTimesD<=(longword)std::numeric_limits<word>::max())// This seems to be the correct bound.
	  {//FAST VERSION
//		  std::cerr<<s.v<<" "<<t.v<<" "<<denominatorDivisor.v<<" "<<denominatorDivisor.shift<<"\n";
//			std::cerr<<"FAST\n";
		  if(denominatorDivisor.shift==0)return quickNoShiftBounded(aa,bb,s,t,denominatorDivisor,c);
//		  std::cerr<<"TYPE1\n";
		  for(int i=0;i<c;i++)
		  {
			  aa[i].v=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v+
					  (t.v*denominatorDivisor.multiplicativeInverse)*bb[i].v)>>denominatorDivisor.shift;
			  min.v=MIN(min.v,aa[i].v);
			  max.v=MAX(max.v,aa[i].v);
		  }
		  if(-max<=min)min=-max;
//	      std::cerr<<"RETDOTDIV\n";
		  return min;
	  }
	  else
	  {
		  /* The next case would be to check if the result is known to fit in 32bit.
		   * In that case intermediate results would be handled in 64 bit,
		   * but the final result would not have to be checked for 32 bit overflows.
		   */
//		  if(positiveResultBoundTimesD<(((longword{1}<<(std::numeric_limits<word>::digits-2))*denominatorDivisor.v)))
		  int Dig=std::numeric_limits<word>::digits;
		  if(Dig==0){Dig=127;}//fixes bug in gcc-8.1
		  if(positiveResultBoundTimesD<((longword(denominatorDivisor.v)<<(Dig-2))))
		  {
			  if(0){
				  static int a;
				  a++;
				  if(a%1024==0)
			  if(denominatorDivisor.shift){
				  std::cerr<<toStr(s.v)<<" "<<toStr(t.v)<<" "<<toStr(denominatorDivisor.v)<<"\n";
			  }
			  }
//			  std::cerr<<"TYPE2\n";
//			  std::cerr<<denominatorDivisor.shift<<" "<<toStr(denominatorDivisor.multiplicativeInverse)<<"div"<<toStr(denominatorDivisor.v)<<"\n";
//			  std::cerr<<"s"<<toStr(s.v)<<"t"<<toStr(t.v)<<"\n";
//			  std::cerr<<"aa"<<toStr(aa[0].v)<<"bb"<<toStr(bb[0].v)<<"\n";
			  if(denominatorDivisor.shift)
				  for(int i=0;i<c;i++)
				  {
	//				  aa[i].v=((((word)((/*(uint64_t)*/(((longword)s.v)*((longword)aa[i].v)+((longword)t.v)*((longword)bb[i].v)))>>denominatorDivisor.shift)))*denominatorDivisor.multiplicativeInverse);
					  aa[i].v=((((word)((/*(uint64_t)*/(extMul(s.v,aa[i].v)+extMul(t.v,bb[i].v)))>>denominatorDivisor.shift)))*denominatorDivisor.multiplicativeInverse);
					  if(max<=aa[i])max=aa[i];
					  if(aa[i]<=min)min=aa[i];
				  }
			  else
				  for(int i=0;i<c;i++)
				  {
	//				  aa[i].v=((((word)((/*(uint64_t)*/(((longword)s.v)*((longword)aa[i].v)+((longword)t.v)*((longword)bb[i].v)))>>denominatorDivisor.shift)))*denominatorDivisor.multiplicativeInverse);
					  aa[i].v=(s.v*aa[i].v+t.v*bb[i].v)*denominatorDivisor.multiplicativeInverse;
					  if(max<=aa[i])max=aa[i];
					  if(aa[i]<=min)min=aa[i];
				  }
			  if(-max<=min)min=-max;
		  }
		  /* The last case would be to where we don't know that the results will fit in 32 bit.
		   * Since we need to compute max and min anyway, we can compute these quantities in 64bit
		   * and just check if they fit in 32bit at the end.
		   */
		  else
		  {
			  if(0){
				  std::cerr<<toStr(positiveResultBoundTimesD)<<"\n"<<toStr(denominatorDivisor.v)<<"\n";
				  std::cerr<<"---\n";
				  std::cerr<<toStr(negabs(t).v)<<"\n";
				  std::cerr<<toStr(boundA.v)<<"\n";
				  std::cerr<<toStr(extMul(negabs(t).v,boundA.v))<<"\n";
				  std::cerr<<toStr(negabs(s).v)<<"\n";
				  std::cerr<<toStr(boundB.v)<<"\n";
				  std::cerr<<toStr(extMul(negabs(s).v,boundB.v))<<"\n";
				  std::cerr<<"---\n";

//				  std::cerr<<this->)
			  }
			  int D=std::numeric_limits<word>::digits;
			  if(D==0){D=127;}//fixes bug in gcc-8.1
			  bool doesOverflow=(((word)t.v)==(word{1}<<(D-1)));// What is the purpose of this line. Do we really want to subtract 1? That seems wrong since word is signed. Look at comment below
			  longword min64=0;
			  longword max64=0;
			  for(int i=0;i<c;i++)
			  {
				  // In the unlikely event that all the numbers to be multiplied are -2^31, the following code will overflow. Therefore the overflow flag was set as above.
//				  longword/*int64_t*/ temp=(((longword)s.v)*((longword)aa[i].v)+((longword)t.v)*((longword)bb[i].v))/denominatorDivisor.v;
				  longword/*int64_t*/ temp=(extMul(s.v,aa[i].v)+extMul(t.v,bb[i].v))/denominatorDivisor.v;
				  if(max64<=temp)max64=temp;
				  if(temp<=min64)min64=temp;
				  aa[i].v=word(temp);
			  }
			  if(-max64<=min64)min64=-max64;
			  if(min64<=-std::numeric_limits<word>::max())doesOverflow=true;
			  if(doesOverflow)throw MVMachineIntegerOverflow;
			  min=word(min64);
		  }
	  }
//      std::cerr<<"RETDOTDIV\n";
	  return min;
	}
	static CircuitTableIntPOD quickScaleNoShiftBounded(CircuitTableIntPOD * __restrict__ aa, CircuitTableIntPOD s, CircuitTableIntPOD::Divisor denominatorDivisor,int c, CircuitTableIntPOD boundA)
	{
//		assert(!boundA.isPositive());
		if(0)
		{
			static FrequencyTable A("Multiplier");
			A.record(s.toInt64()&255);
			static FrequencyTable B("Divisor");
			B.record(denominatorDivisor.v&255);
			static FrequencyTable C("AreEqual");
			C.record(s.toInt64()==denominatorDivisor.v);
		}
		CircuitTableIntPOD max{};assert(max.v==0);
	    CircuitTableIntPOD min{};assert(min.v==0);
	    CircuitTableIntPOD ret{};assert(ret.v==0);
	    for(int i=0;i<c;i++)
	      {
	    	aa[i].v=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v);
	    	min.v=MIN(min.v,aa[i].v);
	    	max.v=MAX(max.v,aa[i].v);
	      }
	    if(-max<=min)min=-max;
	    ret=min;

	    {
	    	auto ret2=CircuitTableIntPOD((s.v*denominatorDivisor.multiplicativeInverse)*boundA.v);
	    	ret2=MIN(ret2.v,-ret2.v);
	    	return ret2;
	    	if(ret.v!=ret2.v)
	    	{
	    		std::cerr<<"ret:"<<ret.toString()<<"\tret2:"<<ret2.toString()<<"\n";
	    		assert(0);
	    	}
	    }

	    return ret;
//	    return CircuitTableIntPOD((s.v*denominatorDivisor.multiplicativeInverse)*boundA.v);
	}
	static CircuitTableIntPOD scaleVector(CircuitTableIntPOD * __restrict__ aa, CircuitTableIntPOD s, CircuitTableIntPOD::Divisor denominatorDivisor,int c, CircuitTableIntPOD boundA)__attribute__ ((always_inline))
			//scales the vector aa by s and divides by the divisor. The negative of the return value is an upper bound on the absolute values of the produced entries
	{
		while(((s.v)&1)==0 && denominatorDivisor.shift>0){s.v>>=1;denominatorDivisor.shift--;denominatorDivisor.v>>=1;}

		CircuitTableIntPOD max{};assert(max.v==0);
		CircuitTableIntPOD min{};assert(min.v==0);

		// This is a bound on the absolute value of the product before dividing by the denominator
		// The calculation can overflow, but since we are casting to unsigned long, this is not a problem.
		typename MyMakeUnsigned<longword>::type positiveResultBoundTimesD=extMul(negabs(s).v,boundA.v);

		if(positiveResultBoundTimesD<=(longword)std::numeric_limits<word>::max())
		{//FAST VERSION
		  if(denominatorDivisor.shift==0)return quickScaleNoShiftBounded(aa,s,denominatorDivisor,c,boundA);
		  for(int i=0;i<c;i++)
		  {
			  aa[i].v=((s.v*denominatorDivisor.multiplicativeInverse)*aa[i].v)>>denominatorDivisor.shift;
			  min.v=MIN(min.v,aa[i].v);
			  max.v=MAX(max.v,aa[i].v);
		  }
		  if(-max<=min)min=-max;
		  return min;
		}
	  else
	  {
		  /* The next case would be to check if the result is known to fit in 32bit.
		   * In that case intermediate results would be handled in 64 bit,
		   * but the final result would not have to be checked for 32 bit overflows.
		   */
		  int Dig=std::numeric_limits<word>::digits;
		  if(Dig==0){Dig=127;}//fixes bug in gcc-8.1
		  if(positiveResultBoundTimesD<((longword(denominatorDivisor.v)<<(Dig-2))))
		  {
			  if(denominatorDivisor.shift)
				  for(int i=0;i<c;i++)
				  {
					  aa[i].v=((((word)((/*(uint64_t)*/(extMul(s.v,aa[i].v)))>>denominatorDivisor.shift)))*denominatorDivisor.multiplicativeInverse);
					  if(max<=aa[i])max=aa[i];
					  if(aa[i]<=min)min=aa[i];
				  }
			  else
				  for(int i=0;i<c;i++)
				  {
					  aa[i].v=(s.v*aa[i].v)*denominatorDivisor.multiplicativeInverse;
					  if(max<=aa[i])max=aa[i];
					  if(aa[i]<=min)min=aa[i];
				  }
			  if(-max<=min)min=-max;
		  }
		  /* The last case would be to where we don't know that the results will fit in 32 bit.
		   * Since we need to compute max and min anyway, we can compute these quantities in 64bit
		   * and just check if they fit in 32bit at the end.
		   */
		  else
		  {
			  int D=std::numeric_limits<word>::digits;
			  if(D==0){D=127;}//fixes bug in gcc-8.1
			  bool doesOverflow=false;//(((word)t.v)==(word{1}<<(D-1)));// What is the purpose of this line? t is not defined. Do we really want to subtract 1? That seems wrong since word is signed. Look at comment below
			  longword min64=0;
			  longword max64=0;
			  for(int i=0;i<c;i++)
			  {
				  longword/*int64_t*/ temp=(extMul(s.v,aa[i].v))/denominatorDivisor.v;
				  if(max64<=temp)max64=temp;
				  if(temp<=min64)min64=temp;
				  aa[i].v=word(temp);
			  }
			  if(-max64<=min64)min64=-max64;
			  if(min64<=-std::numeric_limits<word>::max())doesOverflow=true;
			  if(doesOverflow)throw MVMachineIntegerOverflow;
			  min=word(min64);
		  }
	  }
	  return min;
	}
  static bool isField() // This is here to make old gcd computations work
  {
    return false;
  }
#if 0
  static CircuitTableIntPOD gcd(CircuitTableIntPOD a, CircuitTableIntPOD b)
  {
    return ::gcd(a,b);
  }
  #endif
  static CircuitTableIntPOD gcd(CircuitTableIntPOD a, CircuitTableIntPOD b)
  {
    //  std::cerr<<"a"<<a<<"b"<<b<<"\n";
  	if(a.isZero() &&b.isPositive())return b;
  	if(a.isZero() &&b.isNegative())return -b;
  	if(b.isZero() &&a.isPositive())return a;
  	if(b.isZero() &&a.isNegative())return -a;
  	if(a.isNegative())return gcd(-a,- -b);
  	if(b.isNegative())return gcd(- -a,-b);
  	if((b-a).isPositive()){swap(a,b);}
  	while(!b.isZero())
  	{
  		a-=(a/b)*b;
  		swap(a,b);
// 		a-=b;
// 		if((b-a).isPositive()){swap(a,b);}
  	}
  	return a;
  }
  friend CircuitTableIntPOD gcd(CircuitTableIntPOD a, CircuitTableIntPOD b)
  {
    return CircuitTableIntPOD::gcd(a,b);
  }
};




typedef CircuitTableIntPOD<int32_t,int64_t> CircuitTableInt32POD;
typedef CircuitTableIntPOD<int64_t,__int128_t> CircuitTableInt64POD;
typedef CircuitTableIntPOD<__int128_t,my256s> CircuitTableInt128POD;

static bool hasPod(class CircuitTalbeInt32 *c){return true;}
static bool hasPod(class CircuitTalbeInt64 *c){return true;}

template<> inline CircuitTableIntPOD<int32_t,int64_t> CircuitTableIntPOD<int32_t,int64_t>::Double::castToSingle()const//casts and checks precission
{
	if(v>=0x7fffffff || -v>=0x7fffffff) throw MVMachineIntegerOverflow;
	CircuitTableIntPOD ret;
	ret.v=v;
	return ret;
}

template<> inline CircuitTableIntPOD<int64_t,__int128_t> CircuitTableIntPOD<int64_t,__int128_t>::Double::castToSingle()const//casts and checks precission
{
	if(v>=0x7fffffffffffffff || -v>=0x7fffffffffffffff) throw MVMachineIntegerOverflow;
	CircuitTableIntPOD ret;
	ret.v=v;
	return ret;
}

template<> inline CircuitTableIntPOD<__int128_t,my256s> CircuitTableIntPOD<__int128_t,my256s>::Double::castToSingle()const//casts and checks precission
{
  // DANGER !!!
  //if(v>=0x7fffffffffffffffffffffffffffffff || -v>=0x7fffffffffffffffffffffffffffffff) throw MVMachineIntegerOverflow;
	CircuitTableIntPOD ret;
	ret.v=__int128(v);
	return ret;
}

static_assert(std::is_standard_layout<CircuitTableInt32POD>::value,"");
static_assert(std::is_trivial<CircuitTableInt32POD>::value,"");
//static_assert(std::is_pod<CircuitTableInt32POD>::value,"");

class CircuitTableInt32:public CircuitTableInt32POD
{
public:
	typedef CircuitTableInt32POD POD;
	static bool POD2;
	CircuitTableInt32()noexcept{v=0;}
	CircuitTableInt32(CircuitTableInt32POD const &m){v=m.v;}
	CircuitTableInt32(int32_t val){v=val;}
	CircuitTableInt32(std::string const&s){std::istringstream a(s);a>>v;}
};

class CircuitTableInt64:public CircuitTableInt64POD
{
public:
	typedef CircuitTableInt64POD POD;
	static bool POD2;
	CircuitTableInt64()noexcept{v=0;}
	CircuitTableInt64(CircuitTableInt64POD const &m){v=m.v;}
	CircuitTableInt64(int64_t val){v=val;}
	CircuitTableInt64(std::string const&s){std::istringstream a(s);a>>v;}
};

class CircuitTableInt128:public CircuitTableInt128POD
{
public:
	typedef CircuitTableInt128POD POD;
	static bool POD2;
	CircuitTableInt128()noexcept{v=0;}
	CircuitTableInt128(CircuitTableInt128POD const &m){v=m.v;}
	CircuitTableInt128(__int128_t val){v=val;}
	CircuitTableInt128(std::string const&s){
	  int64_t proxy;
	  std::istringstream a(s); a>>proxy;
	  v = proxy;
	  /*assert(0);*/}
};
}


#endif /* GFANLIB_CIRCUITTABLETYPEINT_H_ */
