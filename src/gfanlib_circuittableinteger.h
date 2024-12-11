/*
 * gfanlib_circuittableinteger.h
 *
 *  Created on: 25 Feb 2020
 *      Author: anders
 */

/*
 * This is an implementation of an integer that can replace CircuitTableInt32 for
 * debugging purposes.
 */

#ifndef SRC_GFANLIB_CIRCUITTABLEINTEGER_H_
#define SRC_GFANLIB_CIRCUITTABLEINTEGER_H_


#include <cstdint>
#include <exception>
#include <sstream>
#include <assert.h>

#include "gfanlib_z.h"



namespace gfan{

class CircuitTableInteger{
	Integer v;
 public:
	CircuitTableInteger(Integer const &v_): //protect
		v(v_)
 	 {
 	 }
	CircuitTableInteger() //protect
 	 {
 	 }
	CircuitTableInteger(int32_t val):
		v(val)
	{
	}
        CircuitTableInteger(std::string const&s):  
	        v(s)
	{
	}
	// Explicit casting between CircuitTableInts
	// Inspired by https://stackoverflow.com/questions/714213/c-template-casting
	explicit CircuitTableInteger(CircuitTableInt32 const &a):
		v(a.v)
	{
	}
	operator gfan::CircuitTableInt32()const
	{
		if(!v.fitsInInt())throw IntegerConversionException;
		return v.toInt();
	}
	explicit CircuitTableInteger(CircuitTableInt64 const &a):
		v(a.v)
	{
	}
/*	operator gfan::CircuitTableInt64()const  //Currently functions are missing for implementing this.
	{
		if(!v.fitsInInt())throw IntegerConversionException;
		return v.toInt();
	}*/

	class Divisor{
	public:
		Integer v;
		Divisor(CircuitTableInteger const &a):
			v(a.v)
		{
			assert(!v.isZero());
		}
	};
	class Double{
	public:
		Integer v;
		Double(){};
		Double(Integer a):v(a){};
		Double &operator+=(Double a){v+=a.v;return *this;}
		Double &operator-=(Double a){v-=a.v;return *this;}
		CircuitTableInteger castToSingle()const;
		bool isZero()const{return v.sign()==0;}
		bool isNegative()const{return v.sign()<0;}
		Double operator-()const{Double ret;ret.v=-v;return ret;}
		friend Double operator-(Double const &a, Double const &b){return Double(a.v-b.v);}
		friend Double operator+(Double const &a, Double const &b){return Double(a.v+b.v);}
		Double &addWithOverflowCheck(Double const &a)
		{
			return *this+=a;;
		}
		std::string toString()const{std::stringstream s;s<<v;return s.str();}
	};
 private:
public:
	friend CircuitTableInteger operator+(CircuitTableInteger const &a, CircuitTableInteger const &b){CircuitTableInteger ret;ret.v=a.v+b.v;return ret;}
	friend CircuitTableInteger operator-(CircuitTableInteger const &a, CircuitTableInteger const &b){CircuitTableInteger ret;ret.v=a.v-b.v;return ret;}
	friend CircuitTableInteger operator*(CircuitTableInteger const &a, CircuitTableInteger const &b){CircuitTableInteger ret;ret.v=a.v*b.v;return ret;}
	friend CircuitTableInteger operator/(CircuitTableInteger const &a, CircuitTableInteger const &b){CircuitTableInteger ret;ret.v=a.v/b.v;return ret;}
 public:
	// In the code products of CircuitTableInt32POD objects are summed. To avoid overflows one of the factors must "fit in half" and the "number of summands" may not exceed a certain number. The bounds are specified in the following functions:
	bool fitsInHalf()const{return true;}
	static bool isSuitableSummationNumber(int numberOfSummands){return true;}
  	class CircuitTableInteger &operator=(int32_t a){v=a;return *this;}
	CircuitTableInteger operator-()const{CircuitTableInteger ret;ret.v=-v;return ret;}
	CircuitTableInteger &operator-=(CircuitTableInteger a){v-=a.v;return *this;}
	CircuitTableInteger &operator+=(CircuitTableInteger a){v+=a.v;return *this;}
	CircuitTableInteger &operator*=(CircuitTableInteger a){v*=a.v;return *this;}
	friend bool operator<(CircuitTableInteger const &a, CircuitTableInteger const &b){return a.v<b.v;}
	friend bool operator<=(CircuitTableInteger const &a, CircuitTableInteger const &b){return (a.v<b.v)||(a.v==b.v);}
	friend bool operator==(CircuitTableInteger const &a, CircuitTableInteger const &b){return a.v==b.v;}
	bool isZero()const{return v.sign()==0;}
	bool isOne()const{return v==1;}
	bool isNonZero()const{return v.sign()!=0;}
	bool isNegative()const{return v.sign()<0;}
	bool isPositive()const{return v.sign()>0;}
	int sign()const{return isNegative()?-1:isPositive();}
	friend int determinantSign1(CircuitTableInteger const &a, CircuitTableInteger const &c, CircuitTableInteger const &b, CircuitTableInteger const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		Integer r=(a.v)*(c.v)-(b.v)*(d.v);
		return r.sign();
	}
	friend int determinantSign3(CircuitTableInteger const &a, CircuitTableInteger const &c, CircuitTableInteger const &b, CircuitTableInteger const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		Integer r=(a.v)*(c.v)-(b.v)*(d.v);
		return r.sign();
	}
	friend int determinantSign2(CircuitTableInteger const &a, CircuitTableInteger const &c, CircuitTableInteger const &b, CircuitTableInteger const &d)//NOTICE MIXED ORDER. MUST WORK BY EXTENDING
	{
		Integer r=(a.v)*(c.v)-(b.v)*(d.v);
		return r.sign();
	}
	std::string toString()const{std::stringstream s;s<<(v.fitsInInt()?std::to_string(v.toInt()):std::string("TOOBIG"));return s.str();}///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int64_t toInt64()const{return v.toInt();}//WHAT SHOULD HAPPEN TO THIS FUNCTION? FOR EXAMPLE IT SEEMS TO RETURN ONLY 32 BITS.
	friend std::ostream &operator<<(std::ostream &f, CircuitTableInteger const &a){f<<a.v;return f;}
	Double extend()const{Double ret;ret.v=v;return ret;}
	CircuitTableInteger &maddWithOverflowChecking(CircuitTableInteger const &a, CircuitTableInteger const&b){v=v+a.v*b.v;return *this;}
	CircuitTableInteger &msubWithOverflowChecking(CircuitTableInteger const &a, CircuitTableInteger const&b){v=v-a.v*b.v;return *this;}
	CircuitTableInteger &subWithOverflowChecking(CircuitTableInteger const &a){Double t=this->extend();t-=a.extend();*this=t.castToSingle();return *this;}
	CircuitTableInteger &addWithOverflowChecking(CircuitTableInteger const &a){Double t=this->extend();t+=a.extend();*this=t.castToSingle();return *this;}
	CircuitTableInteger negatedWithOverflowChecking()const{return (-extend()).castToSingle();}
	friend Double extendedMultiplication(CircuitTableInteger const &a, CircuitTableInteger const &b){return Double(a.v*b.v);}
	friend Double extendedMultiplication(CircuitTableInteger const &a, int32_t b){return Double(a.v*Integer(b));}//to be removed?
	friend CircuitTableInteger min(CircuitTableInteger const &a, CircuitTableInteger const &b){return (b.v<a.v)?b:a;}
	friend CircuitTableInteger max(CircuitTableInteger const &a, CircuitTableInteger const &b){return (a.v<b.v)?b:a;}
	friend CircuitTableInteger negabs(CircuitTableInteger const &a){return min(a,-a);}
	friend CircuitTableInteger abs(CircuitTableInteger const &a){return max(a,a.negatedWithOverflowChecking());}
	friend CircuitTableInteger dotDiv(CircuitTableInteger const &s, CircuitTableInteger const &a, CircuitTableInteger const &t, CircuitTableInteger const &b, CircuitTableInteger::Divisor const &q)
	{
		CircuitTableInteger ret;
		ret.v=(((((((s.v)*(a.v)+(t.v)*(b.v))))))/q.v);
		return ret;
	}
	friend void dotDivAssign(CircuitTableInteger &s, CircuitTableInteger const &a, CircuitTableInteger const &t, CircuitTableInteger const &b, CircuitTableInteger::Divisor const &q)//as above but assigning to first argument. This helps the vectorizer.
	{
		s.v=(((((((s.v)*(a.v)+(t.v)*(b.v))))))/q.v);
	}
	static int MIN(int32_t a, int32_t b)
	{
		return (a<b)?a:b;
	}
	static int MAX(int32_t a, int32_t b)
	{
		return (a>b)?a:b;
	}
	static CircuitTableInteger computeNegativeBound(CircuitTableInteger * __restrict__ Ai, int w)
	{
		CircuitTableInteger M;assert(M.v.isZero());
		CircuitTableInteger m;assert(m.v.isZero());
		for(int j=0;j<w;j++)
		{
			m=min(m,Ai[j]);
			M=max(M,Ai[j]);
		}
		return min(m,-M);
	}
	static CircuitTableInteger dotDivVector(CircuitTableInteger * __restrict__ aa, CircuitTableInteger * __restrict__ bb, CircuitTableInteger s, CircuitTableInteger t, CircuitTableInteger::Divisor denominatorDivisor,int c, CircuitTableInteger boundA, CircuitTableInteger boundB)__attribute__ ((always_inline))//assigning to first argument. The return negative of the return value is an upper bound on the absolute values of the produced entries
	{
		CircuitTableInteger ret;assert(ret.v.isZero());
		for(int i=0;i<c;i++)
		  {
			  CircuitTableInteger temp=(s*aa[i]+t*bb[i]).v/denominatorDivisor.v;
			  ret=max(max(ret,temp),-temp);
			  aa[i]=temp;
		  }
		return ret;
	}
	static CircuitTableInteger scaleVector(CircuitTableInteger * __restrict__ aa, CircuitTableInteger s, CircuitTableInteger::Divisor denominatorDivisor,int c, CircuitTableInteger boundA)__attribute__ ((always_inline))//assigning to first argument. The return negative of the return value is an upper bound on the absolute values of the produced entries
	{
		CircuitTableInteger ret;assert(ret.v.isZero());
		for(int i=0;i<c;i++)
		  {
			  CircuitTableInteger temp=(s*aa[i]).v/denominatorDivisor.v;
			  ret=max(max(ret,temp),-temp);
			  aa[i]=temp;
		  }
		return ret;
	}
	static bool isField()
	{
		return false;
	}
  static CircuitTableInteger gcd(CircuitTableInteger a, CircuitTableInteger b)
	  {
		  // This code appears also for CircuitTableIntPOD and could be a template instead
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
//	  		a-=b;
//	  		if((b-a).isPositive()){swap(a,b);}
	  	}
	  	return a;
	  }
  friend CircuitTableInteger gcd(CircuitTableInteger a, CircuitTableInteger b)
  {
    return CircuitTableInteger::gcd(a,b);
  }
};

static bool hasPod(class CircuitTalbeInteger *c){return true;}//?????????????????????????????????????????

inline CircuitTableInteger CircuitTableInteger::Double::castToSingle()const//casts and checks precission
{
	return CircuitTableInteger(this->v);
}

}


#endif /* SRC_GFANLIB_CIRCUITTABLEINTEGER_H_ */
