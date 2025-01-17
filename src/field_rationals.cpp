#include "field_rationals.h"

#include <memory>
#include <assert.h>
#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>
#include <iostream>

#include "lp_cdd.h"
#include "printer.h"

#include "timer.h"
#include "log.h"

#include "field_zmodpz.h"
#include "rational.h"
//static Timer rationalTimer("Rational",1);

#include <thread>
#include <mutex>

int FieldElementRationalsLiving;

static mutex mtx;
class A{
public:
A()
{
	//mtx.lock();
}
~A()
{
	//mtx.unlock();
}
};


class FieldElementRational : public FieldElementImplementation
{
 public:
  mpq_t value;
  FieldElementRational(FieldImplementation &a):FieldElementImplementation(a)
    {
	  A aaaa;
      FieldElementRationalsLiving++;
      mpq_init(value);
    }
  FieldElementRational(FieldImplementation &a,int n_):FieldElementImplementation(a)
    {
	  A aaaa;
      FieldElementRationalsLiving++;
      mpz_init_set_si(mpq_numref(value), n_);
      mpz_init_set_ui(mpq_denref(value), 1);
    }
  FieldElementRational(FieldImplementation &a, mpq_t *n_):FieldElementImplementation(a)
  {
	  A aaaa;
    FieldElementRationalsLiving++;
    mpq_init(value);
    mpq_set(value,*n_);
  }
  virtual ~FieldElementRational()
    {
	  A aaaa;
      FieldElementRationalsLiving--;
      mpq_clear(value);
    }
  FieldElementRational& operator=(const FieldElementRational& a)
    {
	  A aaaa;
      assert(0);
      const FieldElementRational *A=(const FieldElementRational*)&a;
      if (this != A) {
        mpq_clear(value);
        mpz_init_set(mpq_numref(value), mpq_numref(a.value));
        mpz_init_set(mpq_denref(value), mpq_denref(a.value));
      }
      return *this;
    }

  mpq_t const *getGmpRationalTemporaryPointer()const
    {
     return &value;
    }
  bool isInteger()const
  {
	  ::A aaaa;
    return mpz_cmp_si(mpq_denref(value),1)==0;
  }
  FieldElementImplementation *multiplierGivingInteger()const
  {
	  ::A aaaa;
	  FieldElementRational *r= new FieldElementRational(*getField());

      mpz_set(mpq_numref(r->value), mpq_denref(value));

	  return r;
  }
  void operator*=(const FieldElementImplementation &a)
    {
      const FieldElementRational *A=(const FieldElementRational*)&a;
      assert(A);
	  ::A aaaa;

      //      TimerScope ts(&rationalTimer);
      mpq_mul(value,value,A->value);
    }
  void operator+=(const FieldElementImplementation &a)
    {
      const FieldElementRational *A=(const FieldElementRational*)&a;
      assert(A);
	  ::A aaaa;

      //      TimerScope ts(&rationalTimer);
      mpq_add(value,value,A->value);

//    mpq_canonicalize(value);
    }
  void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)
    {
      const FieldElementRational *A=(const FieldElementRational*)&a;
      assert(A);
      const FieldElementRational *B=(const FieldElementRational*)&b;
      assert(B);

	  ::A aaaa;
      //      TimerScope ts(&rationalTimer);
      mpq_t temp;
      mpq_init(temp);
      mpq_mul(temp,A->value,B->value);
      mpq_add(value,value,temp);
      mpq_clear(temp);
    }


  FieldElementRational *one() const;
  bool isZero() const
    {
      return mpq_sgn(value)==0;
    }

  FieldElementRational *sum(const FieldElementImplementation &b)const
    {
      //      TimerScope ts(&rationalTimer);
      const FieldElementRational *B=(const FieldElementRational*)&b;
      FieldElementRational *r= new FieldElementRational(*getField());
      // fprintf(Stderr,"NEW\n");
	  A aaaa;
      mpq_add(r->value,value,B->value);
      return r;
    }
  FieldElementRational *difference(const FieldElementImplementation &b)const
    {
      //      TimerScope ts(&rationalTimer);
      const FieldElementRational *B=(const FieldElementRational*)&b;
      FieldElementRational *r= new FieldElementRational(*getField());
      // fprintf(Stderr,"NEW\n");
	  A aaaa;
      mpq_sub(r->value,value,B->value);
      return r;
    }
  FieldElementRational *negation()const
    {
      FieldElementRational *r= new FieldElementRational(*getField());
	  A aaaa;
      // fprintf(Stderr,"NEW\n");
      mpq_neg(r->value,value);

      return r;
    }
  FieldElementImplementation *inverse()const
  {
    FieldElementRational *r= new FieldElementRational(*getField());
    // fprintf(Stderr,"NEW\n");

    if(isZero())
      {
	AsciiPrinter P(Stderr);
	P.printString("Error inverting FieldElement: ");
		//P.printFieldElement(*this);
	P.printString(toString());
	P.printString("\n");
	assert(0);
      }

	  A aaaa;
    mpq_inv(r->value,value);
    return r;
  }

  int sign()const
  {
    return mpq_sgn(value);
  }

  static int val(mpz_t m, int p)
  {
    int ret=0;
    if(p==2)
      {
        int p2Val=mpz_scan1(m,0);
//        assert(0);
//        assert(p2Val==ret);
        return p2Val;
      }
    while(mpz_divisible_ui_p(m,p))
      {
        mpz_divexact_ui(m,m,p);
        ret++;
      }

    return ret;
  }
  int pAdicValuation(int p)const
  {
    mpq_t temp;
    mpq_init(temp);
    mpq_set(temp,value);
    int ret=val(mpq_numref(temp),p)-val(mpq_denref(temp),p);
    mpq_clear(temp);
    return ret;
  }

  FieldElement pAdicRemainder(Field const &ZModPZ)const
  {
#if USESHAREDPTR
	  int p=(dynamic_cast<const FieldZModPZImplementation *>(ZModPZ.implementingObject.get()))->getCharacteristic();
#else
	  int p=(dynamic_cast<const FieldZModPZImplementation *>(ZModPZ.implementingObject))->getCharacteristic();
#endif
    mpz_t temp;
    mpz_init(temp);
    mpz_set(temp,mpq_numref(value));
    int v=pAdicValuation(p);
    while(v>0){mpz_divexact_ui(temp,temp,p);v--;}
    FieldElement A=ZModPZ.zHomomorphism(mpz_fdiv_ui(temp,p));

    mpz_set(temp,mpq_denref(value));
    while(v<0){mpz_divexact_ui(temp,temp,p);v++;}
    FieldElement B=ZModPZ.zHomomorphism(mpz_fdiv_ui(temp,p));
    mpz_clear(temp);
    //    FieldElement A=ZModPZ.zHomomorphism(mpz_fdiv_ui(mpq_numref(value),p));
//    FieldElement B=ZModPZ.zHomomorphism(mpz_fdiv_ui(mpq_denref(value),p));
    assert(!B.isZero());
    return A*(B.inverse());
  }

  FieldElement modularRepresentative(Field const &ZModPZ)const
  {
	  assert(isInteger());
#if USESHAREDPTR
	  int p=(dynamic_cast<const FieldZModPZImplementation *>(ZModPZ.implementingObject.get()))->getCharacteristic();
#else
	  int p=(dynamic_cast<const FieldZModPZImplementation *>(ZModPZ.implementingObject))->getCharacteristic();
#endif
	  mpz_t temp;
	  mpz_init(temp);
	  mpz_set(temp,mpq_numref(value));
	  mpz_mod_ui(temp,temp,p);
	  FieldElement A=ZModPZ.zHomomorphism(mpz_fdiv_ui(temp,p));
	  mpz_clear(temp);
	  return A;
  }

  static string LaTeXTranslator(const string &s)
  {
    int startIndex=0;
    string sign;
    if(s[0]=='-')
      {
	sign=string("-");
	startIndex=1;
      }
    int slashIndex=-1;
    for(int i=startIndex;i<s.length();i++)if(s[i]=='/')slashIndex=i;
    if(slashIndex==-1)
      return string(s);

    return sign+string("{").append(s,startIndex,slashIndex-startIndex)+string("\\over ").append(s,slashIndex+1,s.length()-slashIndex-1)+string("}");
  }

  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false, bool latexMode=false) const
  {
    FieldElementRational *tempOne=one();
    FieldElementRational *temp=difference(*tempOne);
    bool isOne=temp->isZero();
    //fprintf(Stderr,"DELETE\n");
    delete temp;
    temp=sum(*tempOne);
    bool isMinusOne=temp->isZero();
    //fprintf(Stderr,"DELETE\n");
    delete temp;
    //fprintf(Stderr,"DELETE\n");
    delete tempOne;

    if(!writeIfOne && isOne)
      {
	if(alwaysWriteSign)return std::string("+");
	return std::string("");
      }
    if(!writeIfOne && isMinusOne)
      return std::string("-");
    static char s[1290*1000];
    //    mpq_get_str(s,10,value); //// CHECK BUFFER SIZE!!!!
    // Changed to make code gmp 3.1.1 compatible
    mpz_get_str(s,10,mpq_numref(value)); //CHECK BUFFER SIZE!!!!
    string S(s);
    if(mpz_cmp_ui(mpq_denref(value),1)!=0)
      {
	mpz_get_str(s,10,mpq_denref(value)); //CHECK BUFFER SIZE!!!!
	S=S+string("/")+string(s);
      }

    if(latexMode)S=LaTeXTranslator(S);

    if(alwaysWriteSign && mpq_sgn(value)!=-1)
      return std::string("+")+S;
    return S;
  }

  FieldElementRational *copy()const
  {
    FieldElementRational *r= new FieldElementRational(*getField());
    //      fprintf(Stderr,"NEW\n");

	  A aaaa;
    mpq_clear(r->value);
    mpz_init_set(mpq_numref(r->value), mpq_numref(value));
    mpz_init_set(mpq_denref(r->value), mpq_denref(value));

    return r;
  }
};


bool FieldRationalsImplementation::isRationals()const
{
  return true;
}

int FieldRationalsImplementation::getCharacteristic()const
{
	return 0;
}

FieldRationalsImplementation::FieldRationalsImplementation()
{
}

std::string FieldRationalsImplementation::toString()const
{
  return std::string("Q");
}

FieldElementImplementation *FieldRationalsImplementation::zHomomorphismImplementation(int n)
{
/*  if(n==0) //NEED TO CHANGE THIS BACK SOMEHOW
    {
      static FieldElementImplementation *p;
      if(p==0)p=new FieldElementRational(*this,0);
      p->refCount++;
      return p;
    }
  else
    if(n==1)
      {
	static FieldElementImplementation *p;
	if(p==0)p=new FieldElementRational(*this,1);
	p->refCount++;
	return p;
      }
    else
      if(n==-1)
	{
	  static FieldElementImplementation *p;
	  if(p==0)p=new FieldElementRational(*this,-1);
	  p->refCount++;
	  return p;
	}*/
  FieldElementImplementation *ret=new FieldElementRational(*this,n);
  //      fprintf(Stderr,"NEW\n");
  //  ret->refCount++;
  return ret;
}

FieldElement FieldRationalsImplementation::zHomomorphism(int n)
{
  //      fprintf(Stderr,"NEW\n");
  //  return FieldElement(new FieldElementRational(*this,n));
  return FieldElement(zHomomorphismImplementation(n));
}

const char *FieldRationalsImplementation::name()
{
  return "GmpRationals";
}

/*FieldRationals::FieldRationals():
  Field(new FieldRationalsImplementation())
{
  /*  fprintf(Stderr,"Adding field rationals\n");
  next=list;
  list=this;
  */
/*
  log2 fprintf(Stderr,"Initializing field Rationals\n");
}
*/
				//FieldRationals Q;
Field Q(new FieldRationalsImplementation());

FieldElementRational *FieldElementRational::one() const
{
  //      fprintf(Stderr,"NEW\n");
  return new FieldElementRational(*getField(),1);
}


void printMpq(mpq_t value)
{
  char s[1000];
  char t[1000];
  mpz_get_str(s,10,mpq_numref(value)); //CHECK BUFFER SIZE!!!!
  mpz_get_str(t,10,mpq_denref(value)); //CHECK BUFFER SIZE!!!!
  fprintf(stderr,"%s/%s ",s,t);
}


IntegerVector primitiveVectorOld(vector<FieldElement> const &v)
{
  int n=v.size();

  mpq_t *point = new mpq_t [n];
  for(int i=0;i<n;i++)mpq_init(point[i]);

  for(int i=0;i<n;i++)
    {
      mpq_set(point[i],*(v[i].getGmpRationalTemporaryPointer()));
    }

  scaleToIntegerVector(point,n);

  IntegerVector ret=arrayToIntegerVector(point,n);

  for(int i=0;i<n;i++)mpq_clear(point[i]);

  delete [] point;

  return ret;
  //  return IntegerVector(n);
}


IntegerVector primitiveVector(vector<FieldElement> const &v)
{
  int n=v.size();
  IntegerVector ret(n);

  mpz_t lcm;
  mpz_t gcd;
  mpz_init_set_ui(lcm, 1);
  mpz_init_set_ui(gcd, 0);

  for(int j=0;j<n;j++)
    {
      mpq_t const *a=v[j].getGmpRationalTemporaryPointer();
      if(mpz_cmp_si(mpq_denref(*a),1)!=0)
        mpz_lcm(lcm,lcm,mpq_denref(*a));
      if(mpz_sgn(mpq_numref(*a))!=0)
        mpz_gcd(gcd,gcd,mpq_numref(*a));
    }
  if(mpz_sgn(gcd)!=0)//v is non-zero
    {
      if((mpz_cmp_si(lcm,1)==0)&&(mpz_cmp_si(gcd,1)==0)) //gcd=lcm=1
        {
          for(int i=0;i<n;i++)
            {
              mpq_t const *a=v[i].getGmpRationalTemporaryPointer();
              if((!mpz_fits_sint_p(mpq_numref(*a))))
                {
                  fprintf(stderr,"INTEGER OVERFLOW\n");
                  assert(0);
                }
              ret[i]=mpz_get_si(mpq_numref(*a));
            }
        }
      else
      {
    	  mpz_t tempA;
          mpz_t tempB;
          mpz_init(tempA);
          mpz_init(tempB);
          for(int i=0;i<n;i++)
            {
              mpq_t const *a=v[i].getGmpRationalTemporaryPointer();
              mpz_set(tempA,mpq_denref(*a));
              mpz_set(tempB,mpq_numref(*a));
              mpz_mul(tempA,gcd,tempA);
              mpz_mul(tempB,lcm,tempB);
              mpz_divexact(tempA,tempB,tempA);

              if((!mpz_fits_sint_p(tempA)))
                {
                  fprintf(stderr,"INTEGER OVERFLOW\n");
                  debug<<"Converting (";
                  for(int i=0;i<v.size();i++)
                    {
                      if(i!=0)debug<<",";
                      debug<<v[i];
                    }
                    debug<<") to primitive vector\n";
                  assert(0);
                }
              ret[i]=mpz_get_si(tempA);
            }
          mpz_clear(tempB);
          mpz_clear(tempA);
        }
    }
  mpz_clear(gcd);
  mpz_clear(lcm);

  return ret;
}


mpq_t *fieldElementToGmp(FieldElement const &c)
{
#if USESHAREDPTR
  FieldElementImplementation *i=c.implementingObject.get();
#else
  FieldElementImplementation *i=c.implementingObject;
#endif
  FieldElementRational *i2=dynamic_cast<FieldElementRational*>(i);

  assert(i2);

  return &(i2->value);
}


double fieldElementToFloatingPoint(FieldElement const&c)
{
	return mpq_get_d(*fieldElementToGmp(c));
}


FieldElement fieldElementFromGmp(mpq_t *c)
{
  FieldElement ret=FieldElement(new FieldElementRational(*(Q.implementingObject),c));
  return ret;
}


FieldElement gcd(FieldElement const &a, FieldElement const &b, FieldElement &s, FieldElement &t)
{
  mpz_t G;
  mpz_t S;
  mpz_t T;

  mpz_init(G);
  mpz_init(S);
  mpz_init(T);
  mpz_gcdext(G, S, T, mpq_numref(*fieldElementToGmp(a)), mpq_numref(*fieldElementToGmp(b)));

  mpq_t Grat;
  mpq_init(Grat);
  mpq_set_z(Grat, G);

  mpq_t Srat;
  mpq_init(Srat);
  mpq_set_z(Srat, S);

  mpq_t Trat;
  mpq_init(Trat);
  mpq_set_z(Trat, T);

  FieldElement ret=fieldElementFromGmp(&Grat);
  s=fieldElementFromGmp(&Srat);
  t=fieldElementFromGmp(&Trat);

  mpq_clear(Trat);
  mpq_clear(Srat);
  mpq_clear(Grat);
  mpz_clear(T);
  mpz_clear(S);
  mpz_clear(G);

  return ret;
}

IntegerVector toIntegerVector(vector<FieldElement> const &v)
{
  int n=v.size();

  mpq_t *point = new mpq_t [n];
  for(int i=0;i<n;i++)mpq_init(point[i]);
  for(int i=0;i<n;i++)
    {
      mpq_set(point[i],*(v[i].getGmpRationalTemporaryPointer()));
    }

  IntegerVector ret=arrayToIntegerVector(point,n);

  for(int i=0;i<n;i++)mpq_clear(point[i]);

  delete [] point;

  return ret;
}


FieldElement fieldElementFromGmpZ(const mpz_t *c)
{
  mpq_t temp;
  mpq_init(temp);
  mpz_set(mpq_numref(temp),*c);
  mpz_set_si(mpq_denref(temp),1);
  mpq_canonicalize(temp);
  FieldElement ret=fieldElementFromGmp(&temp);
  mpq_clear(temp);
  return ret;
}

FieldElement integerDivision(FieldElement const &a, FieldElement const &b, FieldElement *remainder)
{
  mpz_t q,r,B;
  mpz_init(q);
  mpz_init(r);
  mpz_init(B);
  assert(a.isInteger());
  assert(b.isInteger());
  assert(!b.isZero());
  mpz_set(B,mpq_numref(*b.getGmpRationalTemporaryPointer()));
  int s=mpz_sgn(B);
  mpz_abs(B,B);
  mpz_fdiv_qr(q,r,mpq_numref(*a.getGmpRationalTemporaryPointer()),B);
  mpz_mul_si(q,q,s);
//    debug<<fieldElementFromGmpZ(mpq_numref(*a.getGmpRationalTemporaryPointer()));
//  debug<<" "<<fieldElementFromGmpZ(&B)<<"\n\n";

  if(remainder)*remainder=fieldElementFromGmpZ(&r);
  FieldElement ret=fieldElementFromGmpZ(&q);

//  debug<<"r="<<fieldElementFromGmpZ(&r)<<"   q="<<ret<<"  a="<<a<<"  b="<<b<<"\n";
  assert((ret*b-(a-fieldElementFromGmpZ(&r))).isZero());

  mpz_clear(B);
  mpz_clear(r);
  mpz_clear(q);

  return ret;
}

// Fix this. Should not depend on linalg.h.
#include "linalg.h"

int toInteger(FieldElement const &a)
{
//	cerr<<"toInteger1\n";
  FieldVector A(a.getField(),1);
//	cerr<<"toInteger2\n";
  A[0]=a;
//	cerr<<"toInteger3\n";
  IntegerVector temp=fieldVectorToIntegerVector(A);
//  cerr<<"AFADSAFSDA";
  return temp[0];
}




//------------------------------------




class FieldElementRational2 : public FieldElementImplementation
{
 public:
  Rational value;
  FieldElementRational2(FieldImplementation &a):FieldElementImplementation(a)
    {
      FieldElementRationalsLiving++;
    }
  FieldElementRational2(FieldImplementation &a,int n_):FieldElementImplementation(a),
		  value(n_)
    {
      FieldElementRationalsLiving++;
    }
  FieldElementRational2(FieldImplementation &a, mpq_t *n_):FieldElementImplementation(a)
  {
	  assert(0);//not yet implemented
	  FieldElementRationalsLiving++;
//    mpq_init(value);






	  //    mpq_set(value,*n_);
  }




  virtual ~FieldElementRational2()
    {
      FieldElementRationalsLiving--;
    }
  FieldElementRational2& operator=(const FieldElementRational2& a)
    {
      assert(0);
      return *this;
    }

  mpq_t const *getGmpRationalTemporaryPointer()const
    {
	  assert(0);//not yet implemented
	return 0;
//      return &value;
    }
  bool isInteger()const
  {
	  assert(0);
	  return 0;
//    return mpz_cmp_si(mpq_denref(value),1)==0;
  }
  void operator*=(const FieldElementImplementation &a)
    {
      const FieldElementRational2 *A=(const FieldElementRational2*)&a;
      assert(A);

//      cerr<<toString()<<"*"<<a.toString()<<"=";

      value*=A->value;
//      cerr<<toString()<<"\n";
      //      TimerScope ts(&rationalTimer);
//      mpq_mul(value,value,A->value);
    }
  void operator+=(const FieldElementImplementation &a)
    {
      const FieldElementRational2 *A=(const FieldElementRational2*)&a;
      assert(A);
assert(0);
      cerr<<toString()<<"+"<<a.toString()<<"=";
      value+=A->value;
      cerr<<toString()<<"\n";
      //      TimerScope ts(&rationalTimer);
//      mpq_add(value,value,A->value);

//    mpq_canonicalize(value);
    }
  void madd(const FieldElementImplementation &a,const FieldElementImplementation &b)
    {
      const FieldElementRational2 *A=(const FieldElementRational2*)&a;
      assert(A);
      const FieldElementRational2 *B=(const FieldElementRational2*)&b;
      assert(B);
assert(0);
      cerr<<toString()<<"+"<<a.toString()<<"*"<<b.toString()<<"=";
      value+=A->value*B->value;
      cerr<<toString()<<"\n";
    }


  FieldElementRational2 *one() const;
  bool isZero() const
    {
	  return value.isZero();
    }



  FieldElementRational2 *sum(const FieldElementImplementation &b)const
    {

	  //      TimerScope ts(&rationalTimer);
      const FieldElementRational2 *B=(const FieldElementRational2*)&b;
      FieldElementRational2 *r= new FieldElementRational2(*getField());
      // fprintf(Stderr,"NEW\n");
//      mpq_add(r->value,value,B->value);
//      cerr<<toString()<<"+"<<b.toString()<<"=";
      r->value=value+B->value;
//      cerr<<toString()<<"\n";
//  	P.printString(this->value.toString());
      return r;
    }
  FieldElementRational2 *difference(const FieldElementImplementation &b)const
    {
      //      TimerScope ts(&rationalTimer);
      const FieldElementRational2 *B=(const FieldElementRational2*)&b;
      FieldElementRational2 *r= new FieldElementRational2(*getField());
      // fprintf(Stderr,"NEW\n");
      r->value=value-B->value;
//      mpq_sub(r->value,value,B->value);
      return r;
    }
  FieldElementRational2 *negation()const
    {
      FieldElementRational2 *r= new FieldElementRational2(*getField());
      // fprintf(Stderr,"NEW\n");
//      mpq_neg(r->value,value);
//      cerr<<"NEW"<<value <<"\n";
      r->value=-value;
//cerr<<"NEW"<<r->value <<"\n";
      return r;







    }
  FieldElementImplementation *inverse()const
  {
    FieldElementRational2 *r= new FieldElementRational2(*getField());
    // fprintf(Stderr,"NEW\n");

    if(isZero())
      {
	AsciiPrinter P(Stderr);
	P.printString("Error inverting FieldElement: ");
	//	P.printFieldElement(*this);
	P.printString(this->value.toString());
	P.printString("\n");
	assert(0);
      }

    r->value=Rational(1)/value;
    //    mpq_inv(r->value,value);
    return r;
  }

  int sign()const
  {
	  return value.sign();
//    return mpq_sgn(value);
  }

  static int val(mpz_t m, int p)
  {
	assert(0);
/*    int ret=0;
    if(p==2)
      {
        int p2Val=mpz_scan1(m,0);
        return p2Val;
      }
    while(mpz_divisible_ui_p(m,p))
      {
        mpz_divexact_ui(m,m,p);
        ret++;
      }

    return ret;
 */
	return 0;
  }
  int pAdicValuation(int p)const
  {
		assert(0);
/*    mpq_t temp;
    mpq_init(temp);
    mpq_set(temp,value);
    int ret=val(mpq_numref(temp),p)-val(mpq_denref(temp),p);
    mpq_clear(temp);
    return ret;*/
	return 0;
  }

  FieldElement pAdicRemainder(Field const &ZModPZ)const
  {
	  assert(0);
/*
	 int p=(dynamic_cast<const FieldZModPZImplementation *>(ZModPZ.implementingObject))->getP();
    mpz_t temp;
    mpz_init(temp);
    mpz_set(temp,mpq_numref(value));
    int v=pAdicValuation(p);
    while(v>0){mpz_divexact_ui(temp,temp,p);v--;}
    FieldElement A=ZModPZ.zHomomorphism(mpz_fdiv_ui(temp,p));

    mpz_set(temp,mpq_denref(value));
    while(v<0){mpz_divexact_ui(temp,temp,p);v++;}
    FieldElement B=ZModPZ.zHomomorphism(mpz_fdiv_ui(temp,p));
    mpz_clear(temp);
    assert(!B.isZero());
    return A*(B.inverse());
*/
//	  return *this;
	  return ZModPZ.zHomomorphism(0);
  }

  static string LaTeXTranslator(const string &s)
  {
    int startIndex=0;
    string sign;
    if(s[0]=='-')
      {
	sign=string("-");
	startIndex=1;
      }
    int slashIndex=-1;
    for(int i=startIndex;i<s.length();i++)if(s[i]=='/')slashIndex=i;
    if(slashIndex==-1)
      return string(s);

    return sign+string("{").append(s,startIndex,slashIndex-startIndex)+string("\\over ").append(s,slashIndex+1,s.length()-slashIndex-1)+string("}");
  }

  std::string toString(bool writeIfOne=true, bool alwaysWriteSign=false, bool latexMode=false) const
  {

	  if(latexMode) return LaTeXTranslator(value.toString(writeIfOne,alwaysWriteSign));//does translator work with + signs?

//	    stringstream s;
//	  s<<value;

	  return value.toString(writeIfOne,alwaysWriteSign);
//	  assert(0);
/*    FieldElementRational *tempOne=one();
    FieldElementRational *temp=difference(*tempOne);
    bool isOne=temp->isZero();
    //fprintf(Stderr,"DELETE\n");
    delete temp;
    temp=sum(*tempOne);
    bool isMinusOne=temp->isZero();
    //fprintf(Stderr,"DELETE\n");
    delete temp;
    //fprintf(Stderr,"DELETE\n");
    delete tempOne;

    if(!writeIfOne && isOne)
      {
	if(alwaysWriteSign)return std::string("+");
	return std::string("");
      }
    if(!writeIfOne && isMinusOne)
      return std::string("-");
    static char s[1290*1000];
    //    mpq_get_str(s,10,value); //// CHECK BUFFER SIZE!!!!
    // Changed to make code gmp 3.1.1 compatible
    mpz_get_str(s,10,mpq_numref(value)); //CHECK BUFFER SIZE!!!!
    string S(s);
    if(mpz_cmp_ui(mpq_denref(value),1)!=0)
      {
	mpz_get_str(s,10,mpq_denref(value)); //CHECK BUFFER SIZE!!!!
	S=S+string("/")+string(s);
      }

    if(latexMode)S=LaTeXTranslator(S);

    if(alwaysWriteSign && mpq_sgn(value)!=-1)
      return std::string("+")+S;
    return S;

*/
//  return "";
  }

  FieldElementRational2 *copy()const
  {
    FieldElementRational2 *r= new FieldElementRational2(*getField());
    //      fprintf(Stderr,"NEW\n");

/*    mpq_clear(r->value);
    mpz_init_set(mpq_numref(r->value), mpq_numref(value));
    mpz_init_set(mpq_denref(r->value), mpq_denref(value));
*/
    r->value=value;
    return r;
  }
};


FieldElementImplementation *FieldRationals2Implementation::zHomomorphismImplementation(int n)
{
//	cerr<<n<<" ";
/*  if(n==0)
    {
      static FieldElementImplementation *p;
      if(p==0)p=new FieldElementRational2(*this,0);
      p->refCount++;
      return p;
    }
  else
    if(n==1)
      {
	static FieldElementImplementation *p;
	if(p==0)p=new FieldElementRational2(*this,1);
	p->refCount++;
	return p;
      }
    else
      if(n==-1)
	{
	  static FieldElementImplementation *p;
	  if(p==0)p=new FieldElementRational2(*this,-1);
	  p->refCount++;
	  return p;
	}*/
  FieldElementImplementation *ret=new FieldElementRational2(*this,n);
//        fprintf(Stderr,"NEW\n");
  //  ret->refCount++;

//  cerr<<ret->toString();
  return ret;
}

FieldElement FieldRationals2Implementation::zHomomorphism(int n)
{
	return FieldElement(zHomomorphismImplementation(n));
}

const char *FieldRationals2Implementation::name()
{
	return "Rationals2";
}

std::string FieldRationals2Implementation::toString()const
{
	return std::string("QQ");
}

bool FieldRationals2Implementation::isRationals()const
{
	return false;//TO BE CHANGED?
}


int FieldRationals2Implementation::getCharacteristic()const
{
	return 0;
}

FieldRationals2Implementation::FieldRationals2Implementation()
{
}

Field QQ(new FieldRationals2Implementation());

FieldElementRational2 *FieldElementRational2::one() const
{
  //      fprintf(Stderr,"NEW\n");
  return new FieldElementRational2(*getField(),1);
}
