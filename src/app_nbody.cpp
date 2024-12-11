#include "vektor.h"
#include "printer.h"
#include "parser.h"
#include "gfanapplication.h"
#include "nbody.h"
#include "field_rationals.h"
#include "buchberger.h"
#include "field_rationalfunctions2.h"
#include <iostream>

class NBodyApplication : public GFanApplication
{
  IntegerOption NOption;
  SimpleOption optionWithMasses;
  SimpleOption optionSymmetric;
  SimpleOption optionAlsoSymmetric;
  SimpleOption optionDziobek;
  IntegerOption optionDeterminants;
  IntegerOption optionDimension;
  SimpleOption optionMM;
  SimpleOption optionSVariables;
  SimpleOption optionLaurent;
  SimpleOption optionVortices;
  SimpleOption optionKaloshin;
  SimpleOption optionAlsoACSpolys;
  SimpleOption optionAlsoBasicAC;
  SimpleOption optionParameters;
  SimpleOption optionCMPar;
public:
  bool includeInDefaultInstallation()
  {
    return false;
  }
  const char *helpText()
  {
    return "This program computes the Albouy-Chenciner equations for the nbody problem.\n";
  }
  NBodyApplication():
    NOption("-N","Specify number of bodies.",3),
    optionWithMasses("--masses","Include mass variables."),
    optionSymmetric("--symmetric","Produce the symmetric equations"),
    optionAlsoSymmetric("--alsosymmetric","Include the symmetric AC equations"),
    optionDziobek("--dziobek","Produce the Dziobek equations"),
    optionDeterminants("--codim","Produce the determinant equations (deprecated)",-1),
    optionDimension("--cayleymenger","Include Cayley-Menger determinants restricting the configuation to the specified dimension."),
    optionMM("--mm","Add in additional mass equation"),
    optionSVariables("-s","Treat the Sij polynomials as variables"),
    optionLaurent("--laurent","Output Laurent polynomials instead of multiplying"),
    optionVortices("--vortices","For the AC equations, produce equations for the vortex case instead of Newtonian case."),
    optionKaloshin("--kaloshin","Produce the polynomials from the Albouy Kaloshin paper"),
	optionAlsoACSpolys("--alsospolys","Include the s-polynomials of the symmetric AC equations with respect to a degree-compatible order."),
  	optionAlsoBasicAC("--alsobasic","Include the basic AC equations."),
	optionParameters("--parameters","Let the masses be parameters of the field."),
	optionCMPar("--cmpar","Just output a parametrization of the relevant CM component.")
  {
	  optionParameters.hide();
	  optionCMPar.hide();
    registerOptions();
  }
  const char *name()
  {
    return "_nbody";
  }
  int offset(int N, int i, int j,bool second=false)
  {
	  int ret=2*N;
	  if(second)ret+=((N-1)*N)/2;
	  while(i>0){ret+=N-i;i--;j--;}
	  return ret+j-1;
  }
  Polynomial varz(PolynomialRing const &r, int i)
  {
    int N=4;
    return r.ithVariable(i);
  }
  Polynomial varw(PolynomialRing const &r, int i)
  {
    int N=4;
    return r.ithVariable(N+i);
  }
  Polynomial varZ(PolynomialRing const &r, int i, int j, int power=1)
  {
    int N=4;
    int sign=(i<j)?1:-1;
    if(j<i)std::swap(i,j);
    return r.polynomialFromField(Q.zHomomorphism(sign))*r.ithVariable(offset(N,i,j),power);
  }
  Polynomial varW(PolynomialRing const &r, int i, int j, int power=1)
  {
    int N=4;
    int sign=(i<j)?1:-1;
    if(j<i)std::swap(i,j);
    return r.polynomialFromField(Q.zHomomorphism(sign))*r.ithVariable(offset(N,i,j,true),power);
  }
  int main()
  {
    AsciiPrinter PP(stdout);
    if(optionCMPar.getValue())
    {
    	int N=NOption.getValue();
    	int D=2;
    	string seperator="Q[";
    	for(int i=1;i<=N;i++)
    		for(int j=0;j<D;j++)
    		{
    			PP<<seperator;
    			PP<<"p"<<i<<j;
    			seperator=",";
    		}
    	for(int i=1;i<=N;i++)
    		for(int j=1;j<i;j++)
    		{
    			PP<<seperator;
    			PP<<"r"<<j<<i;
    			seperator=",";
    		}
    	seperator="]{\n";
    	for(int i=1;i<=N;i++)
    		for(int j=1;j<i;j++)
    		{
    			PP<<seperator;
    			PP<<"r"<<j<<i<<"^2";
    			seperator="";
    			for(int k=0;k<D;k++)
    			{
    				PP<<seperator<<"-p"<<i<<k<<"^2+2*"<<"p"<<i<<k<<"*"<<"p"<<j<<k<<"-p"<<j<<k<<"^2";
    			}
    			seperator=",\n";
    		}
    	PP<<"}\n";
    	return 0;
    }
#define vara varZ
#define varb varW
    if(optionKaloshin.getValue())
      {
    	bool useValuation=true;
    	//a_ij^2=z_ij=z_i-z_j,b_ij^2=w_ij=w_i-w_j
        PolynomialRing r=useValuation?
        			StringParser("Q[z1,z2,z3,z4,w1,w2,w3,w4,a12,a13,a14,a23,a24,a34,b12,b13,b14,b23,b24,b34,t]").parsePolynomialRing()
					:StringParser("Q[z1,z2,z3,z4,w1,w2,w3,w4,a12,a13,a14,a23,a24,a34,b12,b13,b14,b23,b24,b34]").parsePolynomialRing();
//        PolynomialRing r=StringParser("Q[z1,z2,z3,z4,w1,w2,w3,w4,Z12,Z13,Z14,Z23,Z24,Z34,W12,W13,W14,W23,W24,W34]").parsePolynomialRing();

        IntegerVector valuations=StringParser("(2,3,5,7)").parseIntegerVector();
        PolynomialSet g(r);
        int N=4;

        for(int i=0;i<N;i++)
        	for(int j=0;j<i;j++)
        	{
        		g.push_back(vara(r,i,j)*vara(r,i,j)-varz(r,i)+varz(r,j));
        		g.push_back(varb(r,i,j)*varb(r,i,j)-varw(r,i)+varw(r,j));
//        		g.push_back((varz(r,i)-varz(r,j))*(varw(r,i)-varw(r,j))*(varw(r,i)-varw(r,j))*(varw(r,i)-varw(r,j))*varW(r,i,j)*varW(r,i,j)-r.one());
 //       		g.push_back((varw(r,i)-varw(r,j))*(varz(r,i)-varz(r,j))*(varz(r,i)-varz(r,j))*(varz(r,i)-varz(r,j))*varZ(r,i,j)*varZ(r,i,j)-r.one());
        	}





        /*        for(int i=0;i<N;i++)
          for(int j=0;j<i;j++)
            {
              g.push_back((varz(r,i)-varz(r,j))*(varw(r,i)-varw(r,j))*(varw(r,i)-varw(r,j))*(varw(r,i)-varw(r,j))*varW(r,i,j)*varW(r,i,j)-r.one());
              g.push_back((varw(r,i)-varw(r,j))*(varz(r,i)-varz(r,j))*(varz(r,i)-varz(r,j))*(varz(r,i)-varz(r,j))*varZ(r,i,j)*varZ(r,i,j)-r.one());
            }
            */
        for(int i=0;i<N;i++)
          {
            Polynomial A=r.zero()-varz(r,i);
            Polynomial B=r.zero()-varw(r,i);

            for(int j=0;j<N;j++)
            if(i!=j)
              {
            	if(useValuation)
            	{
                    A+=vara(r,i,j,-1)*varb(r,i,j,-3)*r.ithVariable(r.getNumberOfVariables()-1,valuations[j]);//times m
                    B+=varb(r,i,j,-1)*vara(r,i,j,-3)*r.ithVariable(r.getNumberOfVariables()-1,valuations[j]);//times m
            	}
            	else
            	{
            		A+=vara(r,i,j,-1)*varb(r,i,j,-3);//times m
            		B+=varb(r,i,j,-1)*vara(r,i,j,-3);//times m
            	}
              }
            g.push_back(A);
            g.push_back(B);
          }
        {
        	Polynomial Zsum(r),Wsum(r),Isum(r);
            for(int j=0;j<N;j++)
            {
            	Zsum+=r.ithVariable(r.getNumberOfVariables()-1,valuations[j])*varz(r,j);
            	Wsum+=r.ithVariable(r.getNumberOfVariables()-1,valuations[j])*varw(r,j);
            	Isum+=r.ithVariable(r.getNumberOfVariables()-1,valuations[j])*varz(r,j)*varw(r,j);
            }
            g.push_back(Zsum);
            g.push_back(Wsum);
            g.push_back(Isum);
        }
/*
 * This computes the right hand side expression for I, but we cannot add it as an equation, since we only know that it should have a constant value.
 * 	 	{
        	Polynomial I(r);
        	for(int j=0;j<N;j++)
            	for(int i=0;i<j;i++)
            		I+=varZ(r,i,j)*varW(r,i,j)*r.ithVariable(r.getNumberOfVariables()-1,valuations[i]+valuations[j]);
        	g.push_back(I);
        }
        */
        PP<<g.getRing();
        PP<<g;
        return 0;
      }


    PolynomialSet g=AlbouyChencinerEquations(NOption.getValue(),optionWithMasses.getValue(),optionSymmetric.getValue(),optionSVariables.getValue(),!optionLaurent.getValue(),optionVortices.getValue());

    if(optionAlsoSymmetric.getValue())
      {
	auto G=AlbouyChencinerEquations(NOption.getValue(),optionWithMasses.getValue(),true,optionSVariables.getValue(),!optionLaurent.getValue(),optionVortices.getValue());
	for(auto &f:G)g.push_back(f);
      }

    if(optionAlsoBasicAC.getValue())
    {
    	auto G=basicAlbouyChencinerEquations(NOption.getValue(),optionWithMasses.getValue(),true,optionSVariables.getValue(),!optionLaurent.getValue(),optionVortices.getValue());
    	for(auto &f:G)g.push_back(f);
    }

    if(optionAlsoACSpolys.getValue())
      {
    	auto G=AlbouyChencinerEquations(NOption.getValue(),optionWithMasses.getValue(),true,optionSVariables.getValue(),!optionLaurent.getValue(),optionVortices.getValue());
    	AsciiPrinter(stderr)<<G;
		StandardGradedLexicographicTermOrder T;
		G.mark(T);
		for(auto i=G.begin();i!=G.end();i++)
			for(auto j=G.begin();j!=i;j++)
				g.push_back(sPolynomial(*i,*j));
      }

    if(optionSVariables.getValue())
      {
	PolynomialSet s=SEquations(g.getRing(),NOption.getValue(),optionWithMasses.getValue());
	for(PolynomialSet::const_iterator i=s.begin();i!=s.end();i++)g.push_back(*i);
      }

    if(optionDziobek.getValue())
      {
	PolynomialSet dz=DziobekEquations(g.getRing(),NOption.getValue(),optionWithMasses.getValue(),optionSVariables.getValue(),!optionLaurent.getValue());
	for(PolynomialSet::const_iterator i=dz.begin();i!=dz.end();i++)g.push_back(*i);
      }

    if(optionDeterminants.getValue()!=-1)
      {
    	std::cerr << "Warning: option --codim is deprecated. Use option --cayleymenger instead.\n";
    	PolynomialSet de=nbodyDeterminants(g.getRing(),NOption.getValue(),optionWithMasses.getValue(),1+NOption.getValue()-optionDeterminants.getValue());
    	for(PolynomialSet::const_iterator i=de.begin();i!=de.end();i++)g.push_back(*i);
      }
    if(optionDimension.getValue())
      {
    	PolynomialSet de=nbodyDeterminants(g.getRing(),NOption.getValue(),optionWithMasses.getValue(),optionDimension.getValue()+3);
    	for(PolynomialSet::const_iterator i=de.begin();i!=de.end();i++)g.push_back(*i);
      }

    if(optionMM.getValue())g.push_back(massEquation(g.getRing(),NOption.getValue(),optionWithMasses.getValue()));

    if(optionParameters.getValue())
    	g=makeVariablesParameters(makeVariablesParameters(g.getRing(),NOption.getValue()),g);


    PP<<g.getRing();
    PP<<g;

    return 0;
  }
};

static NBodyApplication theApplication;
