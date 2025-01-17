/*
 * app_anders.cpp
 *
 *  Created on: Aug 11, 2023
 *      Author: au27818
 */

#include  <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "buchberger.h"
#include "wallideal.h"
#include "termorder.h"
#include "gfanapplication.h"
#include "tropical2.h"
#include "tropical.h"
#include "saturation.h"
#include "log.h"

using namespace std;

class AndersApplication : public GFanApplication
{
	IntegerOption start;
	IntegerOption endSkip;
	SimpleOption ignoreNegative;
	SimpleOption outputIdealToBeSaturated;
	SimpleOption homogenize;
public:
	bool includeInDefaultInstallation()
	{
		return false;
	}
  const char *helpText()
  {
    return "This program takes: Ring, ideal, generators, list of vectors. For each vector it compute the sum of the initial ideal and the ideal generated by the initial foms. Then it does a saturation.\n";
  }
  AndersApplication():
	  start("--start","Specify number of vectors to skip at start.",0),
	  endSkip("--endskip","Specify number of vectors to skip at end.",0),
	  ignoreNegative("--ignorenegative","Ignore a vector if its coordinate sum is negative."),
	  outputIdealToBeSaturated("--showideal","Prints the ideal that needs saturation."),
	  homogenize("-h","Homogenise the generators before saturating and use faster saturation.")
  {
    registerOptions();
  }

  const char *name()
  {
    return "_anders";
  }

  int main()
  {
    FileParser P(Stdin);

    PolynomialSet g=P.parsePolynomialSetWithRing();
    PolynomialSet h=P.parsePolynomialSet(g.getRing());
    IntegerVectorList wList=P.parseIntegerVectorList();

    vector<IntegerVector> list;
    for(auto &w:wList)list.push_back(w);


    for(int i=start.getValue();i<list.size()-endSkip.getValue();i++)
    	if(ignoreNegative.getValue()&&list[i].sum()<0)
    	{
        	debug<<"Ignoring vector "<<i<<list[i]<<"\n";
    	}
    	else
    {
    	auto w=list[i];
    	debug<<"Processing vector "<<i<<w<<"\n";
    	WeightReverseLexicographicTermOrder T(w);
    	PolynomialSet g2=g;
    	debug<<"Computing initial ideal.\n";
    	buchberger(&g2,T);
    	debug<<"Done computing initial ideal.\n";
    	g2=initialForms(g2,w);
    	h.mark(T);
    	auto h2=initialForms(h,w);
        for(auto &p:h2)p.saturate();
    	debug<<"Computing saturation of initial ideal.\n";
//    	g2=saturatedIdeal(g2);
    	g2=nonHomogeneousSaturation(g2);
    	debug<<"Done computing saturation of initial ideal.\n";
    	g2.unionSet(h2);

    	auto r2=g2.getRing().withVariablesAppended("H");

    	if(homogenize.getValue())
    		g2=g2.homogenization(r2,0);

    	if(outputIdealToBeSaturated.getValue())
    		pout<<g2.getRing()<<g2<<"\n";

    	debug<<"Computing saturation of sum.\n";

    	if(homogenize.getValue())
    		g2=saturatedIdeal(g2);
    	else
    		g2=nonHomogeneousSaturation(g2);

    	debug<<"Done computing saturation of sum.\n";
    	pout<<g2;
    }
    return 0;
  }
};

static AndersApplication theApplication;



