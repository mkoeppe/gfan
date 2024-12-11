/*
 * packedmonomial.h
 *
 *  Created on: Feb 7, 2014
 *      Author: anders
 */

#ifndef PACKEDMONOMIAL_H_
#define PACKEDMONOMIAL_H_

#include <bitset>
#include <iostream>
#include "monomial.h"
#include "printer.h"

/*
 * This is an attempt to implement packed monomials in a similar way as Singular does.
 * Methods are designed for speed, not flexibility.
 *
 * Only non-negative values can be stored. That saves a few bits.
 */


//first entry in the packed representation is the weight
//However, to look it up, we use the last entry in varaibleData
class PacMan{
public:
	// For each original coordinate
	struct VariableData
	{
//		int64 mask;
		u_int64_t mask;
		int index;
		int shift;
		VariableData(int64 mask_, int index_, int shift_):
			mask(mask_),
			index(index_),
			shift(shift_)
		{
		}
	};
	vector<VariableData> variableData;

	struct WordData
	{
//		int64 overflowMasks;
		u_int64_t overflowMasks;
		int freeBits;
		WordData():
			overflowMasks(0),
			freeBits(64)
		{
		}
	};
	vector<WordData> wordData;

	int nWords;//number of 64 bit words in packed monomial
	int nExps;//number of varaibles in ring
	int64 overflowFlag;//If overflow happens with these variables, this should be set to non-zero
	/*
	 * The point of setting both addOverflowBit and freeOverflowBit is to save one and
	 * having the overflowbit be part of a different coordinate. Overflow can still be checked,
	 * but only with the trick in dividesNoAdditionalBit.
	 */
	VariableData alloc(int bitsNeeded, bool addOverflowBit, bool freeOverflowBit)
	{
		assert((!freeOverflowBit) || addOverflowBit);
		int j;
		for(j=0;j<wordData.size();j++)
		{
			if(wordData[j].freeBits>=bitsNeeded+addOverflowBit)break;
		}
		if(j==wordData.size())
		{
			wordData.push_back(WordData());
		}
		VariableData ret=VariableData(((bitsNeeded<64)?(((int64)1)<<bitsNeeded):0)-1,j,wordData[j].freeBits-bitsNeeded-addOverflowBit);
		if(addOverflowBit)wordData[j].overflowMasks|=((int64)1)<<(wordData[j].freeBits-1);
		wordData[j].freeBits-=bitsNeeded+addOverflowBit-freeOverflowBit;
		return ret;
	}
	static vector<int> bitsNeeded(IntegerVector const &maxEntries)
		{
			vector<int> ret;
			for(int i=0;i<maxEntries.size();i++)
			{
				int l=0;
				int v=maxEntries[i];
				assert(v>=0);
				while(v>0){v=v>>1;l++;}
				ret.push_back(l);
			}
			return ret;
		}
	bool fits(IntegerVector const &v)
	{
		for(int i=0;i<v.size();i++)if(v[i]!=(v[i]&variableData[i].mask)){/*cerr<<"FADFSAFDA"<<variableData[i].mask;*/return false;}
		return true;
	}
	PacMan()
	{
	}
	PacMan(PolynomialRing const &r, vector<int> const &bounds, int nWeightBits, bool deleteOverflowBit=false, int minimalNumberOfWords=4): //bounds are without overflow bit
		overflowFlag(0),
		wordData(minimalNumberOfWords)
	{
		//		debug<<"Nvar:"<<r.getNumberOfVariables()<<"\n";
		assert(r.getNumberOfVariables()==bounds.size());
		VariableData temp=alloc(nWeightBits,false,false);
		for(int i=0;i<bounds.size();i++)
		{
			variableData.push_back(alloc(bounds[i],true,deleteOverflowBit));
		}
		nWords=wordData.size();
		nExps=variableData.size();
		variableData.push_back(temp);
	}
//	void initializePacked(class PackedMonomial &m)const;
	void initializePacked()const;
	int bytesPerMonomial()const
	{
		return nWords*8;
	}
	void print()
	{
		for(int j=0;j<nExps+1;j++)
		{
			for(int i=0;i<nWords;i++)
				std::cerr<<std::bitset<64>((i==variableData[j].index)?variableData[j].mask<<variableData[j].shift:0)<<" ";
			std::cerr<<"\n";
		}
		for(int i=0;i<nWords;i++)
			std::cerr<<std::bitset<64>(wordData[i].overflowMasks)<<" ";
		std::cerr<<"(overflow mask)\n";
	}
};

// Only use this class with "placement new" when NWORDS=0

// skal eksistere i to varianter
// med fast antal words
// med variabelt antal words, hvor stoerrelse laeses fra manager og plads til data opnaas ved "slicing"
template<int NWORDS>
class PackedMonomial{
public:
	u_int64_t data[NWORDS];
//	int64 data[NWORDS];
public:
	u_int16_t hash()const
	{
		u_int64_t s=0;
		for(int i=0;i<NWORDS;i++)
			s+=data[i];
//		std::cerr<<"S:"<<s<<" "<< (((0xa847c7b39b4369aeull*s)>>48) &0xffffull)<<"\n";
		return ((0xa847c7b39b4369aeull*s)>>48) &0xffffull;
	}
	PackedMonomial()
	{
		for(int i=0;i<NWORDS;i++)
			data[i]=0;
	}
	PackedMonomial(int64 weight, IntegerVector const &exponent, PacMan const &man)
	{
		assert(NWORDS==0||NWORDS>=man.nWords);//Assume that manager does not use too many words
		int nWords=NWORDS?NWORDS:man.nWords;
		for(int i=0;i<nWords;i++)data[i]=0;
		int nExps=man.nExps;
		assert(exponent.size()==nExps);
		data[man.variableData[man.nExps].index]|=(weight&man.variableData[man.nExps].mask)<<man.variableData[man.nExps].shift;
		for(int i=0;i<nExps;i++)
			data[man.variableData[i].index]|=(exponent[i]&man.variableData[i].mask)<<man.variableData[i].shift;
	}
	u_int64_t extractWeight(PacMan const &man)const
	{
		return data[man.variableData.back().index]>>man.variableData.back().shift&man.variableData.back().mask;
	}
	IntegerVector extractExponent(PacMan const &man)const
	{
		int nExps=man.nExps;
		IntegerVector ret(nExps);
		for(int i=0;i<nExps;i++)
		{
//			cerr<<std::bitset<64>(data[man.variableData[i].index])<<" "<<man.variableData[i].mask<<" "<<man.variableData[i].shift<<"\n";
			ret[i]=data[man.variableData[i].index]>>man.variableData[i].shift&man.variableData[i].mask;
		}
			return ret;
	}
	void initializeZero(PacMan const &man)
	{
		int nExps=man.nExps;
		for(int i=0;i<nExps;i++)
			data[i]=0;
	}
	void multiplyBy(PackedMonomial const &b, PacMan const &man)
	{
		int i=0;
		int nWords=NWORDS?NWORDS:man.nWords;
		int64 overflow=0;
		do{
			data[i]+=b.data[i];
			overflow|=data[i]&man.wordData[i].overflowMasks;
		}while(++i<nWords);
		assert(overflow==0);
	}
	bool divides(PackedMonomial const &b, PacMan const &man)const
	{
		int i=0;
		int nWords=NWORDS?NWORDS:man.nWords;
		int64 overflow=0;
		do{
			overflow|=(b.data[i]-data[i])&man.wordData[i].overflowMasks;
		}while(++i<nWords);
		return overflow==0;
	}
	bool dividesNoAdditionalBit(PackedMonomial const &b, PacMan const &man)const
	{
//		assert(man.nWords==NWORDS);
//		assert(man.wordData.size()>=NWORDS);
		int i=0;
		int nWords=NWORDS?NWORDS:man.nWords;
		int64 overflow=0;
		do{
			overflow|=((b.data[i]-data[i])^b.data[i]^data[i])&man.wordData[i].overflowMasks;
		}while(++i<nWords);
		return overflow==0;
	}
	bool dividesLCM(PackedMonomial const &b, PackedMonomial const &c, PacMan const &man)const
	{
		int i=0;
		int nWords=NWORDS?NWORDS:man.nWords;
		int64 overflow=-1;
		do{
			overflow&=((man.wordData[i].overflowMasks+c.data[i]-data[i])|(man.wordData[i].overflowMasks+b.data[i]-data[i]))|~man.wordData[i].overflowMasks;
		}while(++i<nWords);
		return overflow==(int64)-1;
	}
	bool strictlyBiggerOnOneCoordinate(PackedMonomial const &b, PackedMonomial const &c, PacMan const &man)const
	{
		int i=0;
		int nWords=NWORDS?NWORDS:man.nWords;
		int64 overflow=-1;
		do{
			overflow&=(man.wordData[i].overflowMasks+c.data[i]-data[i])|(man.wordData[i].overflowMasks+b.data[i]-data[i])|~man.wordData[i].overflowMasks;
		}while(++i<nWords);
		return overflow!=-1;
	}
// Overflows are assumed not to happen.
	void divideBy(PackedMonomial const &b, PacMan const &man)
	{
		int i=0;
		int nWords=NWORDS?NWORDS:man.nWords;
		do{
			data[i]-=b.data[i];
		}while(++i<nWords);
	}
	void print(PacMan const &man)const
	{
		int nWords=NWORDS?NWORDS:man.nWords;
		for(int i=0;i<nWords;i++)
			std::cerr<<std::bitset<64>(data[i])<<" ";
	}
	void print()const
	{
		for(int i=0;i<NWORDS;i++)
			std::cerr<<std::bitset<64>(data[i])<<" ";
	}
/*	bool compareLess(PackedMonomial const &b,PacMan const &man)
	{
		int nWords=NWORDS?NWORDS:man.nWords;
		int answer=0;
		for(int i=0;i<nWords;i++)//check that this loop is unrolled
		{
			if(!answer){if(b.data[i]<data[i])answer=-1;if(b.data[i]>data[i])answer=-1;}//write this without branches
		}
		return i<0;
	}*/
	bool operator<(PackedMonomial const &b)const
	{
		for(int i=0;i<NWORDS;i++)
		{
			if(data[i]<b.data[i])return true;
			if(data[i]>b.data[i])return false;
		}
		return false;
	}
	bool sumCompare(PackedMonomial const &a2, PackedMonomial const &b1, PackedMonomial const &b2)const
	{
		for(int i=0;i<NWORDS;i++)
		{
			if(data[i]+a2.data[i]<b1.data[i]+b2.data[i])return true;
			if(data[i]+a2.data[i]>b1.data[i]+b2.data[i])return false;
		}
		return false;
	}
	bool isEqualToSum(PackedMonomial const &b1, PackedMonomial const &b2)const
	{
//		cerr<<"isE\n";
		for(int i=0;i<NWORDS;i++)
		{
			if(data[i]!=b1.data[i]+b2.data[i])return false;
		}
		return true;
	}
	PackedMonomial sumNoOverflow(PackedMonomial const &a2)const
	{
		PackedMonomial<NWORDS> r;
		for(int i=0;i<NWORDS;i++)
			r.data[i]=data[i]+a2.data[i];
		return r;
	}
	bool operator==(PackedMonomial const &a2)const
		{
		for(int i=0;i<NWORDS;i++)
			if(data[i]!=a2.data[i])return false;
			return true;
		}
};

// For now we just assume that coefficients are one
template<int NWORDS>
class PackedPolynomial
{
	vector<PackedMonomial<NWORDS> > p;

	PackedPolynomial(IntegerVector const &weight, Polynomial const &p_, PacMan const &man)
	{
		for(TermMap::const_iterator i=p_.terms.begin();i!=p_.terms.end();i++)
		{
			p.push_back(PackedMonomial<NWORDS>(dotLong(weight,i->first.exponent),i->first,man));//dot product may overflow
		}
		debug<<"SORTDISABLED";
//		sort(p.begin(),p.end(),[](int a, int b)
//			    {
//					return 1;
//			    });
	}
	void print(PacMan const &man)
	{
//		for(vector<PackedMonomial<NWORDS> >::const_iterator i=p.begin();i!=p.end();i++)i->print();
		for(int i=0;i<p.size();i++)p[i].print();
		debug<<"\n";
	}
};
void packedTest();

#endif /* PACKEDMONOMIAL_H_ */
