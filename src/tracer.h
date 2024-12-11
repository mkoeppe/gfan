/*
 * tracer.h
 *
 *  Created on: 6 Aug 2018
 *      Author: anders
 */

/*
 * The tracer is useful for comparing intermediate results of two slightly different runs of the code.
 * Commands like "trace('A',to_string(a));" can be inserted at various places.
 * When the Tracer object is initialised with 0 it records and when initialised with 1 it compares with a previous run.
 * If any mismatch is discovered an assertion fails.
 */

#ifndef SRC_TRACER_H_
#define SRC_TRACER_H_

#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

void trace(char c, std::string const &s);
class Tracer{
public:
	std::fstream f;
	bool compare;
	int numberOfIgnoredAs;
	int numberOfTrackedAs;
	Tracer(bool compare_, int _numberOfIgnoreAs=0): //alternative to compare is record
		compare(compare_),
		numberOfIgnoredAs(_numberOfIgnoreAs),
		numberOfTrackedAs(0),
		f("TRACE_",std::fstream::binary|(compare_?std::ios_base::in:std::ios_base::out))
	{
	}
	~Tracer()
	{
		trace('Q',"");
	}
};

static Tracer theTracer(1);

void trace(char c, std::string const &s)
{
	if(c=='A')theTracer.numberOfTrackedAs++;
	if(c=='A' && theTracer.numberOfIgnoredAs-->0)return;
	if(!theTracer.compare)
	{
		int l=s.size();
		theTracer.f<<c<<l<<"\n";
		for(int i=0;i<l;i++)theTracer.f<<s[i];
	}
	else
	{
		char C;
		theTracer.f>>C;
		int l;
		theTracer.f>>l;
		std::string S(l,' ');

		char T;theTracer.f.read(&T,1);
		theTracer.f.read(&(S[0]),l);

		if(c!=C)std::cerr<<"Chars do not match\n"<<c<<C;
		if(S.length()!=s.length())std::cerr<<"Lengths do not match\n";
		if(S!=s)std::cerr<<"--------------------------------"<<c<<"\n"<<s<<"--------------------------------"<<C<<"\n"<<S<<"--------------------------------\n";
		if(C!=c || S!=s)std::cerr<<"NumerOfTrackedAs:"<<theTracer.numberOfTrackedAs<<"\n";
		assert(C==c && S==s);
	}
}


#endif /* SRC_TRACER_H_ */
