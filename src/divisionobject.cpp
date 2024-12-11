/*
 * divisionobject.cpp
 *
 *  Created on: Sep 22, 2023
 *      Author: au27818
 */

#include <queue>
#include <set>
#include <memory>
//#include <ranges>
//#include <range/v3/all.hpp>
//#include "log.h"
#include "polynomial.h"
#include "packedmonomial.h"
#include "wallideal.h"
#include "division.h"
#include "gfanlib_frequencytable.h"
#include "gfanlibglue.h"
#include "gfanlib_tableau.h"

#include "divisionobject.h"
#include "field_zmodpz.h"

//namespace gfan{
	std::unique_ptr<class DivisionObject> DivisionObjectFactory(Field const &k, PolynomialSet &g_, IntegerVector const *wHint, TermOrder const &termOrder)
	{
		// This line disables any intelligent DivisionObject:
		return std::make_unique<DivisionObjectNaive>(g_,termOrder);

		if(auto *p=dynamic_cast<FieldZModPZImplementation*>(k.implementingObject))
		{
//			std::cerr<<"Factory of kind 1\n";
//			return std::unique_ptr<DivisionObject>((DivisionObject*)new DivisionObjectFast<4,FieldZModPZImplementation>(g_,wHint,termOrder));
			return std::make_unique<DivisionObjectFast<4,FieldZModPZImplementation>>(g_,wHint,termOrder);
		}
//		return std::unique_ptr<DivisionObject>(new DivisionObjectFast<4,void>(g_,wHint,termOrder));
//		std::cerr<<"Factory of kind 2\n"<<k.toString().c_str()<<"\n";
		return std::make_unique<DivisionObjectFast<4,void>>(g_,wHint,termOrder);
	}

	void DivisionObjectNaive::reduce(Polynomial &f)
	{
		f=division(f,g,termOrder);
	}

//};
