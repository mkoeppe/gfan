/*
 * gfanlibglue.h
 *
 *  Created on: Jun 27, 2023
 *      Author: anders
 */

#ifndef SRC_GFANLIBGLUE_H_
#define SRC_GFANLIBGLUE_H_


#include <assert.h>

#include "gfanlib.h"
#include "vektor.h"

namespace gfan{
template<typename typ>
gfan::Matrix<typ> convertMatrix(IntMatrix const &m)
{
	gfan::Matrix<typ> ret(m.getHeight(),m.getWidth());
	for(int i=0;i<m.getHeight();i++)
		for(int j=0;j<m.getWidth();j++)
			ret[i][j]=typ(m[i][j]);
	return ret;
}

IntMatrix rowsToIntegerMatrix2(IntegerVectorList const &l);//taken from app_mixedvolume.cpp

template<typename typ>
IntegerVector toIntegerVector(gfan::Vector<typ> const &w)
{
	IntegerVector W(w.size());
	for(int i=0;i<W.size();i++)
		W[i]=w[i].toInt64();
	return W;
}


};

#endif /* SRC_GFANLIBGLUE_H_ */
