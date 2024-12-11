/*
 * gfanlibglue.cpp
 *
 *  Created on: Sep 22, 2023
 *      Author: au27818
 */

#include "gfanlibglue.h"

namespace gfan{

IntMatrix rowsToIntegerMatrix2(IntegerVectorList const &l)//taken from app_mixedvolume.cpp
  {
	  assert(l.size());
	  IntMatrix ret(l.size(),l.front().size());
	  IntegerVectorList::const_iterator I=l.begin();
	  for(int i=0;i<ret.getHeight();i++,I++)
		  for(int j=0;j<ret.getWidth();j++)
			  ret[i][j]=(*I)[j];
	  return ret;
  }

};
