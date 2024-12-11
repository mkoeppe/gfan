/*
 * gfanlib_memoryresource.cpp
 *
 *  Created on: 9 Jul 2018
 *      Author: anders
 */

#include "gfanlib_memoryresource.h"

StackResource stackResource(10160000);
StackResource stackResource2(10160000);

std::string memoryResourceToString(const std::experimental::pmr::memory_resource *mr)
{
	if(get_default_resource()==mr)
	{
		return std::string("DEFAULT");
	}
	auto p=dynamic_cast<const StackResource*>(mr);
	if(p)
	{
		return std::string("STACKRESOURCE");//+std::to_string((mr));
	}
	return std::string("UNKNOWN");
}
