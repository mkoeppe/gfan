/*
 * gfanlib_histogram.h
 *
 *  Created on: Sep 29, 2022
 *      Author: anders
 */

#pragma once

#include <cstdint>
#include <mutex>
#include <vector>
#include <string>
#include <numeric>
#include <iostream>
#include <iomanip>

namespace gfan{
	class FrequencyTable
	{
		std::vector<uint64_t> data;
		std::mutex m;
		std::string name;
		bool report;
	public:
		FrequencyTable(std::string name_,bool report_=false):
			report(report_),
			name(name_)
		{
		}
		void record(unsigned int v)
		{
			if(report)
			{
				std::lock_guard<std::mutex> lck{m};
				if(v>=data.size())
					data.resize(v+1);
				data[v]++;
			}
		}
		~FrequencyTable()
		{
			if(report)
			{
				std::cerr<<"FrequencyTable: "<<name<<"\n";
				auto s=std::accumulate(data.begin(),data.end(),(uint64_t)0);
				std::cerr<<"Total:"<<s<<std::endl;
				std::cerr<<"----------------------------"<<std::endl;
				for(unsigned int i=0;i<data.size();i++)
					if(data[i])
						std::cerr
						<<std::setw(8)<<i
						<<std::setw(12)<<data[i]
						<<std::setw(8)<<(100*data[i]/s)<<std::endl;
				std::cerr<<"----------------------------"<<std::endl;
			}
		}
	};

};


